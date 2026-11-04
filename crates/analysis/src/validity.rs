//! Geometry validity checking and repair following the OGC Simple Features Specification (SFS).
//!
//! Provides two geoprocessing capabilities:
//! - **Check Geometry**: Identifies validity issues (self-intersections, unclosed rings,
//!   degenerate rings, holes outside shells, overlapping holes, component overlaps, non-finite
//!   coordinates) and outputs a detailed validation report table.
//! - **Fix Geometries**: Automatically repairs invalid geometries using constrained Delaunay
//!   triangulation ([Hugo Ledoux et al., 2014](http://dx.doi.org/10.1016/j.cageo.2014.01.009)),
//!   closing rings, decomposing self-intersecting polygons into valid components, and
//!   re-nesting holes, while preserving all feature attributes and layer CRS.

use std::collections::HashMap;

use geo::validation::{InvalidMultiPolygon, InvalidPolygon, RingRole};
use geo::{Area, MakeValid, Validation};
use rustgis_data::{Coord, Dataset, Field, FieldType, Geometry, GeometryType, Value};
use serde_json::Value as Json;

use crate::topology::{CrossingCheck, check_self_crossing, segments_cross};
use crate::{Inputs, Output, ToolResult};

/// Status of a feature's geometry validation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValidityStatus {
    Valid,
    Invalid,
    Empty,
    Null,
    Unsupported,
    Error,
}

impl ValidityStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            ValidityStatus::Valid => "Valid",
            ValidityStatus::Invalid => "Invalid",
            ValidityStatus::Empty => "Empty",
            ValidityStatus::Null => "Null",
            ValidityStatus::Unsupported => "Unsupported",
            ValidityStatus::Error => "Error",
        }
    }
}

/// Category of validity violation or geometry state.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ValidityCategory {
    None,
    SelfIntersection,
    UnclosedRing,
    TooFewPoints,
    NonFiniteCoord,
    HoleOutsideShell,
    HolesOverlap,
    ComponentsOverlap,
    EmptyGeometry,
    MissingGeometry,
    UnsupportedType,
    ProcessingError,
}

impl ValidityCategory {
    pub fn as_str(self) -> &'static str {
        match self {
            ValidityCategory::None => "None",
            ValidityCategory::SelfIntersection => "Self-Intersection",
            ValidityCategory::UnclosedRing => "Unclosed Ring",
            ValidityCategory::TooFewPoints => "Too Few Points",
            ValidityCategory::NonFiniteCoord => "Non-Finite Coordinate",
            ValidityCategory::HoleOutsideShell => "Hole Outside Shell",
            ValidityCategory::HolesOverlap => "Holes Overlap",
            ValidityCategory::ComponentsOverlap => "Components Overlap",
            ValidityCategory::EmptyGeometry => "Empty Geometry",
            ValidityCategory::MissingGeometry => "Missing Geometry",
            ValidityCategory::UnsupportedType => "Unsupported Type",
            ValidityCategory::ProcessingError => "Processing Error",
        }
    }
}

/// Detailed record from validating a single feature.
#[derive(Clone, Debug, PartialEq)]
pub struct ValidationRecord {
    pub feature_id: u64,
    pub geometry_type: &'static str,
    pub status: ValidityStatus,
    pub valid: bool,
    pub error_category: ValidityCategory,
    pub error_reason: String,
    pub error_location: Option<Coord>,
}

/// Detailed counts of dropped or degenerate elements during geometry repair.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FixStats {
    pub dropped_non_finite: usize,
    pub dropped_parts: usize,
    pub dropped_holes: usize,
    pub dropped_components: usize,
}

/// The outcome of attempting to repair a geometry.
#[derive(Clone, Debug, PartialEq)]
pub enum FixOutcome {
    /// The geometry was already valid; preserved with identical coordinates.
    Unchanged(Geometry),
    /// The geometry was invalid and successfully repaired.
    Repaired(Geometry),
    /// The geometry could not be repaired (e.g. collapsed to zero area or non-finite).
    Failed(String),
    /// The geometry was empty or null; skipped.
    Skipped(&'static str),
}

fn to_geo_coord(c: Coord) -> geo::Coord<f64> {
    geo::Coord { x: c.x, y: c.y }
}

fn from_geo_coord(c: geo::Coord<f64>) -> Coord {
    Coord { x: c.x, y: c.y }
}

fn ring_role_label(role: RingRole) -> String {
    match role {
        RingRole::Exterior => "exterior ring".to_string(),
        RingRole::Interior(idx) => format!("interior ring {idx}"),
    }
}

fn ring_role_index(role: RingRole) -> usize {
    match role {
        RingRole::Exterior => 0,
        RingRole::Interior(idx) => idx + 1,
    }
}

fn find_rings_intersection(ring1: &[Coord], ring2: &[Coord]) -> Option<Coord> {
    let segs1: Vec<(Coord, Coord)> = ring1.windows(2).filter_map(|w| Some((*w.first()?, *w.get(1)?))).collect();
    let segs2: Vec<(Coord, Coord)> = ring2.windows(2).filter_map(|w| Some((*w.first()?, *w.get(1)?))).collect();
    for a in &segs1 {
        for b in &segs2 {
            if let Some(c) = segments_cross(*a, *b) {
                return Some(c);
            }
        }
    }
    ring2.first().copied().or_else(|| ring1.first().copied())
}

/// Converts RustGIS polygon rings into a `geo::MultiPolygon<f64>` in native coordinates.
fn to_geo_multipolygon(polys: &[Vec<Vec<Coord>>]) -> Option<geo::MultiPolygon<f64>> {
    let mut geo_polys = Vec::new();
    for rings in polys {
        let mut it = rings.iter();
        let Some(exterior_coords) = it.next() else { continue };
        if exterior_coords.is_empty() {
            continue;
        }
        let exterior = geo::LineString::new(exterior_coords.iter().copied().map(to_geo_coord).collect());
        let interiors = it.map(|hole| geo::LineString::new(hole.iter().copied().map(to_geo_coord).collect())).collect();
        geo_polys.push(geo::Polygon::new(exterior, interiors));
    }
    if geo_polys.is_empty() {
        return None;
    }
    Some(geo::MultiPolygon::new(geo_polys))
}

/// Converts a `geo::MultiPolygon<f64>` back to RustGIS `Geometry::Polygon` in native coordinates.
fn from_geo_multipolygon(mp: &geo::MultiPolygon<f64>) -> Geometry {
    let polys: Vec<Vec<Vec<Coord>>> =
        mp.0.iter()
            .map(|p| {
                let ext = p.exterior().0.iter().copied().map(from_geo_coord).collect();
                let mut rings = vec![ext];
                for hole in p.interiors() {
                    rings.push(hole.0.iter().copied().map(from_geo_coord).collect());
                }
                rings
            })
            .collect();
    Geometry::Polygon(polys)
}

/// Validates a single feature against OGC SFS validity rules.
pub fn validate_geometry(geom: Option<&Geometry>, feature_id: u64, layer_kind: GeometryType) -> ValidationRecord {
    let Some(g) = geom else {
        return ValidationRecord {
            feature_id,
            geometry_type: layer_kind.label(),
            status: ValidityStatus::Null,
            valid: false,
            error_category: ValidityCategory::MissingGeometry,
            error_reason: "Feature has missing (null) geometry".into(),
            error_location: None,
        };
    };

    match g {
        Geometry::Point(c) => validate_point(*c, feature_id),
        Geometry::MultiPoint(pts) => validate_multipoint(pts, feature_id),
        Geometry::Line(parts) => validate_line(parts, feature_id),
        Geometry::Polygon(polys) => validate_polygon(polys, feature_id),
    }
}

fn validate_point(c: Coord, feature_id: u64) -> ValidationRecord {
    if !c.is_finite() {
        return ValidationRecord {
            feature_id,
            geometry_type: "Point",
            status: ValidityStatus::Invalid,
            valid: false,
            error_category: ValidityCategory::NonFiniteCoord,
            error_reason: "Point contains non-finite coordinate (NaN or Inf)".into(),
            error_location: None,
        };
    }
    ValidationRecord {
        feature_id,
        geometry_type: "Point",
        status: ValidityStatus::Valid,
        valid: true,
        error_category: ValidityCategory::None,
        error_reason: "Geometry is valid".into(),
        error_location: None,
    }
}

fn validate_multipoint(pts: &[Coord], feature_id: u64) -> ValidationRecord {
    if pts.is_empty() {
        return ValidationRecord {
            feature_id,
            geometry_type: "MultiPoint",
            status: ValidityStatus::Empty,
            valid: false,
            error_category: ValidityCategory::EmptyGeometry,
            error_reason: "MultiPoint has no coordinates".into(),
            error_location: None,
        };
    }
    for (i, c) in pts.iter().enumerate() {
        if !c.is_finite() {
            return ValidationRecord {
                feature_id,
                geometry_type: "MultiPoint",
                status: ValidityStatus::Invalid,
                valid: false,
                error_category: ValidityCategory::NonFiniteCoord,
                error_reason: format!("MultiPoint coordinate at index {i} is non-finite (NaN or Inf)"),
                error_location: None,
            };
        }
    }
    ValidationRecord {
        feature_id,
        geometry_type: "MultiPoint",
        status: ValidityStatus::Valid,
        valid: true,
        error_category: ValidityCategory::None,
        error_reason: "Geometry is valid".into(),
        error_location: None,
    }
}

fn validate_line(parts: &[Vec<Coord>], feature_id: u64) -> ValidationRecord {
    if parts.is_empty() || parts.iter().all(|p| p.is_empty()) {
        return ValidationRecord {
            feature_id,
            geometry_type: "Polyline",
            status: ValidityStatus::Empty,
            valid: false,
            error_category: ValidityCategory::EmptyGeometry,
            error_reason: "Polyline has no vertices or parts".into(),
            error_location: None,
        };
    }

    for (part_idx, part) in parts.iter().enumerate() {
        for (coord_idx, c) in part.iter().enumerate() {
            if !c.is_finite() {
                return ValidationRecord {
                    feature_id,
                    geometry_type: "Polyline",
                    status: ValidityStatus::Invalid,
                    valid: false,
                    error_category: ValidityCategory::NonFiniteCoord,
                    error_reason: format!("Polyline part {part_idx} coordinate {coord_idx} is non-finite"),
                    error_location: None,
                };
            }
        }

        let mut distinct = part.clone();
        distinct.dedup();
        if distinct.len() < 2 {
            return ValidationRecord {
                feature_id,
                geometry_type: "Polyline",
                status: ValidityStatus::Invalid,
                valid: false,
                error_category: ValidityCategory::TooFewPoints,
                error_reason: format!("Polyline part {part_idx} has fewer than 2 distinct points"),
                error_location: part.first().copied(),
            };
        }

        match check_self_crossing(part, false) {
            CrossingCheck::Crossing(cross_loc) => {
                return ValidationRecord {
                    feature_id,
                    geometry_type: "Polyline",
                    status: ValidityStatus::Invalid,
                    valid: false,
                    error_category: ValidityCategory::SelfIntersection,
                    error_reason: format!("Polyline part {part_idx} self-intersects at ({:.6}, {:.6})", cross_loc.x, cross_loc.y),
                    error_location: Some(cross_loc),
                };
            }
            CrossingCheck::BudgetExceeded => {
                return ValidationRecord {
                    feature_id,
                    geometry_type: "Polyline",
                    status: ValidityStatus::Error,
                    valid: false,
                    error_category: ValidityCategory::ProcessingError,
                    error_reason: format!("Polyline part {part_idx}: work budget exceeded during self-intersection check"),
                    error_location: part.first().copied(),
                };
            }
            CrossingCheck::None => {}
        }
    }

    ValidationRecord {
        feature_id,
        geometry_type: "Polyline",
        status: ValidityStatus::Valid,
        valid: true,
        error_category: ValidityCategory::None,
        error_reason: "Geometry is valid".into(),
        error_location: None,
    }
}

fn validate_polygon(polys: &[Vec<Vec<Coord>>], feature_id: u64) -> ValidationRecord {
    if polys.is_empty() || polys.iter().all(|rings| rings.iter().all(|ring| ring.is_empty())) {
        return ValidationRecord {
            feature_id,
            geometry_type: "Polygon",
            status: ValidityStatus::Empty,
            valid: false,
            error_category: ValidityCategory::EmptyGeometry,
            error_reason: "Polygon has no rings or coordinates".into(),
            error_location: None,
        };
    }

    // Step 1: Check ring-level properties (finiteness, closure, point counts, self-crossings)
    for (poly_idx, rings) in polys.iter().enumerate() {
        for (ring_idx, ring) in rings.iter().enumerate() {
            let role = if ring_idx == 0 { "exterior ring" } else { "interior ring (hole)" };

            // 1a. Finiteness
            for (coord_idx, c) in ring.iter().enumerate() {
                if !c.is_finite() {
                    return ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::NonFiniteCoord,
                        error_reason: format!("{role} in polygon {poly_idx} coordinate {coord_idx} is non-finite"),
                        error_location: None,
                    };
                }
            }

            // 1b. Closure
            if ring.len() >= 2 && ring.first() != ring.last() {
                return ValidationRecord {
                    feature_id,
                    geometry_type: "Polygon",
                    status: ValidityStatus::Invalid,
                    valid: false,
                    error_category: ValidityCategory::UnclosedRing,
                    error_reason: format!("{role} in polygon {poly_idx} is not closed: first vertex != last vertex"),
                    error_location: ring.last().copied(),
                };
            }

            // 1c. Minimum vertices: at least 4 coordinates (3 distinct + closure)
            let mut distinct = ring.clone();
            distinct.dedup();
            if distinct.len() >= 2 && distinct.first() == distinct.last() {
                distinct.pop();
            }
            if ring.len() < 4 || distinct.len() < 3 {
                return ValidationRecord {
                    feature_id,
                    geometry_type: "Polygon",
                    status: ValidityStatus::Invalid,
                    valid: false,
                    error_category: ValidityCategory::TooFewPoints,
                    error_reason: format!(
                        "{role} in polygon {poly_idx} has too few distinct points (requires >= 3, found {})",
                        distinct.len()
                    ),
                    error_location: ring.first().copied(),
                };
            }

            // 1d. Ring self-crossing
            match check_self_crossing(ring, true) {
                CrossingCheck::Crossing(cross_loc) => {
                    return ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::SelfIntersection,
                        error_reason: format!("{role} in polygon {poly_idx} self-intersects at ({:.6}, {:.6})", cross_loc.x, cross_loc.y),
                        error_location: Some(cross_loc),
                    };
                }
                CrossingCheck::BudgetExceeded => {
                    return ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Error,
                        valid: false,
                        error_category: ValidityCategory::ProcessingError,
                        error_reason: format!("{role} in polygon {poly_idx}: work budget exceeded during self-intersection check"),
                        error_location: ring.first().copied(),
                    };
                }
                CrossingCheck::None => {}
            }
        }
    }

    // Step 2: Use geo::Validation for containment and topological relationships
    let Some(geo_mp) = to_geo_multipolygon(polys) else {
        return ValidationRecord {
            feature_id,
            geometry_type: "Polygon",
            status: ValidityStatus::Error,
            valid: false,
            error_category: ValidityCategory::ProcessingError,
            error_reason: "Failed to construct geometry representation for validation".into(),
            error_location: None,
        };
    };

    if let Err(err) = geo_mp.check_validation() {
        return interpret_geo_multipolygon_error(err, polys, feature_id);
    }

    ValidationRecord {
        feature_id,
        geometry_type: "Polygon",
        status: ValidityStatus::Valid,
        valid: true,
        error_category: ValidityCategory::None,
        error_reason: "Geometry is valid".into(),
        error_location: None,
    }
}

fn interpret_geo_multipolygon_error(err: InvalidMultiPolygon, polys: &[Vec<Vec<Coord>>], feature_id: u64) -> ValidationRecord {
    match err {
        InvalidMultiPolygon::InvalidPolygon(geom_idx, poly_err) => {
            let p_idx = geom_idx.0;
            match poly_err {
                InvalidPolygon::SelfIntersection(role) => {
                    let role_str = match role {
                        RingRole::Exterior => "exterior ring",
                        RingRole::Interior(_) => "interior ring (hole)",
                    };
                    let loc = polys.get(p_idx).and_then(|rings| {
                        let r_idx = match role {
                            RingRole::Exterior => 0,
                            RingRole::Interior(i) => i + 1,
                        };
                        rings.get(r_idx).and_then(|r| match check_self_crossing(r, true) {
                            CrossingCheck::Crossing(c) => Some(c),
                            _ => r.first().copied(),
                        })
                    });
                    ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::SelfIntersection,
                        error_reason: format!("{role_str} in polygon {p_idx} self-intersects"),
                        error_location: loc,
                    }
                }
                InvalidPolygon::TooFewPointsInRing(role) => {
                    let role_str = match role {
                        RingRole::Exterior => "exterior ring",
                        RingRole::Interior(_) => "interior ring (hole)",
                    };
                    let loc = polys.get(p_idx).and_then(|rings| {
                        let r_idx = match role {
                            RingRole::Exterior => 0,
                            RingRole::Interior(i) => i + 1,
                        };
                        rings.get(r_idx).and_then(|r| r.first().copied())
                    });
                    ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::TooFewPoints,
                        error_reason: format!("{role_str} in polygon {p_idx} has too few points"),
                        error_location: loc,
                    }
                }
                InvalidPolygon::NonFiniteCoord(_, _) => ValidationRecord {
                    feature_id,
                    geometry_type: "Polygon",
                    status: ValidityStatus::Invalid,
                    valid: false,
                    error_category: ValidityCategory::NonFiniteCoord,
                    error_reason: format!("polygon {p_idx} contains non-finite coordinates"),
                    error_location: None,
                },
                InvalidPolygon::InteriorRingNotContainedInExteriorRing(RingRole::Interior(hole_idx)) => {
                    let loc = polys.get(p_idx).and_then(|rings| rings.get(hole_idx + 1)).and_then(|h| h.first().copied());
                    ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::HoleOutsideShell,
                        error_reason: format!("interior ring {hole_idx} in polygon {p_idx} is not contained inside exterior ring"),
                        error_location: loc,
                    }
                }
                InvalidPolygon::InteriorRingNotContainedInExteriorRing(_) => ValidationRecord {
                    feature_id,
                    geometry_type: "Polygon",
                    status: ValidityStatus::Invalid,
                    valid: false,
                    error_category: ValidityCategory::HoleOutsideShell,
                    error_reason: format!("interior ring in polygon {p_idx} is not contained inside exterior ring"),
                    error_location: None,
                },
                InvalidPolygon::IntersectingRingsOnAnArea(r1, r2) => {
                    let r1_label = ring_role_label(r1);
                    let r2_label = ring_role_label(r2);
                    let r1_idx = ring_role_index(r1);
                    let r2_idx = ring_role_index(r2);
                    let loc = polys.get(p_idx).and_then(|rings| {
                        let ring1 = rings.get(r1_idx)?;
                        let ring2 = rings.get(r2_idx)?;
                        find_rings_intersection(ring1, ring2)
                    });
                    ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::HolesOverlap,
                        error_reason: format!("{r1_label} and {r2_label} in polygon {p_idx} overlap on an area"),
                        error_location: loc,
                    }
                }
                InvalidPolygon::IntersectingRingsOnALine(r1, r2) => {
                    let r1_label = ring_role_label(r1);
                    let r2_label = ring_role_label(r2);
                    let r1_idx = ring_role_index(r1);
                    let r2_idx = ring_role_index(r2);
                    let loc = polys.get(p_idx).and_then(|rings| {
                        let ring1 = rings.get(r1_idx)?;
                        let ring2 = rings.get(r2_idx)?;
                        find_rings_intersection(ring1, ring2)
                    });
                    ValidationRecord {
                        feature_id,
                        geometry_type: "Polygon",
                        status: ValidityStatus::Invalid,
                        valid: false,
                        error_category: ValidityCategory::HolesOverlap,
                        error_reason: format!("rings {r1_label} and {r2_label} in polygon {p_idx} touch or intersect on a line"),
                        error_location: loc,
                    }
                }
            }
        }
        InvalidMultiPolygon::ElementsOverlaps(idx1, idx2) => {
            let loc = polys
                .get(idx1.0)
                .and_then(|p1| {
                    polys.get(idx2.0).and_then(|p2| {
                        let r1 = p1.first()?;
                        let r2 = p2.first()?;
                        find_rings_intersection(r1, r2)
                    })
                })
                .or_else(|| polys.get(idx2.0).and_then(|rings| rings.first()).and_then(|r| r.first().copied()));
            ValidationRecord {
                feature_id,
                geometry_type: "Polygon",
                status: ValidityStatus::Invalid,
                valid: false,
                error_category: ValidityCategory::ComponentsOverlap,
                error_reason: format!("polygons at indices {} and {} in MultiPolygon overlap on an area", idx1.0, idx2.0),
                error_location: loc,
            }
        }
        InvalidMultiPolygon::ElementsTouchOnALine(idx1, idx2) => {
            let loc = polys
                .get(idx1.0)
                .and_then(|p1| {
                    polys.get(idx2.0).and_then(|p2| {
                        let r1 = p1.first()?;
                        let r2 = p2.first()?;
                        find_rings_intersection(r1, r2)
                    })
                })
                .or_else(|| polys.get(idx2.0).and_then(|rings| rings.first()).and_then(|r| r.first().copied()));
            ValidationRecord {
                feature_id,
                geometry_type: "Polygon",
                status: ValidityStatus::Invalid,
                valid: false,
                error_category: ValidityCategory::ComponentsOverlap,
                error_reason: format!("polygons at indices {} and {} in MultiPolygon touch along a line", idx1.0, idx2.0),
                error_location: loc,
            }
        }
    }
}

/// Automatically repairs an invalid geometry while preserving its spatial meaning.
pub fn fix_geometry(geom: &Geometry) -> FixOutcome {
    let (outcome, _) = fix_geometry_with_stats(geom);
    outcome
}

/// Automatically repairs an invalid geometry and returns repair statistics.
pub fn fix_geometry_with_stats(geom: &Geometry) -> (FixOutcome, FixStats) {
    match geom {
        Geometry::Point(c) => (fix_point(*c), FixStats::default()),
        Geometry::MultiPoint(pts) => fix_multipoint_with_stats(pts),
        Geometry::Line(parts) => fix_line_with_stats(parts),
        Geometry::Polygon(polys) => fix_polygon_with_stats(polys),
    }
}

fn fix_point(c: Coord) -> FixOutcome {
    if c.is_finite() {
        FixOutcome::Unchanged(Geometry::Point(c))
    } else {
        FixOutcome::Failed("point coordinates are non-finite (NaN or Inf)".into())
    }
}

fn fix_multipoint_with_stats(pts: &[Coord]) -> (FixOutcome, FixStats) {
    let mut stats = FixStats::default();
    if pts.is_empty() {
        return (FixOutcome::Skipped("MultiPoint has no coordinates"), stats);
    }
    let clean: Vec<Coord> = pts.iter().copied().filter(|c| c.is_finite()).collect();
    let dropped = pts.len() - clean.len();
    stats.dropped_non_finite = dropped;
    if clean.is_empty() {
        return (FixOutcome::Failed("all points were non-finite (NaN or Inf)".into()), stats);
    }
    if dropped == 0 {
        (FixOutcome::Unchanged(Geometry::MultiPoint(pts.to_vec())), stats)
    } else {
        (FixOutcome::Repaired(Geometry::MultiPoint(clean)), stats)
    }
}

fn fix_line_with_stats(parts: &[Vec<Coord>]) -> (FixOutcome, FixStats) {
    let mut stats = FixStats::default();
    if parts.is_empty() || parts.iter().all(|p| p.is_empty()) {
        return (FixOutcome::Skipped("polyline has no vertices or parts"), stats);
    }

    let mut clean_parts: Vec<Vec<Coord>> = Vec::new();
    let mut modified = false;

    for (part_idx, part) in parts.iter().enumerate() {
        let mut clean: Vec<Coord> = Vec::with_capacity(part.len());
        let mut non_finite_in_part = 0usize;
        for c in part {
            if c.is_finite() {
                clean.push(*c);
            } else {
                non_finite_in_part += 1;
            }
        }
        if non_finite_in_part > 0 {
            stats.dropped_non_finite += non_finite_in_part;
            modified = true;
        }

        let orig_clean_len = clean.len();
        clean.dedup();
        if clean.len() != orig_clean_len {
            modified = true;
        }

        if clean.len() >= 2 {
            match check_self_crossing(&clean, false) {
                CrossingCheck::Crossing(c) => {
                    return (
                        FixOutcome::Failed(format!(
                            "polyline part {part_idx} self-intersects at ({:.6}, {:.6}) and cannot be automatically repaired",
                            c.x, c.y
                        )),
                        stats,
                    );
                }
                CrossingCheck::BudgetExceeded => {
                    return (FixOutcome::Failed("work budget exceeded during self-intersection check".into()), stats);
                }
                CrossingCheck::None => {}
            }
            clean_parts.push(clean);
        } else {
            stats.dropped_parts += 1;
            modified = true;
        }
    }

    if clean_parts.is_empty() {
        return (FixOutcome::Failed("all line parts collapsed to degenerate (fewer than 2 points)".into()), stats);
    }

    if !modified && clean_parts.len() == parts.len() {
        (FixOutcome::Unchanged(Geometry::Line(parts.to_vec())), stats)
    } else {
        (FixOutcome::Repaired(Geometry::Line(clean_parts)), stats)
    }
}

fn clean_ring_coords(ring: &[Coord]) -> (Vec<Coord>, usize) {
    let mut clean: Vec<Coord> = Vec::with_capacity(ring.len());
    let mut dropped_non_finite = 0usize;
    for c in ring {
        if c.is_finite() {
            clean.push(*c);
        } else {
            dropped_non_finite += 1;
        }
    }
    clean.dedup();
    if let Some(first) = clean.first().copied()
        && clean.last() != Some(&first)
    {
        clean.push(first);
    }
    (clean, dropped_non_finite)
}

fn fix_polygon_with_stats(polys: &[Vec<Vec<Coord>>]) -> (FixOutcome, FixStats) {
    let mut stats = FixStats::default();
    if polys.is_empty() || polys.iter().all(|rings| rings.iter().all(|ring| ring.is_empty())) {
        return (FixOutcome::Skipped("polygon has no rings or coordinates"), stats);
    }

    // Fast path: if already valid, keep 100% bitwise coordinates intact.
    let validation = validate_polygon(polys, 0);
    if validation.valid {
        return (FixOutcome::Unchanged(Geometry::Polygon(polys.to_vec())), stats);
    }
    if validation.status == ValidityStatus::Error {
        return (FixOutcome::Failed("work budget exceeded during validation; cannot determine validity safely".into()), stats);
    }

    // Step 1: Pre-clean rings (filter non-finite, dedup consecutive duplicate vertices, ensure closed)
    let mut cleaned_polys: Vec<geo::Polygon<f64>> = Vec::new();
    for rings in polys {
        let mut ring_iter = rings.iter();
        let Some(exterior_ring) = ring_iter.next() else {
            stats.dropped_components += 1;
            continue;
        };

        let (clean_ext, non_finite_ext) = clean_ring_coords(exterior_ring);
        stats.dropped_non_finite += non_finite_ext;

        let mut distinct_ext = clean_ext.clone();
        distinct_ext.dedup();
        if distinct_ext.len() >= 2 && distinct_ext.first() == distinct_ext.last() {
            distinct_ext.pop();
        }

        if clean_ext.len() < 4 || distinct_ext.len() < 3 {
            // Exterior ring collapsed: do not promote a hole to the exterior shell!
            stats.dropped_components += 1;
            continue;
        }

        let mut cleaned_interiors = Vec::new();
        for hole in ring_iter {
            let (clean_hole, non_finite_hole) = clean_ring_coords(hole);
            stats.dropped_non_finite += non_finite_hole;

            let mut distinct_hole = clean_hole.clone();
            distinct_hole.dedup();
            if distinct_hole.len() >= 2 && distinct_hole.first() == distinct_hole.last() {
                distinct_hole.pop();
            }

            if clean_hole.len() >= 4 && distinct_hole.len() >= 3 {
                cleaned_interiors.push(geo::LineString::new(clean_hole.into_iter().map(to_geo_coord).collect()));
            } else {
                stats.dropped_holes += 1;
            }
        }

        let exterior_ls = geo::LineString::new(clean_ext.into_iter().map(to_geo_coord).collect());
        cleaned_polys.push(geo::Polygon::new(exterior_ls, cleaned_interiors));
    }

    if cleaned_polys.is_empty() {
        return (FixOutcome::Failed("all polygon components collapsed to degenerate".into()), stats);
    }

    let input_mp = geo::MultiPolygon::new(cleaned_polys);

    // Step 2: Apply triangulation-based polygon repair (Hugo Ledoux et al., 2014)
    match input_mp.make_valid() {
        Ok(repaired_mp) => {
            if repaired_mp.0.is_empty() {
                return (FixOutcome::Failed("polygon collapsed to empty during repair".into()), stats);
            }
            let area: f64 = repaired_mp.0.iter().map(|p| p.unsigned_area()).sum();
            if area <= 1e-14 {
                return (FixOutcome::Failed("repaired polygon has zero area".into()), stats);
            }

            let new_geom = from_geo_multipolygon(&repaired_mp);
            (FixOutcome::Repaired(new_geom), stats)
        }
        Err(err) => (FixOutcome::Failed(format!("repair failed: {err}")), stats),
    }
}

fn s<'a>(p: &'a Json, k: &str) -> &'a str {
    p.get(k).and_then(Json::as_str).unwrap_or("")
}

fn out_name(p: &Json, default: String) -> String {
    let n = s(p, "output").trim();
    if n.is_empty() { default } else { n.to_string() }
}

fn input<'a>(inp: &'a dyn Inputs, p: &Json, k: &str) -> Result<&'a Dataset, String> {
    if let Some(n) = p.get(k).filter(|v| v.is_number()) {
        let table = n.as_u64().and_then(|id| u32::try_from(id).ok()).and_then(|id| inp.table(id));
        return table.ok_or_else(|| format!("no standalone table with id {n}"));
    }
    let name = s(p, k);
    inp.dataset(name).ok_or_else(|| format!("no layer named '{name}'"))
}

fn rows(inp: &dyn Inputs, p: &Json, k: &str, ds: &Dataset) -> Vec<usize> {
    let sel = inp.selection(s(p, k));
    let len = ds.features.len();
    if sel.is_empty() { (0..len).collect() } else { sel.iter().copied().filter(|r| *r < len).collect() }
}

/// The runner function for the `checkGeometry` geoprocessing tool.
pub fn check_geometry(inp: &dyn Inputs, p: &Json) -> Result<ToolResult, String> {
    let ds = input(inp, p, "input")?;
    let target_rows = rows(inp, p, "input", ds);

    let fields = vec![
        Field::new("feature_id", FieldType::Integer),
        Field::new("geometry_type", FieldType::Text),
        Field::new("status", FieldType::Text),
        Field::new("valid", FieldType::Boolean),
        Field::new("error_category", FieldType::Text),
        Field::new("error_reason", FieldType::Text),
        Field::new("error_x", FieldType::Double),
        Field::new("error_y", FieldType::Double),
    ];

    let output_name = out_name(p, format!("{}_Validation", ds.name));
    let mut out_table = Dataset::new(output_name, GeometryType::None, fields);

    let total = target_rows.len();
    let mut valid_count = 0usize;
    let mut invalid_count = 0usize;
    let mut empty_count = 0usize;
    let mut null_count = 0usize;
    let mut category_counts: HashMap<ValidityCategory, usize> = HashMap::new();

    crate::control::begin(total);
    for (i, &row_idx) in target_rows.iter().enumerate() {
        if crate::control::cancelled() {
            return Err("Check Geometry was cancelled".into());
        }
        crate::control::reached(i);

        let Some(f) = ds.features.get(row_idx) else { continue };
        let record = validate_geometry(f.geometry.as_ref(), f.id, ds.geometry_type);

        match record.status {
            ValidityStatus::Valid => valid_count += 1,
            ValidityStatus::Invalid => {
                invalid_count += 1;
                *category_counts.entry(record.error_category).or_insert(0) += 1;
            }
            ValidityStatus::Empty => empty_count += 1,
            ValidityStatus::Null => null_count += 1,
            ValidityStatus::Unsupported | ValidityStatus::Error => {
                invalid_count += 1;
                *category_counts.entry(record.error_category).or_insert(0) += 1;
            }
        }

        let err_x = record.error_location.map_or(Value::Null, |c| Value::Float(c.x));
        let err_y = record.error_location.map_or(Value::Null, |c| Value::Float(c.y));

        let attrs = vec![
            Value::Int(record.feature_id as i64),
            Value::Text(record.geometry_type.into()),
            Value::Text(record.status.as_str().into()),
            Value::Bool(record.valid),
            Value::Text(record.error_category.as_str().into()),
            Value::Text(record.error_reason),
            err_x,
            err_y,
        ];

        out_table.push_next(None, attrs);
    }

    let mut messages = Vec::new();
    messages.push(format!(
        "Check Geometry complete: {total} feature(s) checked ({valid_count} valid, {invalid_count} invalid, {empty_count} empty, {null_count} null)."
    ));

    if invalid_count > 0 {
        let mut breakdown: Vec<String> = category_counts.iter().map(|(cat, count)| format!("{count} {}", cat.as_str())).collect();
        breakdown.sort();
        messages.push(format!("Identified issues: {}.", breakdown.join(", ")));
    }

    Ok(ToolResult { output: Output::Table(out_table), messages })
}

/// The runner function for the `fixGeometries` geoprocessing tool.
pub fn fix_geometries(inp: &dyn Inputs, p: &Json) -> Result<ToolResult, String> {
    let ds = input(inp, p, "input")?;
    if ds.geometry_type == GeometryType::None {
        return Err(format!("layer '{}' has no geometry (standalone table)", ds.name));
    }
    let target_rows = rows(inp, p, "input", ds);

    let output_name = out_name(p, format!("{}_Fixed", ds.name));
    let mut out_ds = Dataset::new(output_name, ds.geometry_type, ds.fields.clone());
    out_ds.crs = ds.crs;
    out_ds.crs_incomplete = ds.crs_incomplete;
    out_ds.zm_dropped = ds.zm_dropped;

    let total = target_rows.len();
    let mut unchanged_count = 0usize;
    let mut repaired_count = 0usize;
    let mut failed_count = 0usize;
    let mut skipped_count = 0usize;
    let mut failed_details: Vec<String> = Vec::new();

    let mut total_dropped_non_finite = 0usize;
    let mut features_with_dropped_non_finite = 0usize;
    let mut total_dropped_parts = 0usize;
    let mut features_with_dropped_parts = 0usize;
    let mut total_dropped_holes = 0usize;
    let mut features_with_dropped_holes = 0usize;
    let mut total_dropped_components = 0usize;
    let mut features_with_dropped_components = 0usize;

    crate::control::begin(total);
    for (i, &row_idx) in target_rows.iter().enumerate() {
        if crate::control::cancelled() {
            return Err("Fix Geometries was cancelled".into());
        }
        crate::control::reached(i);

        let Some(f) = ds.features.get(row_idx) else { continue };

        let Some(g) = &f.geometry else {
            skipped_count += 1;
            out_ds.push_next(None, f.attrs.clone());
            continue;
        };

        let (outcome, stats) = fix_geometry_with_stats(g);

        if stats.dropped_non_finite > 0 {
            total_dropped_non_finite += stats.dropped_non_finite;
            features_with_dropped_non_finite += 1;
        }
        if stats.dropped_parts > 0 {
            total_dropped_parts += stats.dropped_parts;
            features_with_dropped_parts += 1;
        }
        if stats.dropped_holes > 0 {
            total_dropped_holes += stats.dropped_holes;
            features_with_dropped_holes += 1;
        }
        if stats.dropped_components > 0 {
            total_dropped_components += stats.dropped_components;
            features_with_dropped_components += 1;
        }

        match outcome {
            FixOutcome::Unchanged(orig) => {
                unchanged_count += 1;
                out_ds.push_next(Some(orig), f.attrs.clone());
            }
            FixOutcome::Repaired(repaired) => {
                repaired_count += 1;
                out_ds.push_next(Some(repaired), f.attrs.clone());
            }
            FixOutcome::Failed(reason) => {
                failed_count += 1;
                if failed_details.len() < 5 {
                    failed_details.push(format!("Feature {}: {}", f.id, reason));
                }
                // Preserve feature with original geometry rather than silently dropping
                out_ds.push_next(Some(g.clone()), f.attrs.clone());
            }
            FixOutcome::Skipped(_) => {
                skipped_count += 1;
                out_ds.push_next(None, f.attrs.clone());
            }
        }
    }

    let mut messages = Vec::new();
    messages.push(format!(
        "Fix Geometries complete: {total} feature(s) processed ({unchanged_count} unchanged, {repaired_count} repaired, {failed_count} failed, {skipped_count} skipped)."
    ));

    if failed_count > 0 {
        messages.push(format!("Failed repairs (original geometries preserved): {}", failed_details.join("; ")));
    }

    if total_dropped_non_finite > 0 {
        messages.push(format!(
            "{total_dropped_non_finite} non-finite coordinate(s) dropped across {features_with_dropped_non_finite} feature(s)."
        ));
    }
    if total_dropped_parts > 0 {
        messages.push(format!(
            "{total_dropped_parts} degenerate line part(s) (fewer than 2 points) dropped across {features_with_dropped_parts} feature(s)."
        ));
    }
    if total_dropped_holes > 0 {
        messages.push(format!(
            "{total_dropped_holes} degenerate interior ring(s) (holes) dropped across {features_with_dropped_holes} feature(s)."
        ));
    }
    if total_dropped_components > 0 {
        messages.push(format!(
            "{total_dropped_components} degenerate polygon component(s) dropped across {features_with_dropped_components} feature(s)."
        ));
    }

    Ok(ToolResult { output: Output::Layer(out_ds), messages })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_square(x: f64, y: f64, s: f64) -> Geometry {
        Geometry::Polygon(vec![vec![vec![
            Coord::new(x, y),
            Coord::new(x + s, y),
            Coord::new(x + s, y + s),
            Coord::new(x, y + s),
            Coord::new(x, y),
        ]]])
    }

    #[test]
    fn test_valid_point() {
        let p = Geometry::Point(Coord::new(10.0, 20.0));
        let rec = validate_geometry(Some(&p), 1, GeometryType::Point);
        assert_eq!(rec.status, ValidityStatus::Valid);
        assert!(rec.valid);
        assert_eq!(rec.error_category, ValidityCategory::None);

        let fix = fix_geometry(&p);
        assert!(matches!(fix, FixOutcome::Unchanged(_)));
    }

    #[test]
    fn test_invalid_point_non_finite() {
        let p = Geometry::Point(Coord::new(f64::NAN, 20.0));
        let rec = validate_geometry(Some(&p), 1, GeometryType::Point);
        assert_eq!(rec.status, ValidityStatus::Invalid);
        assert!(!rec.valid);
        assert_eq!(rec.error_category, ValidityCategory::NonFiniteCoord);

        let fix = fix_geometry(&p);
        assert!(matches!(fix, FixOutcome::Failed(_)));
    }

    #[test]
    fn test_valid_linestring() {
        let l = Geometry::Line(vec![vec![Coord::new(0.0, 0.0), Coord::new(1.0, 1.0), Coord::new(2.0, 0.0)]]);
        let rec = validate_geometry(Some(&l), 1, GeometryType::Polyline);
        assert_eq!(rec.status, ValidityStatus::Valid);
        assert!(rec.valid);

        let fix = fix_geometry(&l);
        assert!(matches!(fix, FixOutcome::Unchanged(_)));
    }

    #[test]
    fn test_invalid_linestring_too_few_points() {
        let l = Geometry::Line(vec![vec![Coord::new(1.0, 1.0), Coord::new(1.0, 1.0)]]);
        let rec = validate_geometry(Some(&l), 1, GeometryType::Polyline);
        assert_eq!(rec.status, ValidityStatus::Invalid);
        assert_eq!(rec.error_category, ValidityCategory::TooFewPoints);

        let fix = fix_geometry(&l);
        assert!(matches!(fix, FixOutcome::Failed(_)));
    }

    #[test]
    fn test_valid_polygon() {
        let poly = sample_square(0.0, 0.0, 10.0);
        let rec = validate_geometry(Some(&poly), 1, GeometryType::Polygon);
        assert_eq!(rec.status, ValidityStatus::Valid);
        assert!(rec.valid);

        let fix = fix_geometry(&poly);
        assert!(matches!(fix, FixOutcome::Unchanged(_)));
    }

    #[test]
    fn test_valid_polygon_with_hole() {
        let poly = Geometry::Polygon(vec![vec![
            vec![Coord::new(0.0, 0.0), Coord::new(10.0, 0.0), Coord::new(10.0, 10.0), Coord::new(0.0, 10.0), Coord::new(0.0, 0.0)],
            vec![Coord::new(2.0, 2.0), Coord::new(8.0, 2.0), Coord::new(8.0, 8.0), Coord::new(2.0, 8.0), Coord::new(2.0, 2.0)],
        ]]);
        let rec = validate_geometry(Some(&poly), 1, GeometryType::Polygon);
        assert_eq!(rec.status, ValidityStatus::Valid);
        assert!(rec.valid);

        let fix = fix_geometry(&poly);
        assert!(matches!(fix, FixOutcome::Unchanged(_)));
    }

    #[test]
    fn test_unclosed_ring_and_repair() {
        let unclosed =
            Geometry::Polygon(vec![vec![vec![Coord::new(0.0, 0.0), Coord::new(10.0, 0.0), Coord::new(10.0, 10.0), Coord::new(0.0, 10.0)]]]);
        let rec = validate_geometry(Some(&unclosed), 1, GeometryType::Polygon);
        assert_eq!(rec.status, ValidityStatus::Invalid);
        assert_eq!(rec.error_category, ValidityCategory::UnclosedRing);

        let fix = fix_geometry(&unclosed);
        let FixOutcome::Repaired(repaired) = fix else { panic!("expected Repaired, got {fix:?}") };
        let rec2 = validate_geometry(Some(&repaired), 1, GeometryType::Polygon);
        assert_eq!(rec2.status, ValidityStatus::Valid);
    }

    #[test]
    fn test_self_intersecting_polygon_bowtie_and_repair() {
        let bowtie = Geometry::Polygon(vec![vec![vec![
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 2.0),
            Coord::new(2.0, 0.0),
            Coord::new(0.0, 2.0),
            Coord::new(0.0, 0.0),
        ]]]);
        let rec = validate_geometry(Some(&bowtie), 1, GeometryType::Polygon);
        assert_eq!(rec.status, ValidityStatus::Invalid);
        assert_eq!(rec.error_category, ValidityCategory::SelfIntersection);
        assert!(rec.error_location.is_some());
        let loc = rec.error_location.unwrap();
        assert!((loc.x - 1.0).abs() < 1e-4 && (loc.y - 1.0).abs() < 1e-4);

        let fix = fix_geometry(&bowtie);
        let FixOutcome::Repaired(repaired) = fix else { panic!("expected Repaired, got {fix:?}") };
        let rec2 = validate_geometry(Some(&repaired), 1, GeometryType::Polygon);
        assert_eq!(rec2.status, ValidityStatus::Valid);

        let Geometry::Polygon(polys) = repaired else { panic!("expected Polygon") };
        assert_eq!(polys.len(), 2, "bowtie repairs into two valid triangles");
    }

    #[test]
    fn test_hole_outside_shell_and_repair() {
        let hole_outside = Geometry::Polygon(vec![vec![
            vec![Coord::new(0.0, 0.0), Coord::new(10.0, 0.0), Coord::new(10.0, 10.0), Coord::new(0.0, 10.0), Coord::new(0.0, 0.0)],
            vec![Coord::new(20.0, 20.0), Coord::new(25.0, 20.0), Coord::new(25.0, 25.0), Coord::new(20.0, 25.0), Coord::new(20.0, 20.0)],
        ]]);
        let rec = validate_geometry(Some(&hole_outside), 1, GeometryType::Polygon);
        assert_eq!(rec.status, ValidityStatus::Invalid);
        assert_eq!(rec.error_category, ValidityCategory::HoleOutsideShell);

        let fix = fix_geometry(&hole_outside);
        let FixOutcome::Repaired(repaired) = fix else { panic!("expected Repaired, got {fix:?}") };
        let rec2 = validate_geometry(Some(&repaired), 1, GeometryType::Polygon);
        assert_eq!(rec2.status, ValidityStatus::Valid);
    }

    #[test]
    fn test_empty_and_null_geometry() {
        let null_rec = validate_geometry(None, 1, GeometryType::Polygon);
        assert_eq!(null_rec.status, ValidityStatus::Null);
        assert_eq!(null_rec.error_category, ValidityCategory::MissingGeometry);

        let empty_poly = Geometry::Polygon(vec![]);
        let empty_rec = validate_geometry(Some(&empty_poly), 2, GeometryType::Polygon);
        assert_eq!(empty_rec.status, ValidityStatus::Empty);
        assert_eq!(empty_rec.error_category, ValidityCategory::EmptyGeometry);
    }

    #[test]
    fn test_valid_polygon_with_duplicates_is_unchanged() {
        let dup = Geometry::Polygon(vec![vec![vec![
            Coord::new(0.0, 0.0),
            Coord::new(0.0, 0.0),
            Coord::new(10.0, 0.0),
            Coord::new(10.0, 10.0),
            Coord::new(10.0, 10.0),
            Coord::new(0.0, 10.0),
            Coord::new(0.0, 0.0),
        ]]]);
        let fix = fix_geometry(&dup);
        // Valid geometry preserves bitwise identical coordinates without modification
        assert_eq!(fix, FixOutcome::Unchanged(dup));
    }

    #[test]
    fn test_invalid_polygon_with_duplicates_deduped_on_repair() {
        // An unclosed polygon that also has duplicate consecutive vertices
        let dup_unclosed = Geometry::Polygon(vec![vec![vec![
            Coord::new(0.0, 0.0),
            Coord::new(0.0, 0.0),
            Coord::new(10.0, 0.0),
            Coord::new(10.0, 10.0),
            Coord::new(10.0, 10.0),
            Coord::new(0.0, 10.0),
        ]]]);
        let fix = fix_geometry(&dup_unclosed);
        let FixOutcome::Repaired(repaired) = fix else { panic!("expected Repaired, got {fix:?}") };
        let Geometry::Polygon(polys) = &repaired else { panic!("expected Polygon") };
        assert_eq!(polys.len(), 1);
        let ring = &polys[0][0];
        // After repair: ring is closed and consecutive duplicates are removed (5 vertices: 4 distinct + closure)
        assert_eq!(ring.len(), 5);
        assert_eq!(ring.first(), ring.last());
        let rec = validate_geometry(Some(&repaired), 1, GeometryType::Polygon);
        assert_eq!(rec.status, ValidityStatus::Valid);
    }

    #[test]
    fn test_exterior_ring_collapse_fails_without_promoting_hole() {
        // Exterior ring is collapsed (all points identical [A, A, A, A]), but hole is valid
        let collapsed_exterior = Geometry::Polygon(vec![vec![
            vec![Coord::new(0.0, 0.0), Coord::new(0.0, 0.0), Coord::new(0.0, 0.0), Coord::new(0.0, 0.0)],
            vec![Coord::new(2.0, 2.0), Coord::new(8.0, 2.0), Coord::new(8.0, 8.0), Coord::new(2.0, 8.0), Coord::new(2.0, 2.0)],
        ]]);
        let fix = fix_geometry(&collapsed_exterior);
        // Must NOT return Repaired with the hole promoted to shell!
        assert!(matches!(fix, FixOutcome::Failed(_)), "expected Failed, got {fix:?}");
    }

    #[test]
    fn test_self_intersecting_polyline() {
        // Self-intersecting polyline
        let line = Geometry::Line(vec![vec![Coord::new(0.0, 0.0), Coord::new(2.0, 2.0), Coord::new(0.0, 2.0), Coord::new(2.0, 0.0)]]);
        let rec = validate_geometry(Some(&line), 1, GeometryType::Polyline);
        assert_eq!(rec.status, ValidityStatus::Invalid);
        assert_eq!(rec.error_category, ValidityCategory::SelfIntersection);
        assert!(rec.error_location.is_some());
        let loc = rec.error_location.unwrap();
        assert!((loc.x - 1.0).abs() < 1e-4 && (loc.y - 1.0).abs() < 1e-4);

        let fix = fix_geometry(&line);
        assert!(matches!(fix, FixOutcome::Failed(_)));
    }

    #[test]
    fn test_multipoint_empty_and_lossy_repair() {
        let empty_mp = Geometry::MultiPoint(vec![]);
        let (fix_empty, _) = fix_geometry_with_stats(&empty_mp);
        assert!(matches!(fix_empty, FixOutcome::Skipped(_)));

        let lossy_mp = Geometry::MultiPoint(vec![Coord::new(1.0, 2.0), Coord::new(f64::NAN, 3.0)]);
        let (fix_lossy, stats) = fix_geometry_with_stats(&lossy_mp);
        assert_eq!(stats.dropped_non_finite, 1);
        let FixOutcome::Repaired(Geometry::MultiPoint(clean)) = fix_lossy else { panic!("expected Repaired") };
        assert_eq!(clean, vec![Coord::new(1.0, 2.0)]);
    }

    #[test]
    fn test_line_part_dropped_stats() {
        // Multipart line with one valid part and one collapsed part [A, A]
        let line =
            Geometry::Line(vec![vec![Coord::new(0.0, 0.0), Coord::new(10.0, 10.0)], vec![Coord::new(5.0, 5.0), Coord::new(5.0, 5.0)]]);
        let (fix, stats) = fix_geometry_with_stats(&line);
        assert_eq!(stats.dropped_parts, 1);
        let FixOutcome::Repaired(Geometry::Line(parts)) = fix else { panic!("expected Repaired") };
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0].len(), 2);
    }
}
