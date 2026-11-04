//! Integration tests for Check Geometry and Fix Geometries tools.

use rustgis_analysis::{Inputs, Output};
use rustgis_data::raster::Crs;
use rustgis_data::{Coord, Dataset, Field, FieldType, Geometry, GeometryType, Value};
use serde_json::json;

struct TestLayer {
    ds: Dataset,
    selected: Vec<usize>,
}

impl Inputs for TestLayer {
    fn dataset(&self, name: &str) -> Option<&Dataset> {
        (name == "input" || name == self.ds.name.as_str()).then_some(&self.ds)
    }
    fn selection(&self, _: &str) -> &[usize] {
        &self.selected
    }
}

fn create_layer(name: &str, ty: GeometryType) -> TestLayer {
    TestLayer {
        ds: Dataset::new(name, ty, vec![Field::new("name", FieldType::Text), Field::new("value", FieldType::Integer)]),
        selected: vec![],
    }
}

#[test]
fn test_check_geometry_valid_and_invalid_polygons() {
    let mut l = create_layer("parcels", GeometryType::Polygon);

    // Feature 0: Valid square
    let square = Geometry::Polygon(vec![vec![vec![
        Coord::new(0.0, 0.0),
        Coord::new(10.0, 0.0),
        Coord::new(10.0, 10.0),
        Coord::new(0.0, 10.0),
        Coord::new(0.0, 0.0),
    ]]]);
    l.ds.push_next(Some(square), vec![Value::Text("valid_sq".into()), Value::Int(100)]);

    // Feature 1: Bowtie self-intersecting polygon
    let bowtie = Geometry::Polygon(vec![vec![vec![
        Coord::new(0.0, 0.0),
        Coord::new(2.0, 2.0),
        Coord::new(2.0, 0.0),
        Coord::new(0.0, 2.0),
        Coord::new(0.0, 0.0),
    ]]]);
    l.ds.push_next(Some(bowtie), vec![Value::Text("bowtie".into()), Value::Int(200)]);

    // Feature 2: Unclosed ring
    let unclosed =
        Geometry::Polygon(vec![vec![vec![Coord::new(0.0, 0.0), Coord::new(5.0, 0.0), Coord::new(5.0, 5.0), Coord::new(0.0, 5.0)]]]);
    l.ds.push_next(Some(unclosed), vec![Value::Text("unclosed".into()), Value::Int(300)]);

    // Feature 3: Null geometry
    l.ds.push_next(None, vec![Value::Text("null_geom".into()), Value::Int(400)]);

    // Feature 4: Empty geometry
    l.ds.push_next(Some(Geometry::Polygon(vec![])), vec![Value::Text("empty_geom".into()), Value::Int(500)]);

    // Run checkGeometry
    let res =
        rustgis_analysis::run("checkGeometry", &l, &json!({"input": "parcels", "output": "report"})).expect("checkGeometry should succeed");

    let Output::Table(table) = res.output else { panic!("expected Output::Table") };
    assert_eq!(table.features.len(), 5);

    // Verify row 0: valid square
    assert_eq!(table.features[0].attrs[2], Value::Text("Valid".into()));
    assert_eq!(table.features[0].attrs[3], Value::Bool(true));
    assert_eq!(table.features[0].attrs[4], Value::Text("None".into()));

    // Verify row 1: bowtie self-intersection
    assert_eq!(table.features[1].attrs[2], Value::Text("Invalid".into()));
    assert_eq!(table.features[1].attrs[3], Value::Bool(false));
    assert_eq!(table.features[1].attrs[4], Value::Text("Self-Intersection".into()));
    let Value::Float(err_x) = table.features[1].attrs[6] else { panic!("expected float x") };
    let Value::Float(err_y) = table.features[1].attrs[7] else { panic!("expected float y") };
    assert!((err_x - 1.0).abs() < 1e-4);
    assert!((err_y - 1.0).abs() < 1e-4);

    // Verify row 2: unclosed ring
    assert_eq!(table.features[2].attrs[2], Value::Text("Invalid".into()));
    assert_eq!(table.features[2].attrs[3], Value::Bool(false));
    assert_eq!(table.features[2].attrs[4], Value::Text("Unclosed Ring".into()));

    // Verify row 3: null geometry
    assert_eq!(table.features[3].attrs[2], Value::Text("Null".into()));
    assert_eq!(table.features[3].attrs[4], Value::Text("Missing Geometry".into()));

    // Verify row 4: empty geometry
    assert_eq!(table.features[4].attrs[2], Value::Text("Empty".into()));
    assert_eq!(table.features[4].attrs[4], Value::Text("Empty Geometry".into()));

    // Verify messages summary
    assert!(res.messages[0].contains("5 feature(s) checked"));
    assert!(res.messages[0].contains("1 valid, 2 invalid, 1 empty, 1 null"));
}

#[test]
fn test_fix_geometries_workflow() {
    let mut l = create_layer("parcels", GeometryType::Polygon);
    l.ds.crs = Some(Crs::Epsg(4326));
    l.ds.crs_incomplete = true;
    l.ds.zm_dropped = true;

    // Feature 0: Valid square
    let square = Geometry::Polygon(vec![vec![vec![
        Coord::new(0.0, 0.0),
        Coord::new(10.0, 0.0),
        Coord::new(10.0, 10.0),
        Coord::new(0.0, 10.0),
        Coord::new(0.0, 0.0),
    ]]]);
    l.ds.push_next(Some(square.clone()), vec![Value::Text("valid_sq".into()), Value::Int(100)]);

    // Feature 1: Bowtie
    let bowtie = Geometry::Polygon(vec![vec![vec![
        Coord::new(0.0, 0.0),
        Coord::new(2.0, 2.0),
        Coord::new(2.0, 0.0),
        Coord::new(0.0, 2.0),
        Coord::new(0.0, 0.0),
    ]]]);
    l.ds.push_next(Some(bowtie), vec![Value::Text("bowtie".into()), Value::Int(200)]);

    // Feature 2: Unclosed ring
    let unclosed =
        Geometry::Polygon(vec![vec![vec![Coord::new(0.0, 0.0), Coord::new(5.0, 0.0), Coord::new(5.0, 5.0), Coord::new(0.0, 5.0)]]]);
    l.ds.push_next(Some(unclosed), vec![Value::Text("unclosed".into()), Value::Int(300)]);

    // Feature 3: Null geometry
    l.ds.push_next(None, vec![Value::Text("null_geom".into()), Value::Int(400)]);

    // Run fixGeometries
    let res = rustgis_analysis::run("fixGeometries", &l, &json!({"input": "parcels", "output": "parcels_fixed"}))
        .expect("fixGeometries should succeed");

    let Output::Layer(out) = res.output else { panic!("expected Output::Layer") };
    assert_eq!(out.name, "parcels_fixed");
    assert_eq!(out.crs, Some(Crs::Epsg(4326)));
    assert!(out.crs_incomplete);
    assert!(out.zm_dropped);
    assert_eq!(out.features.len(), 4);

    // Feature 0: unchanged, attributes intact
    assert_eq!(out.features[0].geometry, Some(square));
    assert_eq!(out.features[0].attrs, vec![Value::Text("valid_sq".into()), Value::Int(100)]);

    // Feature 1: repaired into 2 triangles, attributes intact
    assert_eq!(out.features[1].attrs, vec![Value::Text("bowtie".into()), Value::Int(200)]);
    let Some(Geometry::Polygon(bowtie_parts)) = &out.features[1].geometry else { panic!("expected repaired polygon") };
    assert_eq!(bowtie_parts.len(), 2, "bowtie decomposed into 2 triangles");

    // Feature 2: repaired into closed ring, attributes intact
    assert_eq!(out.features[2].attrs, vec![Value::Text("unclosed".into()), Value::Int(300)]);
    let Some(Geometry::Polygon(unclosed_parts)) = &out.features[2].geometry else { panic!("expected repaired polygon") };
    assert_eq!(unclosed_parts[0][0].first(), unclosed_parts[0][0].last());

    // Feature 3: null geometry skipped, attributes intact
    assert_eq!(out.features[3].geometry, None);
    assert_eq!(out.features[3].attrs, vec![Value::Text("null_geom".into()), Value::Int(400)]);

    // Verify messages
    assert!(res.messages[0].contains("Fix Geometries complete: 4 feature(s) processed (1 unchanged, 2 repaired, 0 failed, 1 skipped)"));

    // Now validate the repaired dataset: all repaired features must be valid!
    let fixed_layer = TestLayer { ds: out, selected: vec![] };
    let check_res = rustgis_analysis::run("checkGeometry", &fixed_layer, &json!({"input": "parcels_fixed"}))
        .expect("checkGeometry on repaired layer should succeed");
    let Output::Table(report) = check_res.output else { panic!("expected Output::Table") };

    // Repaired features 0, 1, 2 must all be Valid
    assert_eq!(report.features[0].attrs[2], Value::Text("Valid".into()));
    assert_eq!(report.features[1].attrs[2], Value::Text("Valid".into()));
    assert_eq!(report.features[2].attrs[2], Value::Text("Valid".into()));
}

#[test]
fn test_selection_filtering() {
    let mut l = create_layer("lines", GeometryType::Polyline);
    let line1 = Geometry::Line(vec![vec![Coord::new(0.0, 0.0), Coord::new(1.0, 1.0)]]);
    let line2 = Geometry::Line(vec![vec![Coord::new(2.0, 2.0), Coord::new(3.0, 3.0)]]);
    l.ds.push_next(Some(line1), vec![Value::Text("line1".into()), Value::Int(1)]);
    l.ds.push_next(Some(line2), vec![Value::Text("line2".into()), Value::Int(2)]);

    // Select only feature 1
    l.selected = vec![1];

    let res = rustgis_analysis::run("checkGeometry", &l, &json!({"input": "lines"})).unwrap();
    let Output::Table(table) = res.output else { panic!("expected Output::Table") };
    assert_eq!(table.features.len(), 1);
    assert_eq!(table.features[0].attrs[0], Value::Int(l.ds.features[1].id as i64));
}

#[test]
fn test_unrepairable_geometry_preserves_feature() {
    let mut l = create_layer("points", GeometryType::Point);
    let non_finite_pt = Geometry::Point(Coord::new(f64::NAN, 10.0));
    l.ds.push_next(Some(non_finite_pt.clone()), vec![Value::Text("bad_pt".into()), Value::Int(999)]);

    let res = rustgis_analysis::run("fixGeometries", &l, &json!({"input": "points"})).unwrap();
    let Output::Layer(out) = res.output else { panic!("expected Output::Layer") };

    // Feature is kept, attributes intact, original geometry preserved rather than dropped
    assert_eq!(out.features.len(), 1);
    assert_eq!(out.features[0].attrs, vec![Value::Text("bad_pt".into()), Value::Int(999)]);
    assert!(res.messages[0].contains("1 failed"));
}

#[test]
fn test_lossy_notices_and_clean_notices() {
    // 1. Clean dataset: must produce NO lossy notices
    let mut clean_layer = create_layer("clean", GeometryType::Polygon);
    let clean_square = Geometry::Polygon(vec![vec![vec![
        Coord::new(0.0, 0.0),
        Coord::new(10.0, 0.0),
        Coord::new(10.0, 10.0),
        Coord::new(0.0, 10.0),
        Coord::new(0.0, 0.0),
    ]]]);
    clean_layer.ds.push_next(Some(clean_square), vec![Value::Text("clean".into()), Value::Int(1)]);

    let clean_res = rustgis_analysis::run("fixGeometries", &clean_layer, &json!({"input": "clean"}))
        .expect("fixGeometries on clean input should succeed");
    assert_eq!(clean_res.messages.len(), 1, "clean input should only have summary message");
    assert!(!clean_res.messages[0].contains("dropped"));

    // 2. Lossy dataset: polygon with a NaN coordinate and a degenerate hole, plus a multipart line with a collapsed part
    let mut lossy_poly_layer = create_layer("lossy_poly", GeometryType::Polygon);
    let lossy_poly = Geometry::Polygon(vec![vec![
        // Exterior ring with a NaN coordinate (repairable polygon: 4 valid coords + 1 NaN)
        vec![
            Coord::new(0.0, 0.0),
            Coord::new(10.0, 0.0),
            Coord::new(f64::NAN, 5.0),
            Coord::new(10.0, 10.0),
            Coord::new(0.0, 10.0),
            Coord::new(0.0, 0.0),
        ],
        // Degenerate hole (fewer than 3 distinct points)
        vec![Coord::new(5.0, 5.0), Coord::new(5.0, 5.0), Coord::new(5.0, 5.0), Coord::new(5.0, 5.0)],
    ]]);
    lossy_poly_layer.ds.push_next(Some(lossy_poly), vec![Value::Text("lossy_poly".into()), Value::Int(2)]);

    let poly_res = rustgis_analysis::run("fixGeometries", &lossy_poly_layer, &json!({"input": "lossy_poly"}))
        .expect("fixGeometries on lossy polygon should succeed");
    let all_poly_msgs = poly_res.messages.join("\n");
    assert!(all_poly_msgs.contains("1 non-finite coordinate(s) dropped across 1 feature(s)."));
    assert!(all_poly_msgs.contains("1 degenerate interior ring(s) (holes) dropped across 1 feature(s)."));

    // 3. Lossy line layer: multi-part line with one valid part and one collapsed part
    let mut lossy_line_layer = create_layer("lossy_line", GeometryType::Polyline);
    let lossy_line =
        Geometry::Line(vec![vec![Coord::new(0.0, 0.0), Coord::new(10.0, 10.0)], vec![Coord::new(5.0, 5.0), Coord::new(5.0, 5.0)]]);
    lossy_line_layer.ds.push_next(Some(lossy_line), vec![Value::Text("lossy_line".into()), Value::Int(3)]);

    let line_res = rustgis_analysis::run("fixGeometries", &lossy_line_layer, &json!({"input": "lossy_line"}))
        .expect("fixGeometries on lossy line should succeed");
    let all_line_msgs = line_res.messages.join("\n");
    assert!(all_line_msgs.contains("1 degenerate line part(s) (fewer than 2 points) dropped across 1 feature(s)."));
}

#[test]
fn test_exterior_ring_collapse_fails_in_fix_geometries() {
    let mut l = create_layer("collapsed_exterior", GeometryType::Polygon);
    let poly = Geometry::Polygon(vec![vec![
        vec![Coord::new(0.0, 0.0), Coord::new(0.0, 0.0), Coord::new(0.0, 0.0), Coord::new(0.0, 0.0)],
        vec![Coord::new(2.0, 2.0), Coord::new(8.0, 2.0), Coord::new(8.0, 8.0), Coord::new(2.0, 8.0), Coord::new(2.0, 2.0)],
    ]]);
    l.ds.push_next(Some(poly), vec![Value::Text("collapsed".into()), Value::Int(10)]);

    let res = rustgis_analysis::run("fixGeometries", &l, &json!({"input": "collapsed_exterior"})).unwrap();
    assert!(res.messages[0].contains("1 failed"));
    let Output::Layer(out) = res.output else { panic!("expected Output::Layer") };
    // Original geometry preserved, hole was NOT promoted to exterior!
    let Some(Geometry::Polygon(parts)) = &out.features[0].geometry else { panic!("expected polygon") };
    assert_eq!(parts[0][0][0], Coord::new(0.0, 0.0));
}

#[test]
fn test_polyline_self_intersection_check_and_fix() {
    let mut l = create_layer("crossing_line", GeometryType::Polyline);
    let line = Geometry::Line(vec![vec![Coord::new(0.0, 0.0), Coord::new(2.0, 2.0), Coord::new(0.0, 2.0), Coord::new(2.0, 0.0)]]);
    l.ds.push_next(Some(line.clone()), vec![Value::Text("cross".into()), Value::Int(1)]);

    // Check Geometry detects self-intersection
    let check_res = rustgis_analysis::run("checkGeometry", &l, &json!({"input": "crossing_line"})).unwrap();
    let Output::Table(table) = check_res.output else { panic!("expected table") };
    assert_eq!(table.features[0].attrs[2], Value::Text("Invalid".into()));
    assert_eq!(table.features[0].attrs[4], Value::Text("Self-Intersection".into()));

    // Fix Geometries fails repair and preserves original geometry
    let fix_res = rustgis_analysis::run("fixGeometries", &l, &json!({"input": "crossing_line"})).unwrap();
    assert!(fix_res.messages[0].contains("1 failed"));
    let Output::Layer(out) = fix_res.output else { panic!("expected layer") };
    assert_eq!(out.features[0].geometry, Some(line));
}
