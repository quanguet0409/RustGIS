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
