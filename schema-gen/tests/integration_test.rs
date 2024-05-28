use hem::input::Input;
use schemars::schema_for;

#[test]
fn test_generate_json_schema() {
    let schema = schema_for!(Input);
    assert!(serde_json::to_string_pretty(&schema).is_ok());
}
