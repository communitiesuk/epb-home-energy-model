use hem::input::Input;
use schemars::schema_for;

fn main() {
    let schema = schema_for!(Input);
    println!("{}", serde_json::to_string_pretty(&schema).unwrap());
}
