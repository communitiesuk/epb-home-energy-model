use hem::output::Output;
use hem::read_weather_file::weather_data_to_vec;
use hem::run_project;
use lambda_http::{run, service_fn, tracing, Body, Error, Request, Response};
use parking_lot::Mutex;
use std::io;
use std::io::{BufReader, Cursor, ErrorKind, Write};
use std::str::from_utf8;
use std::sync::Arc;

async fn function_handler(event: Request) -> Result<Response<Body>, Error> {
    // Extract some useful information from the request
    let input = match event.body() {
        Body::Empty => "",
        Body::Text(text) => text.as_str(),
        Body::Binary(_) => unimplemented!(),
    }
    .as_bytes();

    let output = LambdaOutput::new();

    let external_conditions = weather_data_to_vec(BufReader::new(Cursor::new(include_str!(
        "../../src/weather.epw"
    ))))
    .ok();

    run_project(
        input,
        &output,
        external_conditions,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
    )?;

    // Return something that implements IntoResponse.
    // It will be serialized to the right response event automatically by the runtime
    let resp = Response::builder()
        .status(200)
        .header("content-type", "text/plain")
        .body(output.into())
        .map_err(Box::new)?;
    Ok(resp)
}

#[tokio::main]
async fn main() -> Result<(), Error> {
    tracing::init_default_subscriber();

    run(service_fn(function_handler)).await
}

/// This output uses a shared string that individual "file" writers (the FileLikeStringWriter type)
/// can write to - this string can then be used as the response body for the Lambda.
#[derive(Debug)]
struct LambdaOutput(Arc<Mutex<String>>);

impl LambdaOutput {
    fn new() -> Self {
        Self(Arc::new(Mutex::new(String::with_capacity(
            // output is expected to be about 4MB so allocate this up front
            2usize.pow(22),
        ))))
    }
}

impl Output for LambdaOutput {
    fn writer_for_location_key(&self, location_key: &str) -> anyhow::Result<impl Write> {
        Ok(FileLikeStringWriter::new(
            self.0.clone(),
            location_key.to_string(),
        ))
    }
}

impl Output for &LambdaOutput {
    fn writer_for_location_key(&self, location_key: &str) -> anyhow::Result<impl Write> {
        <LambdaOutput as Output>::writer_for_location_key(self, location_key)
    }
}

impl From<LambdaOutput> for Body {
    fn from(value: LambdaOutput) -> Self {
        Arc::try_unwrap(value.0).unwrap().into_inner().into()
    }
}

/// Represents a writer for an individual "file".
struct FileLikeStringWriter {
    string: Arc<Mutex<String>>,
    location_key: String,
    has_output_file_header: bool,
}

impl FileLikeStringWriter {
    fn new(string: Arc<Mutex<String>>, location_key: String) -> Self {
        Self {
            string,
            location_key,
            has_output_file_header: false,
        }
    }
}

impl Write for FileLikeStringWriter {
    /// Writes out bytes to this "file" (part of the wider LambdaOutput string), making sure there is
    /// a human-readable header at the start of the file so a human can know what each part of the output
    /// is sourced from.
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if !self.has_output_file_header {
            let mut output_string = self.string.lock();
            if !output_string.is_empty() {
                output_string.push_str("\n\n");
            }
            output_string
                .push_str(format!("Writing out file '{}':\n\n", self.location_key).as_str());
            self.has_output_file_header = true;
        }
        let utf8 = match from_utf8(buf) {
            Ok(utf8) => utf8,
            Err(_) => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Tried to write out invalid UTF-8.",
                ));
            }
        };
        self.string.lock().push_str(utf8);
        Ok(utf8.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
