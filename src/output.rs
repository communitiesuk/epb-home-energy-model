use formatx::formatx;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub trait Output: Debug {
    fn writer_for_location_key(&self, location_key: &str) -> anyhow::Result<impl Write>;
    /// Whether this output can be considered a no-op and therefore that any code that only writes to the output can be skipped.
    fn is_noop(&self) -> bool {
        false
    }
}

#[derive(Debug)]
pub struct FileOutput {
    directory_path: PathBuf,
    file_template: String,
}

impl FileOutput {
    pub fn new(directory_path: PathBuf, file_template: String) -> Self {
        Self {
            directory_path,
            file_template,
        }
    }
}

impl Output for FileOutput {
    fn writer_for_location_key(&self, location_key: &str) -> anyhow::Result<impl Write> {
        Ok(BufWriter::new(File::create(self.directory_path.join(
            formatx!(&self.file_template, location_key).unwrap(),
        ))?))
    }
}

impl Output for &FileOutput {
    fn writer_for_location_key(&self, location_key: &str) -> anyhow::Result<impl Write> {
        <FileOutput as Output>::writer_for_location_key(self, location_key)
    }
}

/// An output that goes to nowhere/ a "sink"/ /dev/null.
#[derive(Debug, Default)]
pub struct SinkOutput;

impl Output for SinkOutput {
    fn writer_for_location_key(&self, _location_key: &str) -> anyhow::Result<impl Write> {
        Ok(io::sink())
    }

    fn is_noop(&self) -> bool {
        true
    }
}
