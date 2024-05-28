use formatx::formatx;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub trait Output {
    fn writer_for_location_key(&self, location_key: &str) -> anyhow::Result<impl Write>;
}

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
