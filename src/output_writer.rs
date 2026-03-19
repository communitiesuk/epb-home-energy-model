use formatx::formatx;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub trait OutputWriter: Debug + Sync + Send + Sized + Clone {
    fn writer_for_location_key(
        &self,
        location_key: &str,
        file_extension: &str,
    ) -> anyhow::Result<impl Write>;
    /// Whether this output can be considered a no-op and therefore that any code that only writes to the output can be skipped.
    fn is_noop(&self) -> bool {
        false
    }

    /// For writers that use a file template, create a writer with a template provided - else this is a no-op
    fn with_file_template(&self, _file_template: String) -> Self {
        self.clone()
    }
}

#[derive(Clone, Debug)]
pub struct FileOutputWriter {
    directory_path: PathBuf,
    file_template: String,
}

impl FileOutputWriter {
    pub fn new(directory_path: PathBuf, file_template: String) -> Self {
        Self {
            directory_path,
            file_template,
        }
    }
}

impl OutputWriter for FileOutputWriter {
    fn writer_for_location_key(
        &self,
        location_key: &str,
        file_extension: &str,
    ) -> anyhow::Result<impl Write> {
        Ok(BufWriter::new(File::create(self.directory_path.join(
            formatx!(&self.file_template, location_key, file_extension).unwrap(),
        ))?))
    }

    fn with_file_template(&self, file_template: String) -> Self {
        Self {
            directory_path: self.directory_path.clone(),
            file_template,
        }
    }
}

/// An output that goes to nowhere/ a "sink"/ /dev/null.
#[derive(Clone, Debug, Default)]
pub struct SinkOutputWriter;

impl OutputWriter for SinkOutputWriter {
    fn writer_for_location_key(
        &self,
        _location_key: &str,
        _file_extension: &str,
    ) -> anyhow::Result<impl Write> {
        Ok(io::sink())
    }

    fn is_noop(&self) -> bool {
        true
    }
}
