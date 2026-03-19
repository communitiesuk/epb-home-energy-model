use formatx::formatx;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub trait OutputWriter: Debug + Sync + Send {
    fn writer_for_location_key(
        &self,
        location_key: &str,
        file_extension: &str,
        include_filename_prefix: bool,
    ) -> anyhow::Result<impl Write>;
    /// Whether this output can be considered a no-op and therefore that any code that only writes to the output can be skipped.
    fn is_noop(&self) -> bool {
        false
    }
}

#[derive(Debug)]
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
        include_filename_prefix: bool,
    ) -> anyhow::Result<impl Write> {
        let prefix = if include_filename_prefix {
            &self.file_template
        } else {
            ""
        };

        Ok(BufWriter::new(File::create(
            self.directory_path
                .join(formatx!(prefix, location_key, file_extension)?),
        )?))
    }
}

impl OutputWriter for &FileOutputWriter {
    fn writer_for_location_key(
        &self,
        location_key: &str,
        file_extension: &str,
        include_filename_prefix: bool,
    ) -> anyhow::Result<impl Write> {
        <FileOutputWriter as OutputWriter>::writer_for_location_key(
            self,
            location_key,
            file_extension,
            include_filename_prefix,
        )
    }
}

/// An output that goes to nowhere/ a "sink"/ /dev/null.
#[derive(Debug, Default)]
pub struct SinkOutputWriter;

impl OutputWriter for SinkOutputWriter {
    fn writer_for_location_key(
        &self,
        _location_key: &str,
        _file_extension: &str,
        _include_filename_prefix: bool,
    ) -> anyhow::Result<impl Write> {
        Ok(io::sink())
    }

    fn is_noop(&self) -> bool {
        true
    }
}
