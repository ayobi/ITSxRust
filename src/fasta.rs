use anyhow::{Context, Result, anyhow};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub struct FastaRecord {
    pub id: String,
    pub seq: String,
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f =
        File::open(path).with_context(|| format!("Failed to open input: {}", path.display()))?;
    if path.extension().map(|e| e == "gz").unwrap_or(false) {
        let gz = MultiGzDecoder::new(f);
        Ok(Box::new(BufReader::new(gz)))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

pub struct FastaReader {
    r: Box<dyn BufRead>,
    pending_header: Option<String>,
    buf: String,
    done: bool,
}

impl FastaReader {
    pub fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            r: open_maybe_gz(path)?,
            pending_header: None,
            buf: String::new(),
            done: false,
        })
    }

    pub fn next_record(&mut self) -> Result<Option<FastaRecord>> {
        if self.done {
            return Ok(None);
        }

        // find header if we don't already have one
        if self.pending_header.is_none() {
            loop {
                self.buf.clear();
                let n = self.r.read_line(&mut self.buf)?;
                if n == 0 {
                    self.done = true;
                    return Ok(None);
                }
                let line = self.buf.trim_end();
                if line.starts_with('>') {
                    let header = line[1..].trim();
                    let id = header
                        .split_whitespace()
                        .next()
                        .unwrap_or(header)
                        .to_string();
                    self.pending_header = Some(id);
                    break;
                }
            }
        }

        let id = self.pending_header.take().unwrap();
        let mut seq = String::new();

        loop {
            self.buf.clear();
            let n = self.r.read_line(&mut self.buf)?;
            if n == 0 {
                self.done = true;
                break;
            }
            let line = self.buf.trim_end();
            if line.starts_with('>') {
                let header = line[1..].trim();
                let next_id = header
                    .split_whitespace()
                    .next()
                    .unwrap_or(header)
                    .to_string();
                self.pending_header = Some(next_id);
                break;
            }
            seq.push_str(line);
        }

        Ok(Some(FastaRecord { id, seq }))
    }
}

pub fn open_writer(path: &Path) -> Result<BufWriter<File>> {
    let f = File::create(path)
        .with_context(|| format!("Failed to create output: {}", path.display()))?;
    Ok(BufWriter::new(f))
}

pub fn write_fasta_record<W: Write>(w: &mut W, id: &str, seq: &str) -> Result<()> {
    writeln!(w, ">{}", id)?;
    let bytes = seq.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let j = (i + 80).min(bytes.len());
        writeln!(
            w,
            "{}",
            std::str::from_utf8(&bytes[i..j]).map_err(|_| anyhow!("Invalid UTF-8 in seq"))?
        )?;
        i = j;
    }
    Ok(())
}
