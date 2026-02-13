use anyhow::{Context, Result, anyhow};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct FastqRecord {
    pub id: String, // first token after '@'
    pub seq: String,
    pub qual: String,
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

pub struct FastqReader {
    r: Box<dyn BufRead>,
    buf: String,
}

impl FastqReader {
    pub fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            r: open_maybe_gz(path)?,
            buf: String::new(),
        })
    }

    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        self.buf.clear();
        let n = self.r.read_line(&mut self.buf)?;
        if n == 0 {
            return Ok(None);
        }
        let h = self.buf.trim_end().to_string();
        if !h.starts_with('@') {
            return Err(anyhow!("FASTQ header does not start with '@': {}", h));
        }
        let id = h[1..].split_whitespace().next().unwrap_or("").to_string();
        if id.is_empty() {
            return Err(anyhow!("Empty FASTQ id in header: {}", h));
        }

        // seq
        self.buf.clear();
        self.r
            .read_line(&mut self.buf)
            .context("Missing FASTQ sequence line")?;
        let seq = self.buf.trim_end().to_string();

        // plus
        self.buf.clear();
        self.r
            .read_line(&mut self.buf)
            .context("Missing FASTQ '+' line")?;
        let plus = self.buf.trim_end().to_string();
        if !plus.starts_with('+') {
            return Err(anyhow!(
                "FASTQ third line does not start with '+': {}",
                plus
            ));
        }

        // qual
        self.buf.clear();
        self.r
            .read_line(&mut self.buf)
            .context("Missing FASTQ quality line")?;
        let qual = self.buf.trim_end().to_string();

        if seq.len() != qual.len() {
            return Err(anyhow!(
                "FASTQ seq/qual length mismatch for {}: seq={} qual={}",
                id,
                seq.len(),
                qual.len()
            ));
        }

        Ok(Some(FastqRecord { id, seq, qual }))
    }
}

pub fn write_fastq_record<W: Write>(w: &mut W, id: &str, seq: &str, qual: &str) -> Result<()> {
    if seq.len() != qual.len() {
        return Err(anyhow!(
            "Cannot write FASTQ: seq/qual length mismatch for {}",
            id
        ));
    }
    writeln!(w, "@{}", id)?;
    writeln!(w, "{}", seq)?;
    writeln!(w, "+")?;
    writeln!(w, "{}", qual)?;
    Ok(())
}

pub fn open_writer(path: &Path) -> Result<BufWriter<File>> {
    let f = File::create(path)
        .with_context(|| format!("Failed to create output: {}", path.display()))?;
    Ok(BufWriter::new(f))
}
