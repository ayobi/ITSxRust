use anyhow::{Context, Result, anyhow};
use std::path::Path;
use std::process::{Command, Stdio};

pub struct HmmerArgs<'a> {
    pub hmm: &'a Path,
    pub fasta: &'a Path,
    pub tblout: &'a Path,
    pub cpu: usize,

    /// Reporting threshold (nhmmer -E). Controls what gets printed to --tblout.
    pub report_e: f64,

    /// Inclusion threshold (nhmmer --incE). Used for “included” hits (mainly for pipelines),
    /// but not sufficient alone to clean up --tblout.
    pub inc_e: f64,
}

pub fn run_nhmmer(args: &HmmerArgs) -> Result<()> {
    // nhmmer writes a human-readable report to stdout in addition to --tblout.
    // Capturing that stream (Command::output()) buffers the whole report in
    // memory and it scales with the number of *reported* hits: measured at
    // ~464 B/hit, it accounted for 5.98 GB of a 7.07 GB peak RSS on a 230 Mbp
    // input, versus 1.09 GB for the identical parse from a pre-computed tblout.
    // We only ever read --tblout, so discard stdout and keep stderr for errors.
    let output = Command::new("nhmmer")
        .arg("--cpu")
        .arg(args.cpu.to_string())
        .arg("--noali")
        // Reporting threshold: this is the key for cleaning up --tblout.
        .arg("-E")
        .arg(format!("{}", args.report_e))
        // Inclusion threshold: keep same as report threshold by default.
        .arg("--incE")
        .arg(format!("{}", args.inc_e))
        .arg("--tblout")
        .arg(args.tblout)
        .arg(args.hmm)
        .arg(args.fasta)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .with_context(|| "Failed to spawn `nhmmer`. Is your conda env active?")?
        .wait_with_output()
        .with_context(|| "`nhmmer` terminated abnormally")?;

    if !output.status.success() {
        return Err(anyhow!(
            "nhmmer failed (exit={:?})\n--- stderr ---\n{}",
            output.status.code(),
            String::from_utf8_lossy(&output.stderr),
        ));
    }

    Ok(())
}
