use std::collections::HashMap;
use std::fs;
use std::process::Command;
use tempfile::tempdir;

fn count_fasta_headers(path: &std::path::Path) -> usize {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .count()
}

#[test]
fn extract_subset_full_from_existing_tblout() {
    let td = tempdir().unwrap();
    let out = td.path().join("out.fasta");

    let status = Command::new(env!("CARGO_BIN_EXE_itsxrust"))
        .args([
            "extract",
            "--input",
            "data/subset.fasta",
            "--tblout-existing",
            "data/subset.tblout",
            "--output",
            out.to_str().unwrap(),
            "--region",
            "full",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    assert_eq!(count_fasta_headers(&out), 766);
}

fn read_fasta_as_map(path: &std::path::Path) -> HashMap<String, String> {
    let text = fs::read_to_string(path).unwrap();
    let mut map = HashMap::new();
    let mut cur_id: Option<String> = None;
    let mut seq = String::new();

    for line in text.lines() {
        if line.starts_with('>') {
            if let Some(id) = cur_id.take() {
                map.insert(id, seq.clone());
            }
            seq.clear();
            // output header is like "READ|region:..." -> store by READ
            let header = line[1..].split_whitespace().next().unwrap();
            let rid = header.split('|').next().unwrap().to_string();
            cur_id = Some(rid);
        } else {
            seq.push_str(line.trim());
        }
    }
    if let Some(id) = cur_id.take() {
        map.insert(id, seq);
    }
    map
}

fn revcomp_dna(s: &str) -> String {
    fn comp(b: u8) -> u8 {
        match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            _ => b'N',
        }
    }
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        out.push(comp(b));
    }
    String::from_utf8(out).unwrap()
}

#[test]
fn minus_strand_full_matches_expected_revcomp_slice() {
    let td = tempdir().unwrap();
    let out = td.path().join("out.full.fasta");

    // Run extraction (existing tblout, deterministic)
    let status = Command::new(env!("CARGO_BIN_EXE_itsxrust"))
        .args([
            "extract",
            "--input",
            "data/subset.fasta",
            "--tblout-existing",
            "data/subset.tblout",
            "--output",
            out.to_str().unwrap(),
            "--region",
            "full",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // Known '-' strand read and bounds from your earlier validation
    let read_id = "SRR21494940.227";
    let full_start: usize = 32;
    let full_end: usize = 814;

    // Load original and output sequences
    let orig_map = read_fasta_as_map(std::path::Path::new("data/subset.fasta"));
    let out_map = read_fasta_as_map(&out);

    let orig = orig_map
        .get(read_id)
        .expect("read not found in original FASTA");
    let got = out_map
        .get(read_id)
        .expect("read not found in extracted FASTA");

    let l = orig.len();
    // map normalized bounds back to original coords (1-based inclusive)
    let orig_start = l - full_end + 1;
    let orig_end = l - full_start + 1;

    let slice = &orig[(orig_start - 1)..orig_end];
    let expected = revcomp_dna(slice);

    assert_eq!(got.len(), expected.len());
    assert_eq!(got, &expected);
}
