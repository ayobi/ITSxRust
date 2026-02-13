// src/derep.rs
//
// Exact dereplication for reducing redundant HMM searches.
//
// Sequences are grouped by exact match (case-insensitive). Only one representative
// per group is written to the nhmmer target FASTA. After tblout parsing, results are
// projected back to all duplicates via `expand_topk_map`.

use anyhow::Result;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use crate::fasta;
use crate::trim;
use crate::tblout::TopKHits;
use crate::InputFormat;

/// Result of a dereplication pass.
pub struct DerepResult {
    /// representative_id → list of all member IDs (representative is always first).
    pub rep_to_members: HashMap<String, Vec<String>>,
    pub total_seqs: usize,
    pub unique_seqs: usize,
}

/// Read all sequences from `input`, group by exact (case-insensitive) sequence,
/// and write one representative per group to `out_fasta`.
///
/// Returns a `DerepResult` whose `rep_to_members` map is used after tblout parsing
/// to project HMM hits from representatives back to all duplicate members.
pub fn derep_and_write_fasta(
    input: &Path,
    in_fmt: InputFormat,
    out_fasta: &Path,
) -> Result<DerepResult> {
    // Map: uppercase sequence → (representative_id, original_seq, [all member IDs])
    //
    // We store the representative's original-case sequence so nhmmer sees the same
    // bases as in the real input. For typical amplicon datasets the number of unique
    // sequences is a small fraction of total reads, keeping memory manageable.
    let mut groups: HashMap<String, (String, String, Vec<String>)> = HashMap::new();
    let mut total: usize = 0;

    let mut reader = trim::SeqReader::open(input, in_fmt)?;
    while let Some(rec) = reader.next_record()? {
        total += 1;
        let key = rec.seq.to_ascii_uppercase();
        groups
            .entry(key)
            .and_modify(|(_, _, members)| members.push(rec.id.clone()))
            .or_insert_with(|| (rec.id.clone(), rec.seq.clone(), vec![rec.id.clone()]));
    }

    let unique = groups.len();

    // Write one representative per group
    let mut w = fasta::open_writer(out_fasta)?;
    let mut rep_to_members: HashMap<String, Vec<String>> = HashMap::with_capacity(unique);

    for (_key, (rep_id, rep_seq, members)) in &groups {
        fasta::write_fasta_record(&mut w, rep_id, rep_seq)?;
        rep_to_members.insert(rep_id.clone(), members.clone());
    }
    w.flush()?;

    Ok(DerepResult {
        rep_to_members,
        total_seqs: total,
        unique_seqs: unique,
    })
}

/// After tblout parsing, the `topk_map` contains entries only for representative IDs.
/// This function copies each representative's hits to every duplicate member, so the
/// downstream plan computation and trimming sees all original read IDs.
///
/// Hit records for duplicates have their `read_id` updated to the member's ID.
pub fn expand_topk_map(
    topk_map: &mut HashMap<String, TopKHits>,
    derep: &DerepResult,
) {
    let mut additions: Vec<(String, TopKHits)> = Vec::new();

    for (rep_id, members) in &derep.rep_to_members {
        if let Some(rep_hits) = topk_map.get(rep_id) {
            // Skip the first member (the representative itself — already in the map)
            for member_id in members.iter().skip(1) {
                let mut member_hits = rep_hits.clone();
                // Update read_id in every stored hit for consistency
                for bin in &mut member_hits.bins {
                    for hit in bin.iter_mut() {
                        hit.read_id = member_id.clone();
                    }
                }
                additions.push((member_id.clone(), member_hits));
            }
        }
    }

    topk_map.extend(additions);
}