// src/seq.rs
// Small utilities for reverse-complementing DNA and handling FASTQ qualities.

/// Reverse-complement a DNA sequence.
///
/// - Preserves case (A<->T, C<->G, plus common IUPAC ambiguity codes).
/// - Unknown characters are passed through unchanged.
pub fn revcomp_dna(seq: &str) -> String {
    let mut out = String::with_capacity(seq.len());

    for b in seq.as_bytes().iter().rev().copied() {
        let c = match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'U' => 'A',

            b'a' => 't',
            b'c' => 'g',
            b'g' => 'c',
            b't' => 'a',
            b'u' => 'a',

            // IUPAC ambiguity codes
            b'R' => 'Y', // A/G <-> C/T
            b'Y' => 'R',
            b'S' => 'S', // C/G
            b'W' => 'W', // A/T
            b'K' => 'M', // G/T <-> A/C
            b'M' => 'K',
            b'B' => 'V', // C/G/T <-> A/C/G
            b'V' => 'B',
            b'D' => 'H', // A/G/T <-> A/C/T
            b'H' => 'D',
            b'N' => 'N',

            b'r' => 'y',
            b'y' => 'r',
            b's' => 's',
            b'w' => 'w',
            b'k' => 'm',
            b'm' => 'k',
            b'b' => 'v',
            b'v' => 'b',
            b'd' => 'h',
            b'h' => 'd',
            b'n' => 'n',

            _ => b as char,
        };
        out.push(c);
    }

    out
}

/// Reverse a FASTQ quality string (needed when reverse-complementing the sequence).
pub fn reverse_qual(qual: &str) -> String {
    let bytes = qual.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    out.extend(bytes.iter().rev().copied());
    // FASTQ quality strings are ASCII; this should always be valid UTF-8.
    String::from_utf8(out).unwrap_or_else(|_| qual.chars().rev().collect())
}
