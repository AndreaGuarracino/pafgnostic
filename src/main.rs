use clap::Parser;
use flate2::read::GzDecoder;
use regex::Regex;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// The input PAF file (uncompressed or gzipped).
    #[clap(short='p', long, value_parser)]
    paf: String,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let file_path = &args.paf;
    let file = File::open(file_path)?;
    let buf_reader = BufReader::new(file);

    // Check if the file is gzipped based on the filename extension
    let reader: Box<dyn BufRead> = if file_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(buf_reader)))
    } else {
        Box::new(buf_reader)
    };

    let cigar_regex = Regex::new(r"(\d+)([=XID])").unwrap();

    // Print the header
    println!("\
        q.name\tq.start\tq.end\tq.strand\t\
        t.name\tt.start\tt.end\t\
        num.=\tnum.X\t\
        num.I\tnum.D\t\
        unique.=\tunique.X\t\
        unique.I\tunique.D\t\
        longest.=\tshortest.=\tavg.=\t\
        longest.X\tshortest.X\tavg.X\t\
        longest.I\tshortest.I\tavg.I\t\
        longest.D\tshortest.D\tavg.D"
    );

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 11 {
            continue; // Skip lines that don't have enough fields to contain the required data
        }
    
        let query_name = fields[0];
        let query_start = fields[2];
        let query_end = fields[3];
        let strand = fields[4];
        let target_name = fields[5];
        let target_start = fields[7];
        let target_end = fields[8];
    
        if let Some(cigar_field) = fields.iter().find(|&f| f.starts_with("cg:Z:")) {
            let cigar_str = &cigar_field[5..]; // Remove the 'cg:Z:' prefix
            let (matches, mismatches,
                insertions, deletions,
                unique_matches, unique_mismatches,
                unique_insertions, unique_deletions,
                longest_match, shortest_match, avg_match_length,
                longest_mismatch, shortest_mismatch, avg_mismatch_length,
                longest_insertion, shortest_insertion, avg_insertion_length,
                longest_deletion, shortest_deletion, avg_deletion_length) = compute_counts(&cigar_regex, cigar_str);


            // Print each metric along with the new columns, separated by tabs
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}",
                query_name, query_start, query_end, strand, target_name, target_start, target_end,
                matches, mismatches,
                insertions, deletions,
                unique_matches, unique_mismatches,
                unique_insertions, unique_deletions,
                longest_match, shortest_match, avg_match_length,
                longest_mismatch, shortest_mismatch, avg_mismatch_length,
                longest_insertion, shortest_insertion, avg_insertion_length,
                longest_deletion, shortest_deletion, avg_deletion_length
            );
        }
    }
    

    Ok(())
}

fn compute_counts(cigar_regex: &Regex, cigar_str: &str) -> (i64, i64, i64, i64, i64, i64, i64, i64, i64, i64, f64, i64, i64, f64, i64, i64, f64, i64, i64, f64) {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0;
    let mut deletions = 0;

    let mut unique_matches = 0;
    let mut unique_insertions = 0;
    let mut unique_deletions = 0;
    let mut unique_mismatches = 0;
    
    let mut longest_match = 0;
    let mut shortest_match = i64::MAX;
    let mut total_match_stretches = 0;

    let mut longest_mismatch = 0;
    let mut shortest_mismatch = i64::MAX;
    let mut total_mismatch_stretches = 0;

    let mut longest_insertion = 0;
    let mut shortest_insertion = i64::MAX;
    let mut total_insertion_length = 0;

    let mut longest_deletion = 0;
    let mut shortest_deletion = i64::MAX;
    let mut total_deletion_length = 0;

    for cap in cigar_regex.captures_iter(cigar_str) {
        let value: i64 = cap[1].parse().unwrap();
        match &cap[2] {
            "=" => {
                matches += value;
                if value > longest_match { longest_match = value; }
                if value < shortest_match { shortest_match = value; }
                unique_matches += 1;
                total_match_stretches += value;
            },
            "X" => {
                mismatches += value;
                if value > longest_mismatch { longest_mismatch = value; }
                if value < shortest_mismatch { shortest_mismatch = value; }
                unique_mismatches += 1;
                total_mismatch_stretches += value;
            },
            "I" => {
                insertions += value;
                unique_insertions += 1;
                if value > longest_insertion { longest_insertion = value; }
                if value < shortest_insertion { shortest_insertion = value; }
                total_insertion_length += value;
            },
            "D" => {
                deletions += value;
                unique_deletions += 1;
                if value > longest_deletion { longest_deletion = value; }
                if value < shortest_deletion { shortest_deletion = value; }
                total_deletion_length += value;
            },
            _ => {}
        }
    }

    // Adjust for cases where no =XID were found
    if shortest_match == i64::MAX { shortest_match = 0; }
    if shortest_mismatch == i64::MAX { shortest_mismatch = 0; }
    if shortest_insertion == i64::MAX { shortest_insertion = 0; }
    if shortest_deletion == i64::MAX { shortest_deletion = 0; }

    let avg_match_length = if unique_matches > 0 { total_match_stretches as f64 / unique_matches as f64 } else { 0f64 };
    let avg_mismatch_length = if unique_mismatches > 0 { total_mismatch_stretches as f64 / unique_mismatches as f64 } else { 0f64 };
    let avg_insertion_length = if unique_insertions > 0 { total_insertion_length as f64 / unique_insertions as f64 } else { 0f64 };
    let avg_deletion_length = if unique_deletions > 0 { total_deletion_length as f64 / unique_deletions as f64 } else { 0f64 };

    (
        matches, mismatches,
        insertions, deletions,
        unique_matches, unique_mismatches,
        unique_insertions, unique_deletions,
        longest_match, shortest_match, avg_match_length,
        longest_mismatch, shortest_mismatch, avg_mismatch_length,
        longest_insertion, shortest_insertion, avg_insertion_length,
        longest_deletion, shortest_deletion, avg_deletion_length
    )
}
