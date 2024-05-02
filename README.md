# pafgnostic

## Overview
`pafgnostic` supports both uncompressed and gzipped PAF (Pairwise mApping Format) files, and it outputs a comprehensive set of metrics for each alignment, including the number of matches, mismatches, insertions, deletions, unique events, and various lengths and averages for these events.

## Installation

`pafgnostic` is written in Rust, so you can install it using `cargo`:

```shell
git clone https://github.com/AndreaGuarracino/pafgnostic
cd pafgnostic
cargo install --force --path .
```

## Usage
Run `pafgnostic` with a PAF file (compressed or uncompressed) as input:

```shell
pafgnostic --paf input.paf
```

## Output
`pafgnostic` outputs a tab-separated table directly to the console with the following columns:

- Query name, start, end, strand
- Target name, start, end
- Number of matches, mismatches, insertions, deletions
- Unique matches, mismatches, insertions, deletions
- Longest, shortest, and average lengths of matches, mismatches, insertions, deletions
