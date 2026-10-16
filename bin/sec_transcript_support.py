#!/usr/bin/env python3
"""Transcript support of Sec codons.

Reads the Sec codons of a genome (the sec_sites.tsv table of COMPARE_GENES) and spliced alignments of transcripts or
ESTs to that genome (minimap2 PAF with the cg:Z CIGAR tag). A transcript supports a Sec codon when the three
nucleotides of the codon lie inside one aligned exon of its primary alignment. For each Sec codon, writes the number
of supporting transcripts and how far the longest of them extends 3' of the codon (on the strand of the gene), which
shows whether transcripts reach the SECIS region.
"""

import argparse
import bisect
import csv
import re
from collections import defaultdict

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def exons(tstart, cigar):
    """1-based inclusive target intervals of the aligned exons (CIGAR operation N separates exons)."""
    blocks, pos, start = [], tstart, tstart
    for length, op in CIGAR.findall(cigar):
        length = int(length)
        if op in "M=XD":
            pos += length
        elif op == "N":
            if pos > start:
                blocks.append((start + 1, pos))
            pos += length
            start = pos
    if pos > start:
        blocks.append((start + 1, pos))
    return blocks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="sec_sites.tsv from COMPARE_GENES")
    ap.add_argument("--paf", required=True, help="minimap2 alignments (PAF with cg:Z tags)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.sites) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        sites = list(reader)
        fields = reader.fieldnames or []
    by_contig = defaultdict(list)
    for i, s in enumerate(sites):
        by_contig[s["contig"]].append((int(s["start"]), i))
    for v in by_contig.values():
        v.sort()
    support = [0] * len(sites)
    extent = [0] * len(sites)

    with open(args.paf) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[5] not in by_contig:
                continue
            tags = dict(t.split(":", 2)[0::2] for t in f[12:] if t.count(":") >= 2)
            if tags.get("tp", "P") != "P" or "cg" not in tags:
                continue
            tstart, tend = int(f[7]), int(f[8])          # 0-based, end exclusive
            starts = by_contig[f[5]]
            lo = bisect.bisect_left(starts, (tstart + 1, -1))
            hi = bisect.bisect_right(starts, (tend - 2, len(sites)))
            if lo >= hi:
                continue
            blocks = exons(tstart, tags["cg"])
            for start, i in starts[lo:hi]:
                end = start + 2
                if any(a <= start and end <= b for a, b in blocks):
                    support[i] += 1
                    ext = tend - end if sites[i]["strand"] == "+" else (start - 1) - tstart
                    extent[i] = max(extent[i], ext)

    with open(args.output, "w") as out:
        keep = [c for c in ("contig", "start", "end", "strand", "family", "support") if c in fields]
        out.write("\t".join(keep + ["transcripts", "max_3prime_extent"]) + "\n")
        for i, s in enumerate(sites):
            out.write("\t".join([s[c] for c in keep] + [str(support[i]), str(extent[i])]) + "\n")


if __name__ == "__main__":
    main()
