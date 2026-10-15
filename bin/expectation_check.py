#!/usr/bin/env python3
"""Compare the Sec sites found in a genome with an expected selenoproteome.

The expectation file uses the layout of selenoprotein_profiles Expectation_table.csv (one row per lineage, one
column per subfamily; lines starting with # are comments). Subfamily columns are summed per Selenoprofiles family
(GPX1..GPX6 -> GPX, DIO1..3 -> DIO, TXNRD1..3 -> TXNRD, SELENOW1/2 and SELENOV -> SELENOW), because the pipeline
does not assign subfamilies. Observed counts are the Sec sites of the gene comparison table that are labelled
selenocysteine by Selenoprofiles or found by miniprot only.
"""
import argparse
import csv
import re
from collections import Counter


def family(col):
    if col == "SELENOV":
        return "SELENOW"
    if col == "SEPHS2":
        return col
    return re.sub(r"\d+$", "", col)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sites", required=True, help="<sample>.sec_sites.tsv from compare_genes.py")
    ap.add_argument("--expectation", required=True)
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--copies", type=int, default=1, help="Expected copies per gene, e.g. 2 for a diploid assembly")
    ap.add_argument("--hemizygous", default="", help="Comma-separated families counted once whatever --copies (sex-specific regions)")
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    rows = [line for line in open(args.expectation) if line.strip() and not line.startswith("#")]
    table = list(csv.DictReader(rows))
    row = next((r for r in table if r["Lineage"] == args.lineage), None)
    if row is None:
        raise SystemExit(f"Lineage {args.lineage} not in {args.expectation}: {[r['Lineage'] for r in table]}")
    expected = Counter()
    for col, n in row.items():
        if col != "Lineage" and n not in ("", None):
            fam = family(col)
            expected[fam] += int(n) * (1 if fam in args.hemizygous.split(",") else args.copies)
    observed = Counter()
    with open(args.sites) as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["selenoprofiles_label"] == "selenocysteine" or r["support"] == "miniprot":
                observed[r["family"]] += 1
    status = Counter()
    with open(f"{args.prefix}.expectation.tsv", "w") as out:
        out.write("sample\tfamily\texpected\tobserved\tstatus\n")
        for fam in sorted(set(expected) | set(observed)):
            e, o = expected[fam], observed[fam]
            if e == 0 and o == 0:
                continue
            st = "as expected" if e == o else ("missing" if o < e else "extra")
            status[st] += abs(e - o) if st != "as expected" else 1
            out.write(f"{args.sample}\t{fam}\t{e}\t{o}\t{st}\n")
    missing = [f"{f} ({expected[f] - observed[f]})" for f in sorted(expected) if observed[f] < expected[f]]
    extra = [f"{f} ({observed[f] - expected[f]})" for f in sorted(observed) if observed[f] > expected[f]]
    with open(f"{args.prefix}.expectation_summary.tsv", "w") as out:
        out.write("sample\texpected\tobserved\tmissing\textra\n")
        out.write(f"{args.sample}\t{sum(expected.values())}\t{sum(observed.values())}\t{', '.join(missing) or 'none'}\t{', '.join(extra) or 'none'}\n")


if __name__ == "__main__":
    main()
