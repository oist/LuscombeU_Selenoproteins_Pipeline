#!/usr/bin/env python3
"""Rank candidate selenoprotein transcripts: TGA events from readthrough_screen.py combined with the best
SECIS hit in each transcript's 3' UTR (cmsearch --tblout, one or more models).

Writes <prefix>.candidates.tsv: transcripts with a TGA event and a SECIS hit of at least --min-bits in any
model, plus transcripts with a SECIS hit of at least --secis-only-bits and no TGA event (the UTR may belong to
an unannotated gene, or the TGA may lie outside the annotated exons).
"""
import argparse
import csv
import gzip
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def best_hits(path):
    best = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 17 or f[9] != "+":
                continue
            sc = float(f[14])
            if sc > best.get(f[0], (float("-inf"),))[0]:
                best[f[0]] = (sc, f[15], int(f[7]), int(f[8]))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True)
    ap.add_argument("--secis", nargs="+", required=True, help="cmsearch --tblout files, one per model (name taken from the file name)")
    ap.add_argument("--min-bits", type=float, default=15)
    ap.add_argument("--secis-only-bits", type=float, default=30)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    models = {p.split("/")[-1].split(".tbl")[0].split(".")[-1]: best_hits(p) for p in args.secis}
    events = defaultdict(list)
    info = {}
    with open(args.events) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            info[row["transcript"]] = row
            if row["event"] != "none":
                events[row["transcript"]].append(row)
    names = sorted(models)
    out = open(f"{args.prefix}.candidates.tsv", "w")
    out.write("rank_score\ttranscript\tgene\tprotein\tproduct\tcontig\tstrand\tevents\ttga_positions\t"
              + "\t".join(f"secis_{m}_bits\tsecis_{m}_evalue\tsecis_{m}_utr_pos" for m in names) + "\n")
    rows = []
    for tid in set(events) | {t for m in models.values() for t in m}:
        scores = [models[m].get(tid) for m in names]
        top = max((s[0] for s in scores if s), default=None)
        if top is None:
            continue
        ev = events.get(tid, [])
        if (ev and top < args.min_bits) or (not ev and top < args.secis_only_bits):
            continue
        r = info.get(tid, {})
        evs = ";".join(f"{e['event']}({e['extra_codons']} codons)" for e in ev) or "secis_only"
        pos = ";".join(e["tga_positions"] for e in ev)
        cells = []
        for s in scores:
            cells += [f"{s[0]:.1f}", s[1], f"{s[2]}-{s[3]}"] if s else [".", ".", "."]
        rows.append((top, [f"{top:.1f}", tid, r.get("gene", "."), r.get("protein", "."), r.get("product", "."), r.get("contig", "."),
                           r.get("strand", "."), evs, pos or "."] + cells))
    for top, cells in sorted(rows, key=lambda x: -x[0]):
        out.write("\t".join(cells) + "\n")


if __name__ == "__main__":
    main()
