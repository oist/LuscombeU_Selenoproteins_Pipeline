#!/usr/bin/env python3
"""Check how a reference annotation handles each Sec codon predicted by Selenoprofiles.

Classes follow Tico et al. 2026 (PLoS Comput Biol 22:e1013885, Fig 1), at gene level (best transcript per
Sec codon) and transcript level. "Well annotated" is split by whether the annotation declares Sec at that
codon (transl_except) or not (the protein carries X or another residue there).

For every Selenocysteine feature in the Selenoprofiles GTF, annotated CDSs on the same strand that overlap
the prediction are classified:
  Well annotated  the Sec codon is inside the CDS, in frame (sec_residue: Sec declared, or X/other)
  Stop codon      the CDS ends exactly at the Sec codon
  Skipped         annotated CDS on both sides of the Sec codon, but not over it
  Downstream      the annotated CDS lies entirely 3' of the Sec codon (the protein starts after it)
  Upstream        the annotated CDS lies entirely 5' of the Sec codon (and does not end at it)
  Out of frame    inside the CDS, in another reading frame
  Spliced         a CDS boundary falls inside the Sec codon
  Absent          no overlapping CDS; the note says when the codon lies in the UTR of another transcript
"""
import argparse
import gzip
import re
from collections import Counter, defaultdict
from urllib.parse import unquote

# classes of Tico et al. 2026 (PLoS Comput Biol 22:e1013885, Fig 1), best first for the gene-level call
ORDER = ["Well annotated", "Stop codon", "Skipped", "Downstream", "Upstream", "Out of frame", "Spliced", "Absent"]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attrs(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_gff(path):
    tx, exons, cds = {}, defaultdict(list), defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] in ("mRNA", "transcript"):
                tx[a["ID"]] = (f[0], f[6], int(f[3]), int(f[4]), a.get("gene", a.get("Parent", "")), a.get("product", ""))
            elif f[2] == "exon":
                for p in a.get("Parent", "").split(","):
                    exons[p].append((int(f[3]), int(f[4])))
            elif f[2] == "CDS":
                for p in a.get("Parent", "").split(","):
                    cds[p].append((int(f[3]), int(f[4]), a.get("protein_id", a.get("ID", "")), a.get("transl_except", "")))
    return tx, exons, cds


def classify(sec, strand, c):
    """sec: (start, end) of the Sec codon. Returns (class, residue) for one annotated CDS.

    residue (Well annotated only): 'Sec' when a transl_except declares it, otherwise 'X/other', i.e. the
    annotation contains the codon but its protein does not carry Sec there.
    """
    parts = sorted((a, b) for a, b, _, _ in c)
    lo, hi = parts[0][0], parts[-1][1]
    codon = set(range(sec[0], sec[1] + 1))
    covered = {p for p in codon if any(a <= p <= b for a, b in parts)}
    if covered and covered != codon:
        return "Spliced", "."
    first = sec[0] if strand == "+" else sec[1]
    if not covered:
        if sec[1] < lo:
            return ("Upstream" if strand == "-" else "Downstream"), "."
        if sec[0] > hi:
            return ("Downstream" if strand == "-" else "Upstream"), "."
        return "Skipped", "."
    # inside: frame relative to the CDS start
    order = parts if strand == "+" else parts[::-1]
    off = 0
    for a, b in order:
        if a <= first <= b:
            idx = off + (first - a if strand == "+" else b - first)
            break
        off += b - a + 1
    total = sum(b - a + 1 for a, b in parts)
    if idx % 3:
        return "Out of frame", "."
    if idx == total - 3:
        return "Stop codon", "."
    te = unquote(",".join(x[3] for x in c))
    m = re.findall(r"pos:(?:complement\()?(\d+)\.\.(\d+)\)?,aa:Sec", te)
    return "Well annotated", ("Sec" if any(int(a) <= first <= int(b) for a, b in m) else "X/other")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--selenoprofiles-gtf", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    tx, exons, cds = read_gff(args.annotation)
    by_contig = defaultdict(list)
    for tid, t in tx.items():
        if tid in cds:
            by_contig[(t[0], t[1])].append(tid)
    preds, sites = {}, []
    with open_text(args.selenoprofiles_gtf) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            if f[2] == "CDS":
                p = preds.setdefault(gid, [f[0], f[6], int(f[3]), int(f[4])])
                p[2], p[3] = min(p[2], int(f[3])), max(p[3], int(f[4]))
            elif f[2] == "Selenocysteine":
                sites.append((gid.split(":", 1)[1], f[0], int(f[3]), int(f[4]), f[6]))
    out = open(f"{args.prefix}.sec_annotation.tsv", "w")
    out.write("sample\tprediction\tlabel\tcontig\tsec_start\tsec_end\tstrand\tstatus\tsec_residue\ttranscript\tprotein\tproduct\tnote\n")
    tx_out = open(f"{args.prefix}.sec_annotation_transcripts.tsv", "w")
    tx_out.write("sample\tprediction\ttranscript\tprotein\tstatus\tsec_residue\n")
    gene_counts, tx_counts, residue = Counter(), Counter(), Counter()
    for gid, contig, s, e, strand in sites:
        label = gid.rsplit(".", 1)[1]
        _, _, plo, phi = preds.get(gid, (contig, strand, s, e))
        best = ("Absent", ".", ".", ".", ".")
        for tid in by_contig.get((contig, strand), []):
            t = tx[tid]
            clo, chi = min(x[0] for x in cds[tid]), max(x[1] for x in cds[tid])
            if chi < plo or clo > phi:   # the annotated CDS must overlap the prediction
                continue
            st, res = classify((s, e), strand, cds[tid])
            tx_counts[st] += 1
            tx_out.write(f"{args.sample}\t{gid}\t{tid}\t{cds[tid][0][2]}\t{st}\t{res}\n")
            if ORDER.index(st) < ORDER.index(best[0]) or (st == best[0] == "Well annotated" and res == "Sec" and best[1] != "Sec"):
                best = (st, res, tid, cds[tid][0][2], t[5])
        note = "."
        if best[0] == "Absent":
            tx_counts["Absent"] += 1
            host = [tid for tid in by_contig.get((contig, strand), []) if tx[tid][2] <= s and e <= tx[tid][3]]
            if host:
                note = "inside transcript " + ",".join(f"{h} ({cds[h][0][2]})" for h in host[:2])
        gene_counts[best[0]] += 1
        if best[0] == "Well annotated":
            residue[best[1]] += 1
        out.write(f"{args.sample}\t{gid}\t{label}\t{contig}\t{s}\t{e}\t{strand}\t" + "\t".join(best) + f"\t{note}\n")
    with open(f"{args.prefix}.sec_annotation_summary.tsv", "w") as summ:
        cols = ORDER + ["Well annotated, Sec declared", "Well annotated, no Sec declared"] + [f"transcripts: {k}" for k in ORDER]
        vals = [gene_counts[k] for k in ORDER] + [residue["Sec"], residue["X/other"]] + [tx_counts[k] for k in ORDER]
        summ.write("sample\t" + "\t".join(cols) + "\n" + args.sample + "\t" + "\t".join(map(str, vals)) + "\n")


if __name__ == "__main__":
    main()
