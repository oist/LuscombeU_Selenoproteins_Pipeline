#!/usr/bin/env python3
"""Screen annotated transcripts for TGA codons that may encode selenocysteine.

For each transcript with a CDS (GFF3 with mRNA or transcript, exon and CDS features):
  start_extension  an in-frame TGA upstream of the CDS start, with an in-frame ATG further upstream
                   and no TAA/TAG in between (the protein may start before a Sec codon)
  upstream_tga     an in-frame TGA upstream of the CDS start without such an ATG (the 5' exons may be
                   missing or wrong)
  stop_readthrough the CDS ends with TGA and the frame continues to a TAA/TAG in the transcript
  internal_tga     a TGA inside the CDS, translated as something other than Sec
  sec_declared     a transl_except=...aa:Sec attribute on the CDS
Writes <prefix>.readthrough.tsv (one row per transcript and event; event "none" when there is no event) and <prefix>.3utr.fa with the
3' UTR of every transcript after the last in-frame stop (TAA/TAG), for SECIS searches. For annotations without
UTRs, --extend adds genomic sequence at both ends of each transcript (no splicing is assumed there).
"""
import argparse
import gzip
import re
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG"}


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = "".join(buf).upper()
    return seqs


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
                tx[a["ID"]] = {"contig": f[0], "strand": f[6], "gene": a.get("gene", a.get("Parent", "")), "product": a.get("product", "")}
            elif f[2] == "exon" and "Parent" in a:
                for p in a["Parent"].split(","):
                    exons[p].append((int(f[3]), int(f[4])))
            elif f[2] == "CDS" and "Parent" in a:
                for p in a["Parent"].split(","):
                    cds[p].append((int(f[3]), int(f[4]), a.get("protein_id", a.get("ID", "")), a.get("transl_except", "")))
    return tx, exons, cds


def spliced(seq, parts, strand):
    if strand == "+":
        return "".join(seq[a - 1:b] for a, b in sorted(parts))
    return "".join(seq[a - 1:b][::-1].translate(COMP) for a, b in sorted(parts, reverse=True))


def tx_index(parts, strand, pos):
    """Index of genomic position pos in the spliced transcript (5' to 3')."""
    off = 0
    order = sorted(parts) if strand == "+" else sorted(parts, reverse=True)
    for a, b in order:
        if a <= pos <= b:
            return off + (pos - a if strand == "+" else b - pos)
        off += b - a + 1
    return None


def genomic(parts, strand, idx):
    off = 0
    order = sorted(parts) if strand == "+" else sorted(parts, reverse=True)
    for a, b in order:
        n = b - a + 1
        if idx < off + n:
            return a + (idx - off) if strand == "+" else b - (idx - off)
        off += n
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--max-utr", type=int, default=1000, help="Longest 3' UTR written for SECIS searches")
    ap.add_argument("--extend", type=int, default=0,
                    help="Extend each transcript by this many genomic nt at both ends, for annotations without UTRs (e.g. AUGUSTUS)")
    args = ap.parse_args()

    tx, exons, cds = read_gff(args.gff)
    genome = read_fasta(args.genome, {t["contig"] for t in tx.values()})
    out = open(f"{args.prefix}.readthrough.tsv", "w")
    out.write("transcript\tgene\tprotein\tproduct\tcontig\tstrand\tevent\textra_codons\ttga_codons\ttga_positions\tutr_length\n")
    utr = open(f"{args.prefix}.3utr.fa", "w")
    for tid, t in tx.items():
        if tid not in cds or t["contig"] not in genome:
            continue
        ex = sorted(exons.get(tid) or [(a, b) for a, b, _, _ in cds[tid]])
        strand, seq = t["strand"], genome[t["contig"]]
        if args.extend:
            ex[0] = (max(1, ex[0][0] - args.extend), ex[0][1])
            ex[-1] = (ex[-1][0], min(len(seq), ex[-1][1] + args.extend))
        mrna = spliced(seq, ex, strand)
        c = cds[tid]
        pid = c[0][2]
        te = ",".join(sorted({x[3] for x in c if x[3]}))
        lo, hi = min(a for a, b, _, _ in c), max(b for a, b, _, _ in c)
        s_idx = tx_index(ex, strand, lo if strand == "+" else hi)
        e_idx = tx_index(ex, strand, hi if strand == "+" else lo)
        if s_idx is None or e_idx is None:
            continue
        rows = []
        # inside the CDS (only when the CDS length is a multiple of 3; corrected models may not be)
        internal = [i for i in range(s_idx, e_idx - 2, 3) if mrna[i:i + 3] == "TGA"] if (e_idx - s_idx + 1) % 3 == 0 else []
        if "Sec" in te:
            rows.append(("sec_declared", 0, len(internal), internal))
        elif internal:
            rows.append(("internal_tga", 0, len(internal), internal))
        # upstream of the start codon
        i, tgas, best = s_idx - 3, [], None
        while i >= 0:
            cod = mrna[i:i + 3]
            if cod in STOPS:
                break
            if cod == "TGA":
                tgas.append(i)
            elif cod == "ATG" and tgas:
                best = (i, list(tgas))
            i -= 3
        if best:
            rows.append(("start_extension", (s_idx - best[0]) // 3, len(best[1]), best[1]))
        elif tgas:
            rows.append(("upstream_tga", (s_idx - tgas[-1]) // 3, len(tgas), tgas))
        # through the stop codon
        end = e_idx + 1
        if mrna[e_idx - 2:e_idx + 1] == "TGA":
            j, tgas = e_idx + 1, [e_idx - 2]
            while j + 3 <= len(mrna):
                cod = mrna[j:j + 3]
                if cod in STOPS:
                    break
                if cod == "TGA":
                    tgas.append(j)
                j += 3
            if j + 3 <= len(mrna):
                rows.append(("stop_readthrough", (j - e_idx - 1) // 3, len(tgas), tgas))
                end = j + 3
        u = mrna[end:end + args.max_utr]
        utr.write(f">{tid}\n{u}\n")
        if not rows:
            rows.append(("none", 0, 0, []))
        for ev, extra, n, pos in rows:
            gpos = ",".join(str(genomic(ex, strand, p)) for p in pos[:5])
            out.write(f"{tid}\t{t['gene']}\t{pid}\t{t['product']}\t{t['contig']}\t{strand}\t{ev}\t{extra}\t{n}\t{gpos}\t{len(mrna) - end}\n")
    out.close()
    utr.close()


if __name__ == "__main__":
    main()
