#!/usr/bin/env python3
"""MultiQC table of the selenocysteine machinery found in each genome.

Columns: tRNA-Sec (loci supported by two or more tools, from the per-genome tRNA-Sec summaries) and one column per
gene of the Selenoprofiles machinery profile set (SEPHS2, PSTK, SecS, eEFSec, SBP2, and the associated factor SECp43), read from the
Selenoprofiles GTF files (<sample>.selenoprofiles.gtf). Values:
  Sec        SEPHS2 predicted with a selenocysteine codon
  found      predicted (label homologue, or another residue at the selenocysteine position)
  partial    predicted only as pseudo (in-frame stop or frameshift) or unaligned at the selenocysteine position
  -          no prediction
Absences depend on assembly completeness: a gene broken across contigs may be missed. Presences mean homology to the
profile, not orthology: a related paralog (for example KTI12 for PSTK) can be predicted by a profile.
"""
import argparse, collections, csv, os, re

GENES = [("SEPHS2", "SEPHS2"), ("PSTK", "pstk"), ("SEPSECS", "SecS"), ("EEFSEC", "eEFsec"), ("SECISBP2", "SBP2"),
         ("TRNAU1AP", "secp43")]


def labels(gtf):
    c = collections.defaultdict(collections.Counter)
    seen = set()
    for l in open(gtf):
        m = re.search(r'gene_id "([^":]+)"', l)  # "Sec1:<id>" lines mark selenocysteine positions
        if not m or m.group(1) in seen: continue
        seen.add(m.group(1))
        c[m.group(1).split(".")[0]][m.group(1).rsplit(".", 1)[1]] += 1
    return c


def status(fam, c):
    if not c: return "-"
    if fam == "SEPHS2" and c.get("selenocysteine"): return "Sec"
    if any(k not in ("pseudo", "unaligned") for k in c): return "found"
    return "partial"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="sample names, one per line, in report order")
    ap.add_argument("--gtf", nargs="+", required=True)
    ap.add_argument("--trnasec", nargs="+", required=True)
    ap.add_argument("--out", default="machinery_mqc.tsv")
    a = ap.parse_args()
    samples = [l.strip() for l in open(a.samples) if l.strip()]
    gtfs = {os.path.basename(p)[:-len(".selenoprofiles.gtf")]: p for p in a.gtf}
    trna = {}
    for p in a.trnasec:
        for r in csv.DictReader(open(p), delimiter="\t"): trna[r["sample"]] = r.get("loci_supported", "")
    with open(a.out, "w") as o:
        o.write("# id: 'machinery'\n# section_name: 'Selenocysteine machinery'\n"
                "# description: 'tRNA-Sec loci supported by two or more tools, and the Selenoprofiles predictions of "
                "the machinery genes and of the associated factor SECp43 (TRNAU1AP), which is not required for "
                "selenoprotein synthesis in all cells: Sec, SEPHS2 with a selenocysteine codon; found; partial, only pseudo or "
                "unaligned predictions; -, none. Found means homologous to the profile: orthology to the named gene, as "
                "opposed to a related paralog, is not tested. A gene broken by the assembly may be missed.'\n# plot_type: 'table'\n")
        o.write("Sample\ttRNA-Sec\t" + "\t".join(g for g, _ in GENES) + "\n")
        for s in samples:
            c = labels(gtfs[s]) if s in gtfs else {}
            o.write(f"{s}\t{trna.get(s, '')}\t" + "\t".join(status(f, c.get(f, {})) for _, f in GENES) + "\n")


if __name__ == "__main__":
    main()
