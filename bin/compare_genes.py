#!/usr/bin/env python3
"""Compare the Sec sites of Selenoprofiles predictions with those found by miniprot, for one genome.

Sites match when contig, strand and codon start are identical. miniprot sites are restricted to
the families searched by Selenoprofiles. SECIS hits found by cmsearch in the 3' flank of a
prediction are reported for that prediction; they are evidence, not a filter.
"""

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path

GENE_ID = re.compile(r'gene_id "([^"]+)"')


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def expand_profiles(profiles, sets_file):
    sets = {}
    with open(sets_file) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sets[row["set"]] = row["members"].split(",")

    def expand(name, seen):
        if name in sets and name not in seen:
            return {f for member in sets[name] for f in expand(member, seen | {name})}
        return {name}

    return {f for name in profiles.split(",") if name for f in expand(name.strip(), set())}


def prediction_parts(gene_id):
    """SELENOO.7.selenocysteine -> (SELENOO, selenocysteine)."""
    family, _index, label = gene_id.rsplit(".", 2)
    return family, label


def read_selenoprofiles(gtf):
    predictions, sites = {}, {}
    with open_text(gtf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            gene_id = GENE_ID.search(fields[8]).group(1)
            if fields[2] == "Selenocysteine":
                gene_id = gene_id.split(":", 1)[1]
                key = (fields[0], int(fields[3]), fields[6])
                sites[key] = gene_id
            predictions.setdefault(gene_id, prediction_parts(gene_id))
    return predictions, sites


def read_miniprot(tsv, families):
    sites = {}
    with open(tsv) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["family"] not in families:
                continue
            key = (row["contig"], int(row["start"]), row["strand"])
            previous = sites.get(key)
            if previous is None or int(row["score"]) > int(previous["score"]):
                sites[key] = row
    return sites


def read_secis(paths):
    """Best SECIS hit per prediction from cmsearch --tblout files of 3' flanks: (score, model, E-value, included).

    Hits below the inclusion threshold ('?') are kept, so that a lower reporting threshold (-T) shows weak
    candidates; only included hits ('!') count as SECIS support in the summary.
    """
    best = {}
    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < 17 or fields[9] != "+":
                    continue
                gene_id = fields[0].split(":", 1)[1] if ":" in fields[0] else fields[0]
                score = float(fields[14])
                if score > best.get(gene_id, (float("-inf"),))[0]:
                    best[gene_id] = (score, fields[2], fields[15], fields[16] == "!")
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--selenoprofiles-gtf", required=True)
    parser.add_argument("--miniprot-sites", help="Output of sec_check.py")
    parser.add_argument("--secis", nargs="*", default=[], help="cmsearch --tblout files, optionally gzip-compressed")
    parser.add_argument("--profiles", required=True, help="Profiles searched, as given to Selenoprofiles")
    parser.add_argument("--profile-sets", required=True, help="TSV of Selenoprofiles profile sets")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    families = expand_profiles(args.profiles, args.profile_sets)
    predictions, sp_sites = read_selenoprofiles(args.selenoprofiles_gtf)
    mp_sites = read_miniprot(args.miniprot_sites, families) if args.miniprot_sites else {}
    secis = read_secis(args.secis)

    layer = open(f"{args.prefix}.sec_codons.gff3", "w")
    layer.write("##gff-version 3\n# Sec codons predicted by Selenoprofiles (label selenocysteine or other) and/or miniprot;\n"
                "# a layer for gene annotation pipelines (Tico et al. 2026, Fig 5). Coordinates are 1-based, inclusive.\n")
    with open(f"{args.prefix}.sec_sites.tsv", "w") as out:
        out.write("contig\tstart\tend\tstrand\tfamily\tselenoprofiles_id\tselenoprofiles_label\tminiprot_support\tminiprot_query\tsecis_score\tsecis_evalue\tsecis_model\tsupport\n")
        for key in sorted(set(sp_sites) | set(mp_sites)):
            contig, start, strand = key
            gene_id = sp_sites.get(key)
            mp = mp_sites.get(key)
            family, label = predictions[gene_id] if gene_id else (mp["family"], ".")
            secis_score, secis_model, secis_evalue, _ = secis.get(gene_id, (".", ".", ".", False)) if gene_id else (".", ".", ".", False)
            support = "both" if gene_id and mp else ("selenoprofiles" if gene_id else "miniprot")
            out.write("\t".join(map(str, [
                contig, start, start + 2, strand, family, gene_id or ".", label,
                mp["support"] if mp else ".", mp["query"] if mp else ".", secis_score, secis_evalue, secis_model, support,
            ])) + "\n")
            attrs = {"ID": f"sec_{contig}_{start}", "family": family, "support": support, "selenoprofiles_id": gene_id or ".",
                     "selenoprofiles_label": label, "miniprot_query": mp["query"] if mp else ".", "secis_bits": secis_score,
                     "secis_evalue": secis_evalue}
            layer.write(f"{contig}\tselenoproteins\tstop_codon_redefined_as_selenocysteine\t{start}\t{start + 2}\t.\t{strand}\t.\t"
                        + ";".join(f"{k}={str(v).replace(';', ',').replace('=', ' ')}" for k, v in attrs.items()) + "\n")
    layer.close()

    labels = Counter(label for _family, label in predictions.values())
    sec_sp = {k for k, g in sp_sites.items() if predictions[g][1] == "selenocysteine"}
    both = sec_sp & set(mp_sites)
    union = sec_sp | set(mp_sites)
    sec_genes = {g for g, (_f, label) in predictions.items() if label == "selenocysteine"}
    summary = {
        "sample": args.sample,
        "predictions": len(predictions),
        "selenocysteine": labels.get("selenocysteine", 0),
        "cysteine": labels.get("cysteine", 0),
        "pseudo": labels.get("pseudo", 0),
        "other_labels": sum(n for label, n in labels.items() if label not in ("selenocysteine", "cysteine", "pseudo")),
        "sec_families": len({predictions[g][0] for g in sec_genes}),
        "sec_sites_selenoprofiles": len(sec_sp),
        "sec_sites_miniprot": len(mp_sites) if args.miniprot_sites else "NA",
        "sec_sites_both": len(both) if args.miniprot_sites else "NA",
        "sec_sites_jaccard": (f"{len(both) / len(union):.3f}" if union else "NA") if args.miniprot_sites else "NA",
        "secis_supported": len(sec_genes & {g for g, v in secis.items() if v[3]}) if args.secis else "NA",
    }
    with open(f"{args.prefix}.genes_summary.tsv", "w") as out:
        out.write("\t".join(summary) + "\n" + "\t".join(map(str, summary.values())) + "\n")


if __name__ == "__main__":
    main()
