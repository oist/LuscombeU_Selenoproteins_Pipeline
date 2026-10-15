#!/usr/bin/env python3
"""Merge tRNA-Sec predictions of several tools for one genome.

Predictions on the same contig and strand that overlap by at least one nucleotide form one locus.
Each locus records which tools support it. A locus is supported when two or more tools find it
and, if Secmarker reports it, its discriminator base is G73, as expected for functional
tRNA-Sec genes.
"""

import argparse
import gzip
import re


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_secmarker(path):
    hits = []
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attributes = dict(re.findall(r"(\w+):(\S+)", fields[8]))
            hits.append({
                "contig": fields[0], "start": int(fields[3]), "end": int(fields[4]), "strand": fields[6],
                "score": float(attributes.get("infernal_score", "nan")),
                "detail": f"model={attributes.get('model', '.')};discriminator={attributes.get('discr_base', '.')};anticodon={attributes.get('anticodon', '.')}",
                "discriminator": attributes.get("discr_base", "."),
            })
    return hits


def read_cmsearch(path):
    hits = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 17 or fields[16] != "!":
                continue
            start, end = sorted((int(fields[7]), int(fields[8])))
            hits.append({"contig": fields[0], "start": start, "end": end, "strand": fields[9], "score": float(fields[14]), "detail": f"model={fields[2]}", "discriminator": "."})
    return hits


def read_trnascan(path):
    hits = []
    with open_text(path) as handle:
        for line in handle:
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if len(fields) < 9 or not fields[1].isdigit() or fields[4] != "SeC":
                continue
            begin, end = int(fields[2]), int(fields[3])
            strand = "+" if begin <= end else "-"
            start, end = sorted((begin, end))
            hits.append({"contig": fields[0], "start": start, "end": end, "strand": strand, "score": float(fields[8]), "detail": f"anticodon={fields[5]}", "discriminator": "."})
    return hits


def cluster(tool_hits):
    loci = []
    events = sorted(((h["contig"], h["strand"], h["start"], h["end"], tool, h) for tool, hits in tool_hits.items() for h in hits))
    for contig, strand, start, end, tool, hit in events:
        last = loci[-1] if loci else None
        if last and last["contig"] == contig and last["strand"] == strand and start <= last["end"]:
            last["end"] = max(last["end"], end)
        else:
            last = {"contig": contig, "strand": strand, "start": start, "end": end, "tools": {}}
            loci.append(last)
        if tool not in last["tools"] or hit["score"] > last["tools"][tool]["score"]:
            last["tools"][tool] = hit
    return loci


def is_supported(locus):
    """Two or more tools, and the G73 discriminator base when Secmarker reports the locus."""
    secmarker = locus["tools"].get("secmarker")
    return len(locus["tools"]) >= 2 and (secmarker is None or secmarker["discriminator"] == "G")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--secmarker", help="Secmarker GFF")
    parser.add_argument("--cmsearch", help="cmsearch --tblout file of the tRNA-Sec model, optionally gzip-compressed")
    parser.add_argument("--trnascanse", help="tRNAscan-SE tabular output")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    readers = {"secmarker": read_secmarker, "cmsearch": read_cmsearch, "trnascanse": read_trnascan}
    tool_hits = {tool: readers[tool](getattr(args, tool)) for tool in readers if getattr(args, tool)}
    loci = cluster(tool_hits)

    with open(f"{args.prefix}.trnasec_loci.tsv", "w") as out:
        out.write("contig\tstart\tend\tstrand\tn_tools\t" + "\t".join(f"{t}_score\t{t}_detail" for t in tool_hits) + "\n")
        for locus in loci:
            cells = []
            for tool in tool_hits:
                hit = locus["tools"].get(tool)
                cells += [f"{hit['score']:g}", hit["detail"]] if hit else [".", "."]
            out.write("\t".join(map(str, [locus["contig"], locus["start"], locus["end"], locus["strand"], len(locus["tools"])] + cells)) + "\n")

    summary = {"sample": args.sample}
    for tool in readers:
        summary[f"{tool}_hits"] = len(tool_hits[tool]) if tool in tool_hits else "NA"
    summary["loci"] = len(loci)
    summary["loci_multiple_tools"] = sum(1 for locus in loci if len(locus["tools"]) >= 2)
    summary["loci_supported"] = sum(1 for locus in loci if is_supported(locus))
    summary["secmarker_g73"] = sum(1 for h in tool_hits.get("secmarker", []) if h["discriminator"] == "G") if "secmarker" in tool_hits else "NA"
    with open(f"{args.prefix}.trnasec_summary.tsv", "w") as out:
        out.write("\t".join(summary) + "\n" + "\t".join(map(str, summary.values())) + "\n")


if __name__ == "__main__":
    main()
