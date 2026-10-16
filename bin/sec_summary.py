#!/usr/bin/env python3
"""Aggregate per-genome summaries into MultiQC custom-content tables and a Sec trait call.

Sec trait call per genome:
  present       at least one prediction labelled selenocysteine and one supported tRNA-Sec
                locus (two or more tools; G73 discriminator base when Secmarker reports it)
  absent        neither
  inconsistent  only one of the two
Absence calls depend on assembly completeness.
"""

import argparse
import csv
import statistics
from collections import Counter
from pathlib import Path


def read_rows(paths):
    rows = {}
    for path in paths:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows[row["sample"]] = row
    return rows


def write_table(path, section_id, name, description, headers, rows, plot_type="table"):
    with open(path, "w") as out:
        q = lambda text: text.replace("'", "''")  # YAML single-quoted strings
        out.write(f"# id: '{section_id}'\n# section_name: '{q(name)}'\n# description: '{q(description)}'\n# plot_type: '{plot_type}'\n")
        out.write("\t".join(["Sample"] + headers) + "\n")
        for sample, values in rows:
            out.write("\t".join([sample] + [str(values.get(h, "NA")) for h in headers]) + "\n")


def trait(genes, trnasec):
    if genes is None or trnasec is None:
        return "NA"
    has_genes = int(genes["selenocysteine"]) > 0
    has_trna = int(trnasec["loci_supported"]) > 0
    if has_genes and has_trna:
        return "present"
    if not has_genes and not has_trna:
        return "absent"
    return "inconsistent"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", required=True, help="File listing samples in samplesheet order")
    parser.add_argument("--genes", nargs="*", default=[])
    parser.add_argument("--trnasec", nargs="*", default=[])
    parser.add_argument("--assess", nargs="*", default=[], help="selenoprofiles assess per-gene tables named <sample>.genes.tsv")
    parser.add_argument("--sec-annotation", nargs="*", default=[], help="sec_annotation_check.py summaries")
    parser.add_argument("--candidates", nargs="*", default=[], help="rank_readthrough.py tables named <sample>.candidates.tsv")
    parser.add_argument("--expectation", nargs="*", default=[], help="expectation_check.py summaries")
    parser.add_argument("--transcripts", nargs="*", default=[], help="sec_transcript_support.py tables named <sample>.sec_transcripts.tsv")
    args = parser.parse_args()

    samples = [line.strip() for line in open(args.samples) if line.strip()]
    genes = read_rows(args.genes)
    trnasec = read_rows(args.trnasec)
    assess = {}
    for path in args.assess:
        sample = Path(path).name[: -len(".genes.tsv")]
        with open(path) as handle:
            assess[sample] = Counter(row["Type_annotation"] for row in csv.DictReader(handle, delimiter="\t"))

    gene_headers = ["predictions", "selenocysteine", "cysteine", "pseudo", "other_labels", "sec_families"]
    write_table("selenoprotein_genes_mqc.tsv", "selenoprotein_genes", "Selenoprotein genes",
                "Selenoprofiles predictions per label, and number of families with at least one selenocysteine-labelled prediction.",
                gene_headers, [(s, genes[s]) for s in samples if s in genes])

    agreement_headers = ["sec_sites_selenoprofiles", "sec_sites_miniprot", "sec_sites_both", "sec_sites_jaccard", "secis_supported"]
    write_table("sec_site_agreement_mqc.tsv", "sec_site_agreement", "Sec site agreement",
                "Sec codons predicted by Selenoprofiles (selenocysteine label) and by miniprot, their intersection and Jaccard index, and Selenoprofiles selenoproteins with a SECIS hit in their 3&#39; flank.",
                agreement_headers, [(s, genes[s]) for s in samples if s in genes])

    trna_headers = ["secmarker_hits", "cmsearch_hits", "trnascanse_hits", "loci", "loci_multiple_tools", "loci_supported", "secmarker_g73"]
    write_table("trnasec_mqc.tsv", "trnasec", "tRNA-Sec",
                "tRNA-Sec predictions per tool, merged loci, loci found by two or more tools, supported loci (two or more tools and, when Secmarker reports the locus, the G73 discriminator base), and Secmarker hits with G73.",
                trna_headers, [(s, trnasec[s]) for s in samples if s in trnasec])

    if assess:
        statuses = sorted({status for counts in assess.values() for status in counts})
        write_table("annotation_assessment_mqc.tsv", "annotation_assessment", "Annotation assessment",
                    "Status of the reference annotation at each Selenoprofiles prediction, from selenoprofiles assess.",
                    statuses, [(s, assess[s]) for s in samples if s in assess])

    sec_annot = read_rows(args.sec_annotation)
    if sec_annot:
        heads = [h for h in next(iter(sec_annot.values())) if h != "sample"]
        write_table("sec_annotation_mqc.tsv", "sec_annotation", "Sec codons in the reference annotation",
                    "How the given annotation treats each predicted Sec codon: Sec declared, translated as X or stop, outside the CDS, or no model.",
                    heads, [(s, sec_annot[s]) for s in samples if s in sec_annot])
    if args.candidates:
        cand = {}
        for path in args.candidates:
            sample = Path(path).name[: -len(".candidates.tsv")]
            with open(path) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            cand[sample] = {"candidates": len(rows), "with_tga_event": sum(1 for r in rows if r["events"] != "secis_only"),
                            "best": f"{rows[0]['gene']} ({rows[0]['rank_score']} bits)" if rows else "none"}
        write_table("readthrough_mqc.tsv", "readthrough", "Transcript screen",
                    "Annotated transcripts with a TGA that may encode Sec at the CDS edges or inside, ranked by the SECIS score of their 3' UTR.",
                    ["candidates", "with_tga_event", "best"], [(s, cand[s]) for s in samples if s in cand])

    expect = read_rows(args.expectation)
    if expect:
        write_table("expectation_mqc.tsv", "expectation", "Expected selenoproteome",
                    "Sec sites per family compared with the expectation for the lineage; missing families may be known blind spots of the predictors.",
                    ["expected", "observed", "missing", "extra"], [(s, expect[s]) for s in samples if s in expect])

    if args.transcripts:
        tx = {}
        for path in args.transcripts:
            sample = Path(path).name[: -len(".sec_transcripts.tsv")]
            with open(path) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            supported = [r for r in rows if int(r["transcripts"]) > 0]
            tx[sample] = {"sec_sites": len(rows), "with_transcripts": len(supported),
                          "transcripts_total": sum(int(r["transcripts"]) for r in rows),
                          "median_3prime_extent": statistics.median(int(r["max_3prime_extent"]) for r in supported) if supported else "NA"}
        write_table("transcript_support_mqc.tsv", "transcript_support", "Transcript support",
                    "Sec codons covered by the given transcripts: number of Sec codons, codons inside an aligned exon of at least one transcript, "
                    "supporting alignments in total, and the median extension of the longest supporting transcript 3' of the codon (nt).",
                    ["sec_sites", "with_transcripts", "transcripts_total", "median_3prime_extent"], [(s, tx[s]) for s in samples if s in tx])

    with open("sec_trait_mqc.tsv", "w") as out:
        out.write("# id: 'sec_trait'\n# plot_type: 'generalstats'\n# pconfig:\n#     sec_trait:\n#         title: 'Sec trait'\n#         description: 'Sec trait call from selenoprotein genes and tRNA-Sec'\n")
        out.write("Sample\tsec_trait\n")
        for sample in samples:
            out.write(f"{sample}\t{trait(genes.get(sample), trnasec.get(sample))}\n")


if __name__ == "__main__":
    main()
