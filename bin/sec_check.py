#!/usr/bin/env python3
"""Report genomic TGA codons that miniprot aligns to selenocysteine (U) residues of the queries.

A Sec site is supported when the query U is aligned to a TGA codon ("aligned"), or when U is
the first query residue after the end of the alignment and the next genomic codon is TGA
("next_codon"): miniprot usually stops before a C-terminal Sec because of its in-frame stop
penalty. Codons split by an intron are not detected.

Output columns (1-based, inclusive coordinates of the TGA codon):
contig, start, end, strand, family, query, support, score, identity
"""

import argparse
import gzip
import re
import sys

CS_TOKEN = re.compile(r"(:\d+)|(\*[acgtn]+[A-Z*])|(\+[A-Z*]+)|(-[acgtn]+)|(~[acgtn]{2}\d+[acgtn]{2})")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path, wanted=None):
    sequences, name, chunks = {}, None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif wanted is None or name in wanted:
                chunks.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        sequences[name] = "".join(chunks)
    return sequences


def parse_paf(gff_path):
    """Yield the alignment records embedded as ##PAF lines in miniprot GFF3 output."""
    with open_text(gff_path) as handle:
        for line in handle:
            if not line.startswith("##PAF\t"):
                continue
            fields = line.rstrip("\n").split("\t")[1:]
            tags = dict(field.split(":", 2)[0::2] for field in fields[12:] if field.count(":") >= 2)
            yield {
                "query": fields[0],
                "qlen": int(fields[1]),
                "qstart": int(fields[2]),
                "qend": int(fields[3]),
                "strand": fields[4],
                "contig": fields[5],
                "tstart": int(fields[7]),
                "tend": int(fields[8]),
                "identity": int(fields[9]) / max(int(fields[10]), 1),
                "score": int(tags.get("AS", 0)),
                "cs": tags.get("cs", ""),
            }


def aligned_sec_codons(record, query_seq):
    """Return 0-based half-open genomic intervals of TGA codons aligned to query U."""
    sites = []
    qpos = record["qstart"]
    forward = record["strand"] == "+"
    gpos = record["tstart"] if forward else record["tend"]

    def consume(n):
        nonlocal gpos
        interval = (gpos, gpos + n) if forward else (gpos - n, gpos)
        gpos = gpos + n if forward else gpos - n
        return interval

    for match in CS_TOKEN.finditer(record["cs"]):
        token = match.group(0)
        if token[0] == ":":
            n = int(token[1:])
            consume(3 * n)
            qpos += n
        elif token[0] == "*":
            codon = token[1:-1]
            interval = consume(len(codon))
            if codon == "tga" and qpos < len(query_seq) and query_seq[qpos] == "U":
                sites.append(interval)
            qpos += 1
        elif token[0] == "+":
            qpos += len(token) - 1
        elif token[0] == "-":
            consume(len(token) - 1)
        else:
            consume(int(token[3:-2]))
    return sites, qpos


def next_codon(record, genome):
    contig = genome.get(record["contig"])
    if contig is None:
        return None, None
    if record["strand"] == "+":
        start = record["tend"]
        codon = contig[start : start + 3].upper()
    else:
        start = record["tstart"] - 3
        codon = contig[max(start, 0) : record["tstart"]].upper()[::-1].translate(COMPLEMENT)
    if start < 0 or len(codon) != 3:
        return None, None
    return codon, (start, start + 3)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True, help="miniprot GFF3 output with ##PAF lines")
    parser.add_argument("--queries", required=True, help="Protein queries given to miniprot")
    parser.add_argument("--genome", required=True, help="Genome FASTA, optionally gzip-compressed")
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    queries = read_fasta(args.queries)
    records = list(parse_paf(args.gff))
    needs_sequence = {
        r["contig"] for r in records if r["qend"] < r["qlen"] and queries.get(r["query"], "")[r["qend"] : r["qend"] + 1] == "U"
    }
    genome = read_fasta(args.genome, wanted=needs_sequence) if needs_sequence else {}

    best = {}
    for record in records:
        query_seq = queries.get(record["query"])
        if query_seq is None:
            sys.exit(f"Query {record['query']} not found in {args.queries}")
        sites, qpos = aligned_sec_codons(record, query_seq)
        found = [(interval, "aligned") for interval in sites]
        if qpos == record["qend"] and qpos < len(query_seq) and query_seq[qpos] == "U":
            codon, interval = next_codon(record, genome)
            if codon == "TGA":
                found.append((interval, "next_codon"))
        family = record["query"].split("|", 1)[0]
        for (start, end), support in found:
            key = (record["contig"], start, record["strand"], family)
            row = (record["contig"], start + 1, end, record["strand"], family, record["query"], support, record["score"], f"{record['identity']:.4f}")
            if key not in best or record["score"] > best[key][7]:
                best[key] = row

    with open(args.output, "w") as out:
        out.write("contig\tstart\tend\tstrand\tfamily\tquery\tsupport\tscore\tidentity\n")
        for row in sorted(best.values(), key=lambda r: (r[0], r[1], r[3], r[4])):
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
