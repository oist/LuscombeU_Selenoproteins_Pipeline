#!/usr/bin/env python3
"""Write the Sec-containing sequences of Selenoprofiles profile alignments as protein queries.

Each profile is an alignment <family>.fa with a <family>.fa.config file. Sequences are
ungapped, upper-cased and named <family>|<sequence id>; only sequences containing U are kept.
"""

import argparse
import sys
from pathlib import Path


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def profiles_folder(path):
    path = Path(path)
    nested = path / "selenoprotein_profiles"
    return nested if nested.is_dir() else path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", help="Selenoprofiles data directory or profiles folder")
    parser.add_argument("output", help="Output protein FASTA")
    args = parser.parse_args()

    folder = profiles_folder(args.profiles)
    configs = sorted(folder.glob("*.fa.config"))
    if not configs:
        sys.exit(f"No profile found in {folder}")
    n_written = 0
    with open(args.output, "w") as out:
        for config in configs:
            family = config.name[: -len(".fa.config")]
            seen = set()
            for seq_id, aligned in read_fasta(folder / f"{family}.fa"):
                sequence = aligned.replace("-", "").replace(".", "").upper()
                if "U" not in sequence or sequence in seen:
                    continue
                seen.add(sequence)
                out.write(f">{family}|{seq_id}\n{sequence}\n")
                n_written += 1
    print(f"{n_written} queries from {len(configs)} profiles", file=sys.stderr)


if __name__ == "__main__":
    main()
