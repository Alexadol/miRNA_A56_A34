#!/usr/bin/env python3
"""
assign_names.py

Assigns systematic names to the annotated MIRNA loci of pea, following the
nomenclature used by miRBase, and writes the naming table.

Rules applied
-------------
Prefix. Every name begins with psa-, the three-letter code for Pisum sativum.
Precursor loci are written MIR in capitals (psa-MIR156a) and the mature
products in mixed case (psa-miR156a-5p); this file gives both.

Letters. Within a family, loci are ordered by chromosome and then by start
coordinate, and are lettered a, b, c in that order. Ordering by position rather
than by abundance keeps the names stable when read counts change.

Identical products. Where two or more precursors yield the same functional
sequence, they are not given separate letters. They share one letter and are
distinguished by a numeric suffix, psa-MIR156a-1 and psa-MIR156a-2, which is
the suffix miRBase reserves for exactly this case. Distinct sequences within a
family take distinct letters.

Single-locus families. A family represented by one locus takes no letter,
following miRBase practice: psa-MIR162, not psa-MIR162a.

Arm suffix. The mature name carries -5p or -3p according to the position of the
functional arm in the hairpin. Loci with no functional arm assigned receive the
precursor name only.

Unassigned loci. Loci with no homology to any database are numbered
psa-MIRN1 upward, again by chromosome and position. These names are provisional
and intended only for internal reference and for the manuscript; miRBase
assigns the definitive numbers on deposition, and the provisional names should
not be cited as if they were registered. A column records whether each of these
loci has nonetheless been reported in a published pea dataset, since those are
not novel even though they are absent from the databases.

Input
-----
  --functional   annotation_functional.tsv from assign_functional_arm.py
  --coordinates  annotation_curated.tsv, for chrom, start, end and strand
  --reported     optional, a file listing one cluster identifier per line that
                 was recovered from a published pea dataset

Output
------
  A tab-separated table with the columns

    name, precursor_name, cluster, family, chrom, start, end, strand,
    functional_arm, functional_seq, evidence_level, status

  followed by a commented legend.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

LETTERS = "abcdefghijklmnopqrstuvwxyz"

LEGEND = """
# ---------------------------------------------------------------------------
# Column definitions
#
# name              Name of the mature product, with the arm suffix -5p or -3p.
#                   Absent where no functional arm could be assigned.
# precursor_name    Name of the precursor locus, written MIR in capitals.
# cluster           Locus identifier assigned by ShortStack, retained so that
#                   the names can be traced back to the primary output.
# family            microRNA family. Empty for loci with no database homology.
# chrom start end   Position of the precursor, 1-based inclusive, on the strand
# strand            given in the strand column.
# functional_arm    Position in the hairpin of the arm taken as the functional
# functional_seq    product, and its sequence. See annotation_functional.tsv
#                   for how the assignment was made.
# evidence_level    PARE, homology or none, as in annotation_functional.tsv.
# status            named          assigned to a described family;
#                   unassigned     no database homology, provisional psa-MIRN
#                                  name, reported in a published pea dataset;
#                   unassigned_new no database homology and not found in any
#                                  published pea dataset.
#
# Letters within a family run in order of chromosome and start coordinate.
# Precursors yielding an identical functional sequence share a letter and are
# distinguished by a numeric suffix, which is the miRBase convention for that
# case. Families represented by a single locus take no letter.
#
# psa-MIRN names are provisional. Definitive numbers are assigned by miRBase on
# deposition and should replace these before publication.
# ---------------------------------------------------------------------------
"""


def chrom_key(c):
    m = re.search(r"(\d+)", c or "")
    return (0, int(m.group(1))) if m and c.startswith("chr") else (1, 0, c)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--functional", required=True)
    ap.add_argument("--coordinates", required=True)
    ap.add_argument("--reported")
    ap.add_argument("--prefix", default="psa")
    ap.add_argument("--out", default="mirna_names.tsv")
    args = ap.parse_args()

    with open(args.functional, encoding="utf-8-sig") as fh:
        fn = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    with open(args.coordinates, encoding="utf-8-sig") as fh:
        co = {r["cluster"]: r for r in
              csv.DictReader((l.replace("\r", "") for l in fh), delimiter="\t")}
    reported = set()
    if args.reported:
        reported = {l.strip() for l in open(args.reported) if l.strip()}
    print(f"loci: {len(fn)}; reported in published datasets: {len(reported)}", file=sys.stderr)

    for r in fn:
        c = co.get(r["cluster"], {})
        r["chrom"], r["start"] = c.get("chrom", ""), int(c.get("start", 0) or 0)
        r["end"], r["strand"] = c.get("end", ""), c.get("strand", "")

    def pos(r):
        return (chrom_key(r["chrom"]), r["start"])

    named, p = [], args.prefix

    # --- loci assigned to a family
    byfam = defaultdict(list)
    for r in fn:
        if r["family"]:
            byfam[r["family"]].append(r)

    for fam, loci in byfam.items():
        loci.sort(key=pos)
        core = fam[3:] if fam.startswith("miR") else fam
        # group by functional sequence, keeping first appearance order
        groups, seen = [], {}
        for r in loci:
            k = r["functional_seq"] or f"__{r['cluster']}"
            if k not in seen:
                seen[k] = len(groups)
                groups.append([])
            groups[seen[k]].append(r)
        single = len(groups) == 1
        for gi, grp in enumerate(groups):
            letter = "" if single else LETTERS[gi] if gi < 26 else f"z{gi}"
            for j, r in enumerate(grp):
                suffix = "" if len(grp) == 1 else f"-{j+1}"
                pre = f"{p}-MIR{core}{letter}{suffix}"
                mat = (f"{p}-miR{core}{letter}{suffix}-{r['functional_arm']}"
                       if r["functional_arm"] else "")
                named.append((mat, pre, r, "named"))

    # --- loci with no database homology
    rest = sorted((r for r in fn if not r["family"]), key=pos)
    for i, r in enumerate(rest, 1):
        pre = f"{p}-MIRN{i}"
        mat = f"{p}-miRN{i}-{r['functional_arm']}" if r["functional_arm"] else ""
        st = "unassigned" if r["cluster"] in reported else "unassigned_new"
        named.append((mat, pre, r, st))

    named.sort(key=lambda x: (chrom_key(x[2]["chrom"]), x[2]["start"]))

    hdr = ["name", "precursor_name", "cluster", "family", "chrom", "start", "end",
           "strand", "functional_arm", "functional_seq", "evidence_level", "status"]
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        for mat, pre, r, st in named:
            w.writerow([mat, pre, r["cluster"], r["family"], r["chrom"], r["start"],
                        r["end"], r["strand"], r["functional_arm"],
                        r["functional_seq"], r["evidence_level"], st])
        fh.write(LEGEND)

    from collections import Counter
    print(f"\nnamed loci: {sum(1 for x in named if x[3] == 'named')}", file=sys.stderr)
    print(f"provisional psa-MIRN names: {sum(1 for x in named if x[3].startswith('unassigned'))}",
          file=sys.stderr)
    print(f"without a mature name (no functional arm): "
          f"{sum(1 for x in named if not x[0])}", file=sys.stderr)
    dups = [k for k, v in Counter(x[1] for x in named).items() if v > 1]
    print(f"duplicate precursor names: {len(dups)} {dups}", file=sys.stderr)
    print(f"\nwritten: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
