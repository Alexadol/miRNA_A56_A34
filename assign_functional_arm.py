#!/usr/bin/env python3
"""
assign_functional_arm.py

Assigns the functional arm of each MIRNA locus under a tiered evidence scheme
and writes a compact annotation table.

Why a tiered scheme
-------------------
None of the conventional criteria identifies which strand of the duplex is
loaded into ARGONAUTE and acts on targets:

  position (5p/3p)   carries no information about function; the functional
                     product of a family may derive from either arm;
  abundance          reflects steady-state accumulation of each product rather
                     than loading, and both strands of a duplex can be
                     retained;
  database labels    are derived from the same abundance criterion and, being
                     compiled from many studies, are not applied uniformly.

The functional arm is therefore assigned in descending order of confidence:

  PARE       the arm matches a small RNA whose cleavage activity was validated
             by parallel analysis of RNA ends in Medicago truncatula
             (Zhai et al. 2011, Genes Dev 25:2540, Table S4A, stringent set);
  homology   the arm is not represented in the degradome data, and the arm
             carrying the family-conserved sequence is assigned by homology to
             miRBase and PmiREN;
  none       no homology is detected and no arm is assigned.

Relative abundance is recorded but is not used as a criterion.

Matching rule
-------------
A sequence matches a validated one when the two are identical, or when one is
contained within the other, which accommodates the differing arm boundaries of
the pipelines involved. Exact matches take precedence: where both arms of a
locus match validated sequences, an exact match on one arm overrides a
containment match on the other, and the locus is left unassigned if neither
match is exact.

Input
-----
  --annotation  annotation_curated.tsv, the output of curate_families.py.
                Must contain the columns cluster, family, mature_arm,
                mature_seq, star_arm, star_seq, canonical_strand, and the
                homology fields mature_hit, mature_db, mature_mm and their
                star counterparts.
  --pare        Table S4 of Zhai et al. 2011, read from the sheet named
                A.Stringent with the header on the fourth row. Requires pandas,
                and xlrd 2.0.1 or later for the .xls format.

Output
------
  A tab-separated table with one row per locus and the columns

    cluster, family,
    mature_arm, mature_seq, star_arm, star_seq,
    functional_arm, functional_strand, functional_seq,
    evidence_level, evidence_source

  followed by a commented legend defining every column, including how
  ShortStack derives the mature and star designations. Lines of the legend
  begin with a hash and are ignored by readers configured to skip comments.

Example
-------
  python3 assign_functional_arm.py \
      --annotation annotation_curated.tsv \
      --pare Table_S4_PARE_validation.xls \
      --out annotation_functional.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

HEADER = ["cluster", "family",
          "mature_arm", "mature_seq", "star_arm", "star_seq",
          "functional_arm", "functional_strand", "functional_seq",
          "evidence_level", "evidence_source"]

LEGEND = """
# ---------------------------------------------------------------------------
# Column definitions
#
# cluster            Locus identifier assigned by ShortStack.
# family             microRNA family, from homology to miRBase and PmiREN,
#                    curated as described in Methods.
#
# mature_arm         The arm of the precursor that ShortStack designates as the
# star_arm           mature microRNA, and the complementary arm, designated the
#                    star strand. ShortStack makes this call by abundance: the
#                    single most frequent distinct sequence at the locus
#                    (MajorRNA) is taken as the mature product, and the arm on
#                    which it lies becomes the mature arm. The designation is
#                    therefore a statement about steady-state accumulation, not
#                    about which strand is loaded into ARGONAUTE.
#
# 5p / 3p            Position of the arm within the hairpin, 5-prime or
#                    3-prime of the loop, read on the annotated strand. This is
#                    a coordinate, not a functional attribute: the functional
#                    product of a family may derive from either arm.
#
# functional_arm     The arm taken as the functional product of the locus, with
# functional_strand  its position, its ShortStack designation, and its sequence.
# functional_seq     Empty where no assignment could be made.
#
# evidence_level     How the functional arm was assigned.
#   PARE             The arm matches a small RNA whose cleavage activity was
#                    validated by parallel analysis of RNA ends in Medicago
#                    truncatula (Zhai et al. 2011, Genes Dev 25:2540,
#                    Table S4A, stringent set). Highest confidence.
#   homology         Not represented in the degradome data; the arm carrying
#                    the family-conserved sequence was assigned by homology to
#                    miRBase and PmiREN.
#   none             No homology detected; no arm assigned.
#
# evidence_source    For PARE, the name of the validated sequence and whether
#                    the match was exact or by containment of one sequence in
#                    the other, which accommodates differing arm boundaries.
#                    For homology, the database entry matched, its source, and
#                    the number of mismatches.
#
# Relative abundance was recorded during annotation but was not used as a
# criterion for assigning the functional arm.
# ---------------------------------------------------------------------------
"""


def load_pare(path):
    """Sequence -> set of names from Table S4A."""
    try:
        import pandas as pd
    except ImportError:
        sys.exit("pandas is required; xlrd>=2.0.1 is also needed for .xls")
    df = pd.read_excel(path, sheet_name="A.Stringent", header=3)
    df = df.dropna(subset=["Sequence"])
    pare = {}
    for _, r in df.iterrows():
        pare.setdefault(str(r["Sequence"]).strip().upper(), set()).add(str(r["Annotation"]).strip())
    print(f"PARE-validated sequences: {len(pare)}", file=sys.stderr)
    return pare


def pare_hit(seq, pare):
    """Exact match or containment of one sequence within the other."""
    s = (seq or "").strip().upper().replace("U", "T")
    if len(s) < 18:
        return None
    best = None
    for k, names in pare.items():
        if s == k:
            return ("exact", k, names)
        if (s in k or k in s) and best is None:
            best = ("contained", k, names)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--pare", required=True)
    ap.add_argument("--out", default="annotation_functional.tsv")
    args = ap.parse_args()

    pare = load_pare(args.pare)

    with open(args.annotation, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader((l.replace("\r", "") for l in fh), delimiter="\t"))
    print(f"input loci: {len(rows)}", file=sys.stderr)

    out, stats, flips = [], {"PARE": 0, "homology": 0, "none": 0}, []
    both = 0
    for r in rows:
        hm = pare_hit(r["mature_seq"], pare)
        hs = pare_hit(r["star_seq"], pare)

        # both arms in the degradome: take the exact match, otherwise leave unresolved
        if hm and hs:
            both += 1
            if hm[0] == "exact" and hs[0] != "exact":
                hs = None
            elif hs[0] == "exact" and hm[0] != "exact":
                hm = None

        if hm and not hs:
            strand, arm, seq, lvl = "mature", r["mature_arm"], r["mature_seq"], "PARE"
            pm = f"{sorted(hm[2])[0]} ({hm[0]})"
        elif hs and not hm:
            strand, arm, seq, lvl = "star", r["star_arm"], r["star_seq"], "PARE"
            pm = f"{sorted(hs[2])[0]} ({hs[0]})"
        elif r["mature_hit"] or r["star_hit"]:
            cs = r["canonical_strand"]
            if cs == "mature":
                strand, arm, seq = "mature", r["mature_arm"], r["mature_seq"]
                hit, db, mm = r.get("mature_hit", ""), r.get("mature_db", ""), r.get("mature_mm", "")
            elif cs == "star":
                strand, arm, seq = "star", r["star_arm"], r["star_seq"]
                hit, db, mm = r.get("star_hit", ""), r.get("star_db", ""), r.get("star_mm", "")
            else:
                strand, arm, seq, hit, db, mm = "", "", "", "", "", ""
            lvl = "homology"
            pm = f"{hit} [{db}] {mm}mm" if hit else ""
        else:
            strand, arm, seq, lvl, pm = "", "", "", "none", ""

        stats[lvl] += 1
        if lvl == "PARE" and strand and r["canonical_strand"] != strand:
            flips.append((r["cluster"], r["family"], r["canonical_strand"], strand))

        out.append([r["cluster"], r["family"],
                    r["mature_arm"], r["mature_seq"], r["star_arm"], r["star_seq"],
                    arm, strand, seq, lvl, pm])

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEADER)
        w.writerows(out)
        fh.write(LEGEND)

    print(f"\nevidence levels: {stats}", file=sys.stderr)
    print(f"loci where both arms matched the degradome: {both}", file=sys.stderr)
    print(f"\nreassigned relative to abundance- and homology-based calls: {len(flips)}", file=sys.stderr)
    for c, f, old, new in flips:
        print(f"  {c:<15} {f or '-':<9} {old or '-':<7} -> {new}", file=sys.stderr)
    nf = sum(1 for o in out if o[7] == "star")
    nm = sum(1 for o in out if o[7] == "mature")
    print(f"\nfunctional arm: mature {nm}, star {nf}, unassigned {len(out)-nm-nf}", file=sys.stderr)
    if nm + nf:
        print(f"  functional strand is the minor arm at {nf} of {nm+nf} loci "
              f"({100*nf/(nm+nf):.0f} %)", file=sys.stderr)
    print(f"\nwritten: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
