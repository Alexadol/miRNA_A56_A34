#!/usr/bin/env python3
"""
annotate_mirna_arms.py

Annotate both arms of each MIRNA locus against reference databases (miRBase,
PmiREN, sRNAanno) and determine the conserved arm by match strength rather
than by abundance.

Input
-----
  --loci        ShortStack FASTA with three records per locus (hairpin,
                mature, star), header coordinates in BED convention.
  --db          One or more reference FASTA files of mature miRNAs; the
                database name is taken from the file name (*_db columns).
  --results     (optional) ShortStack Results.txt, to add per-locus reads.

Output
------
  TSV, one row per locus. Columns are described in HEADER below.

Matching
--------
  Gapless sliding alignment within --max-offset; mismatches counted over an
  overlap of at least --min-overlap. Best hit has the fewest mismatches,
  ties broken by the longest overlap. confidence from the mismatch count:
  0 exact, 1-2 high, 3 medium, 4 low, no hit no_match.

Conserved-arm assignment
------------------------
  --canonical-rule quality (default): compare by (mismatches, then overlap
  length). --canonical-rule length: by overlap length only, reproducing the
  earlier pipeline. No match -> none; a tie -> equal.

Example
-------
  python3 annotate_mirna_arms.py --loci mir.fasta --db mature.fa \\
      pmiren_legumes_mature.fa --results Results.txt --out annotation.tsv
"""

import argparse
import csv
import re
import sys
from pathlib import Path

HEADER = [
    "cluster", "chrom", "start", "end", "strand", "precursor_len",
    "mature_arm", "mature_seq", "mature_len",
    "star_arm", "star_seq", "star_len",
    "reads",
    "mature_hit", "mature_db", "mature_mm", "mature_align_len", "mature_offset", "mature_confidence",
    "star_hit", "star_db", "star_mm", "star_align_len", "star_offset", "star_confidence",
    "canonical_arm", "canonical_strand", "canonical_seq", "family",
    "mature_alt_hits", "star_alt_hits",
]

CONFIDENCE = {0: "exact", 1: "high", 2: "high", 3: "medium", 4: "low"}


def read_fasta(path):
    """Return a list of (header, sequence); upper-cased, U -> T."""
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line.upper().replace("U", "T"))
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def parse_loci(path):
    """Parse the ShortStack FASTA into a dict of loci."""
    coord_re = re.compile(r"^(?P<name>[^:]+)::(?P<chrom>\S+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)")
    loci = {}
    for hdr, seq in read_fasta(path):
        m = coord_re.match(hdr)
        if not m:
            print(f"skipping header without coordinates: {hdr}", file=sys.stderr)
            continue
        name = m.group("name")
        rec = dict(chrom=m.group("chrom"), start=int(m.group("start")),
                   end=int(m.group("end")), strand=m.group("strand"), seq=seq)
        if name.endswith(".mature"):
            loci.setdefault(name[:-7], {})["mature"] = rec
        elif name.endswith(".star"):
            loci.setdefault(name[:-5], {})["star"] = rec
        else:
            loci.setdefault(name, {})["hairpin"] = rec
    return loci


def load_dbs(paths):
    """Return a list of (db_name, entry_name, sequence)."""
    db = []
    for p in paths:
        label = Path(p).stem
        n = 0
        for hdr, seq in read_fasta(p):
            db.append((label, hdr.split()[0], seq))
            n += 1
        print(f"  database {label}: {n} entries", file=sys.stderr)
    return db


def best_match(query, db, max_mm, max_offset, min_overlap):
    """
    Gapless sliding alignment.
    Returns (hit_name, db_label, mismatches, align_len, offset, alt_hits).
    """
    best = None
    alts = []
    for label, name, ref in db:
        if abs(len(ref) - len(query)) > max_offset + 3:
            continue
        for off in range(-max_offset, max_offset + 1):
            qs, rs = (0, off) if off >= 0 else (-off, 0)
            overlap = min(len(query) - qs, len(ref) - rs)
            if overlap < min_overlap:
                continue
            mm = 0
            for i in range(overlap):
                if query[qs + i] != ref[rs + i]:
                    mm += 1
                    if mm > max_mm:
                        break
            if mm > max_mm:
                continue
            cand = (mm, -overlap, name, label, off)
            if best is None or cand < best:
                best = cand
            alts.append((mm, name))
    if best is None:
        return None
    mm, neg_len, name, label, off = best
    alt = sorted({f"{n}({m}mm)" for m, n in alts if n != name})[:4]
    return name, label, mm, -neg_len, off, "; ".join(alt)


def family_of(hit_name):
    """ath-miR397a -> miR397 ; Gma-miR1507c -> miR1507 ; psa-miRN12 -> ''"""
    if not hit_name:
        return ""
    m = re.search(r"miR(\d+)", hit_name, re.I)
    return f"miR{m.group(1)}" if m else ""


def load_reads(path):
    """Per-locus read count from ShortStack Results.txt."""
    reads = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        namecol = "Name" if "Name" in rd.fieldnames else rd.fieldnames[1]
        readcol = "Reads" if "Reads" in rd.fieldnames else None
        if readcol is None:
            print("no Reads column in Results.txt, read counts skipped", file=sys.stderr)
            return reads
        for r in rd:
            try:
                reads[r[namecol].strip()] = int(r[readcol])
            except (ValueError, KeyError):
                pass
    return reads


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--db", required=True, nargs="+")
    ap.add_argument("--results")
    ap.add_argument("--out", default="mirna_final_annotation.tsv")
    ap.add_argument("--max-mismatches", type=int, default=4)
    ap.add_argument("--max-offset", type=int, default=3)
    ap.add_argument("--min-overlap", type=int, default=18)
    ap.add_argument("--canonical-rule", choices=["quality", "length"], default="quality")
    args = ap.parse_args()

    print("reading loci", file=sys.stderr)
    loci = parse_loci(args.loci)
    print(f"  loci: {len(loci)}", file=sys.stderr)
    incomplete = [k for k, v in loci.items() if not {"hairpin", "mature", "star"} <= set(v)]
    if incomplete:
        print(f"  WARNING: incomplete records {len(incomplete)}: {incomplete[:8]}", file=sys.stderr)

    print("reading databases", file=sys.stderr)
    db = load_dbs(args.db)
    print(f"  total reference sequences: {len(db)}", file=sys.stderr)

    reads = load_reads(args.results) if args.results else {}

    rows = []
    for i, (cl, parts) in enumerate(sorted(loci.items(), key=lambda x: int(x[0].split("_")[1])), 1):
        if i % 10 == 0:
            print(f"  processed {i}/{len(loci)}", file=sys.stderr)
        hp = parts.get("hairpin")
        mat, star = parts.get("mature"), parts.get("star")
        if not (hp and mat and star):
            continue

        # arm is determined by its position in the hairpin, accounting for strand
        if hp["strand"] == "+":
            m_arm = "5p" if mat["start"] < star["start"] else "3p"
        else:
            m_arm = "5p" if mat["start"] > star["start"] else "3p"
        s_arm = "3p" if m_arm == "5p" else "5p"

        hits = {}
        for key, rec in (("mature", mat), ("star", star)):
            hits[key] = best_match(rec["seq"], db, args.max_mismatches, args.max_offset, args.min_overlap)

        def unpack(h):
            if h is None:
                return "", "", "", "", "", "no_match", ""
            name, label, mm, alen, off, alt = h
            return name, label, mm, alen, off, CONFIDENCE.get(mm, "low"), alt

        mh = unpack(hits["mature"])
        sh = unpack(hits["star"])

        # conserved arm
        hm, hs = hits["mature"], hits["star"]
        if hm is None and hs is None:
            canon_strand, canon_arm, canon_seq = "none", "?", ""
        elif hs is None:
            canon_strand, canon_arm, canon_seq = "mature", m_arm, mat["seq"]
        elif hm is None:
            canon_strand, canon_arm, canon_seq = "star", s_arm, star["seq"]
        else:
            if args.canonical_rule == "quality":
                km, ks = (hm[2], -hm[3]), (hs[2], -hs[3])
            else:
                km, ks = (-hm[3],), (-hs[3],)
            if km < ks:
                canon_strand, canon_arm, canon_seq = "mature", m_arm, mat["seq"]
            elif ks < km:
                canon_strand, canon_arm, canon_seq = "star", s_arm, star["seq"]
            else:
                canon_strand, canon_arm, canon_seq = "equal", "?", ""

        fam = family_of(mh[0]) or family_of(sh[0])

        rows.append([
            cl, hp["chrom"], hp["start"] + 1, hp["end"], hp["strand"], len(hp["seq"]),
            m_arm, mat["seq"], len(mat["seq"]),
            s_arm, star["seq"], len(star["seq"]),
            reads.get(cl, ""),
            mh[0], mh[1], mh[2], mh[3], mh[4], mh[5],
            sh[0], sh[1], sh[2], sh[3], sh[4], sh[5],
            canon_arm, canon_strand, canon_seq, fam,
            mh[6], sh[6],
        ])

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEADER)
        w.writerows(rows)

    import collections
    cs = collections.Counter(r[26] for r in rows)
    cf = collections.Counter(r[18] for r in rows)
    disc = cs["star"]
    print(f"\nrows written: {len(rows)} -> {args.out}", file=sys.stderr)
    print(f"arm consensus: {dict(cs)}", file=sys.stderr)
    print(f"  discordant (conserved sequence on the minor arm): "
          f"{disc} of {len(rows)} ({100*disc/len(rows):.0f} %)", file=sys.stderr)
    print(f"confidence by mature: {dict(cf)}", file=sys.stderr)
    print(f"without family: {sum(1 for r in rows if not r[28])}", file=sys.stderr)


if __name__ == "__main__":
    main()
