#!/usr/bin/env python3
"""
curate_families.py

Manual curation of family assignments on top of the automatic annotation.

Why
---
The automatic choice "best hit = fewest mismatches" over the PmiREN database
systematically loses conserved families to species-specific entries: an exact
hit to Cre-miR13821 from chickpea overrides an equally exact hit to
ath-miR168a-3p from miRBase. This edit restores the conserved family where it
is supported by miRBase, and marks as unannotated the loci with no conserved
alternative.

Each decision below was made from the data: the best hit in the combined
database, the alternative hits, and the result of a separate miRBase-only run
were compared.

Input:  annotation_final.tsv (output of annotate_mirna_arms.py)
Output: annotation_curated.tsv plus an edit log
"""

import csv
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Loci where the conserved family is supported by miRBase and must win.
# Format: cluster -> (family, miRBase hit, mm in miRBase, reason)
# ---------------------------------------------------------------------------
RESTORE = {
    "Cluster_5972":  ("miR168",  "ath-miR168a-3p", 0, "exact miRBase hit, overridden by Cre-miR13821 (0mm)"),
    "Cluster_7606":  ("miR393",  "osa-miR393b-5p", 0, "exact miRBase hit, overridden by Amo-miR15506a (0mm)"),
    "Cluster_6190":  ("miR2119", "mtr-miR2119",    0, "exact miRBase hit, overridden by Mru-miR15401a (0mm)"),
    "Cluster_15810": ("miR396",  "ath-miR396a-3p", 0, "exact miRBase hit, overridden by Gma-miR13486a (0mm)"),
    "Cluster_21785": ("miR319",  "gma-miR319d",    0, "exact miRBase hit, overridden by Ahy-miR15458a (0mm)"),
    "Cluster_21979": ("miR166",  "osa-miR166h-5p", 0, "exact miRBase hit, overridden by Ahy-miR15217a (0mm)"),
    "Cluster_29571": ("miR171",  "osa-miR171i-3p", 0, "exact miRBase hit, overridden by Car-miR15289b (0mm)"),
    "Cluster_9874":  ("miR164",  "csi-miR164a-3p", 1, "conserved family at 1 mismatch vs Amo-miR13269a (0mm)"),
    "Cluster_13079": ("miR164",  "stu-miR164-3p",  1, "conserved family at 1 mismatch vs Vsa-miR13161 (0mm)"),
    "Cluster_19878": ("miR156",  "csi-miR156g-3p", 2, "conserved family, both hits at 2 mismatches"),
    "Cluster_21708": ("miR171",  "ptc-miR171l-5p", 2, "conserved family vs Mru-miR15293 (0mm); Car-miR170 (2mm) same group"),
    "Cluster_25422": ("miR5232", "mtr-miR5232",    2, "conserved family vs Car-miR14005 (1mm)"),
}

# ---------------------------------------------------------------------------
# Loci with no conserved alternative: the family is removed and the locus
# becomes unannotated. A hit to a species-specific entry of another
# species is not a family assignment.
# ---------------------------------------------------------------------------
UNASSIGN = {
    "Cluster_10284": "only hit Mru-miR15633a, species-specific; no miRBase hit",
    "Cluster_13318": "hit on the star arm only, no mature hit",
    "Cluster_18767": "only hit Vun-miR12789 (2mm), species-specific; no miRBase hit",
    "Cluster_946":   "5S rRNA fragment by Rfam, to be excluded from the set",
}

# ---------------------------------------------------------------------------
# Loci where a high-numbered family is also supported by miRBase, i.e. the
# assignment is correct and needs no edit. Listed explicitly to separate
# checked cases from unchecked ones.
# ---------------------------------------------------------------------------
CONFIRMED = {
    "Cluster_7374":  ("miR4414",  "mtr-miR4414a-5p, 3mm"),
    "Cluster_8785":  ("miR4415",  "gma-miR4415a-3p, 1mm"),
    "Cluster_9675":  ("miR5770",  "gma-miR5770b, 0mm"),
    "Cluster_20052": ("miR10999", "mdm-miR10999a, 4mm, weak"),
    "Cluster_31533": ("miR6425",  "ptc-miR6425a-3p, 3mm"),
    "Cluster_31534": ("miR6425",  "ptc-miR6425a-3p, 3mm"),
    "Cluster_31535": ("miR6425",  "ptc-miR6425a-3p, 3mm"),
}

# Merge related families conventionally treated as one.
MERGE = {
    "miR170": "miR171",   # miR170 and miR171 target GRAS/SCL
    "miR165": "miR166",   # miR165 and miR166 target HD-ZIP III
}

# Loci excluded from the final set entirely.
EXCLUDE = {"Cluster_946": "5S rRNA by Rfam cmscan"}

# ---------------------------------------------------------------------------
# Second edit block: reassigning the conserved arm.
#
# Two systematic problems of the automatic run:
#
# 1. miRBase mature.fa contains ALL organisms, not only plants.
#    Twelve loci got their best hit to an animal miRNA (mouse, chicken,
#    frog, nematode), all at 4 mismatches. Animal names use a hyphen
#    (mmu-miR-669b-5p), plant names do not (ath-miR168a-3p), which is
#    used for filtering.
#
# 2. No database priority, so an exact hit to a curated miRBase entry is
#    indistinguishable from an exact hit to a species-specific PmiREN
#    entry, and the arms are called equal. This inflated the equal class
#    to 41 % of reads.
# ---------------------------------------------------------------------------

ANIMAL_RE = re.compile(r"miR-\d")          # animal miRBase nomenclature


def is_animal(hit):
    return bool(hit) and bool(ANIMAL_RE.search(hit))


def is_mirbase(hit):
    """PmiREN entries carry a species prefix with | ; miRBase entries do not."""
    return bool(hit) and "|" not in hit


def arm_rank(hit, mm, alen):
    """
    Arm sort key: smaller is better.

    Order: mismatches, then source, then overlap length. Source is used ONLY
    as a tiebreaker at equal mismatch count: an exact hit to a curated miRBase
    entry outweighs an exact hit to a species-specific PmiREN entry. Source
    must not precede mismatches, or a 3-mismatch miRBase hit would beat a
    0-mismatch PmiREN hit and the assignment would flip without cause.
    """
    if not hit or is_animal(hit):
        return (99, 9, 0)
    source = 0 if is_mirbase(hit) else 1
    try:
        mm = int(mm)
    except (TypeError, ValueError):
        mm = 99
    try:
        alen = int(alen)
    except (TypeError, ValueError):
        alen = 0
    return (mm, source, -alen)


def main():
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "annotation_final.tsv")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "annotation_curated.tsv")

    with open(src, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader((l.replace("\r", "") for l in fh), delimiter="\t"))
    fields = list(rows[0].keys()) + ["family_source", "curation_note"]

    log, kept = [], []
    for r in rows:
        cl = r["cluster"]
        before = r["family"]

        if cl in EXCLUDE:
            log.append(f"EXCLUDED   {cl:<15} {before or '(none)':<10} :: {EXCLUDE[cl]}")
            continue

        if cl in RESTORE:
            fam, hit, mm, why = RESTORE[cl]
            r["family"] = fam
            r["mature_hit"] = hit
            r["mature_db"] = "miRBase"
            r["mature_mm"] = str(mm)
            r["mature_confidence"] = {0: "exact", 1: "high", 2: "high", 3: "medium"}.get(mm, "low")
            r["family_source"] = "curated:miRBase"
            r["curation_note"] = why
            log.append(f"RESTORED  {cl:<15} {before or '(none)':<10} -> {fam:<10} :: {why}")
        elif cl in UNASSIGN:
            r["family"] = ""
            r["family_source"] = "curated:unassigned"
            r["curation_note"] = UNASSIGN[cl]
            log.append(f"REMOVED    {cl:<15} {before or '(none)':<10} -> (none)     :: {UNASSIGN[cl]}")
        elif cl in CONFIRMED:
            fam, why = CONFIRMED[cl]
            r["family_source"] = "confirmed:miRBase"
            r["curation_note"] = f"confirmed by miRBase ({why})"
            if r["family"] != fam:
                log.append(f"WARNING    {cl:<15} expected {fam}, file has {r['family']}")
        else:
            r["family_source"] = "auto"
            r["curation_note"] = ""

        # --- reassign the conserved arm ---
        for side in ("mature", "star"):
            if is_animal(r[f"{side}_hit"]):
                log.append(f"ANIMAL     {cl:<15} {side}: {r[f'{side}_hit']} "
                           f"({r[f"{side}_mm"]}mm) dropped")
                for f_ in ("hit", "db", "mm", "align_len", "offset"):
                    r[f"{side}_{f_}"] = ""
                r[f"{side}_confidence"] = "no_match"

        km = arm_rank(r["mature_hit"], r["mature_mm"], r["mature_align_len"])
        ks = arm_rank(r["star_hit"], r["star_mm"], r["star_align_len"])
        was = r["canonical_strand"]
        if km == (99, 9, 0) and ks == (99, 9, 0):
            new_strand, new_arm = "none", "?"
        elif km < ks:
            new_strand, new_arm = "mature", r["mature_arm"]
        elif ks < km:
            new_strand, new_arm = "star", r["star_arm"]
        else:
            new_strand, new_arm = "equal", "?"
        if new_strand != was:
            why = ("miRBase priority at equal mismatch count"
                   if km[0] == ks[0] else "after dropping animal hits")
            log.append(f"ARM        {cl:<15} {was:<7} -> {new_strand:<7} :: {why}")
            r["curation_note"] = (r["curation_note"] + "; " if r["curation_note"] else "") + \
                                 f"conserved arm reassigned ({why})"
        r["canonical_strand"], r["canonical_arm"] = new_strand, new_arm
        r["canonical_seq"] = {"mature": r["mature_seq"], "star": r["star_seq"]}.get(new_strand, "")

        if r["family"] in MERGE:
            merged = MERGE[r["family"]]
            log.append(f"MERGE      {cl:<15} {r['family']:<10} -> {merged}")
            r["curation_note"] = (r["curation_note"] + "; " if r["curation_note"] else "") + \
                                 f"family {r["family"]} merged into {merged}"
            r["family"] = merged

        kept.append(r)

    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(kept)

    print("\n".join(log))
    import collections
    fam = collections.Counter(r["family"] for r in kept)
    cs = collections.Counter(r["canonical_strand"] for r in kept)
    num = lambda f: int(re.sub(r"\D", "", f)) if f and re.sub(r"\D", "", f) else 0
    hi = sorted({r["family"] for r in kept if num(r["family"]) > 2500})
    tot = sum(int(r["reads"]) for r in kept if r["reads"])
    disc = sum(int(r["reads"]) for r in kept if r["canonical_strand"] == "star" and r["reads"])
    print(f"\n{'='*60}")
    print(f"loci in final set: {len(kept)}")
    print(f"families: {len([f for f in fam if f])}; without family: {fam.get('', 0)}")
    print(f"remaining species-specific (>2500): {len(hi)} {hi}")
    print(f"canonical_strand: {dict(cs)}")
    print(f"discordant: {cs['star']} of {len(kept)} ({100*cs['star']/len(kept):.0f} %)")
    print(f"reads total {tot}, in discordant {disc} ({100*disc/tot:.1f} %)")
    unresolved = [r for r in kept if r["canonical_strand"] == "equal"]
    if unresolved:
        print(f"\nunresolved ties ({len(unresolved)}), both hits from PmiREN, "
              f"resolve manually:")
        for r in sorted(unresolved, key=lambda x: -int(x["reads"] or 0)):
            print(f"  {r['cluster']:<15} {r['family'] or '-':<9} reads={r['reads']:<8} "
                  f"m={r['mature_hit'][:30]:<30} ({r['mature_mm']}mm)  "
                  f"s={r['star_hit'][:30]} ({r['star_mm']}mm)")
    print(f"\nwritten: {out}")


if __name__ == "__main__":
    main()
