#!/usr/bin/env python3
"""
build_pmiren_db.py

Build reference FASTA files from the unpacked PmiREN 3.1 archive, removing
entries that would corrupt family assignment.

Why the filters
---------------
PmiREN holds about 87,000 loci from 340 taxa, most of them species-specific.
Searching the whole database, an exact hit is almost always found in some
obscure entry and beats the true conserved family: miR168 becomes miR13054
from Chinese cabbage, miR393 becomes miR15506 from Ammopiptanthus.

Filters, both on by default:
  --drop-mirn        remove entries with miRN in the name (species-specific
                     candidates that carry no family information).
  --min-species N    keep only families found in at least N species. A true
                     conserved family is present in many taxa; an annotation
                     artefact in one. Default 3.

Pea (Pisum_sativum, Lathyrus_oleraceus) is always excluded from the reference
and written to separate files: comparing an annotation with another annotation
of the same species is uninformative.

Usage
-----
  tar -xzf FTP_download.tar.gz
  python3 build_pmiren_db.py FTP_download -o ~/mirna_db --legumes-only
"""

import argparse
import collections
import re
import sys
from pathlib import Path

PEA = {"Pisum_sativum", "Lathyrus_oleraceus"}

LEGUME_PREFIXES = (
    "Acacia_", "Ammopiptanthus_", "Arachis_", "Astragalus_", "Cajanus_",
    "Cicer_", "Glycine_", "Lathyrus_", "Lotus_", "Macrotyloma_", "Medicago_",
    "Phaseolus_", "Pisum_", "Pongamia_", "Psophocarpus_", "Vicia_", "Vigna_",
)


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line.upper().replace("U", "T"))
    if name is not None:
        yield name, "".join(buf)


def family_of(name):
    """Gma-miR156a -> miR156 ; Mru-miRN85 -> None"""
    m = re.search(r"miR(\d+)", name, re.I)
    return f"miR{m.group(1)}" if m else None


def base_species(dirname):
    """
    Collapse cultivars, subspecies and varieties to the species, so that the
    --min-species threshold counts true species.

    PmiREN keeps cultivars as separate directories: Glycine_max-Amsoy,
    Glycine_max-Dong_Nong_No.50, Lactuca_sativa_cv._B198, Oryza_sativa_subsp._indica.
    Without collapsing, a family annotated in three soybean cultivars reaches
    three "species" and passes the filter although it is not conserved.

        Glycine_max-Dong_Nong_No.50     -> Glycine_max
        Lactuca_sativa_cv._B198         -> Lactuca_sativa
        Oryza_sativa_subsp._indica      -> Oryza_sativa
        Brassica_rapa_var._parachinensis-> Brassica_rapa
    """
    s = dirname.split("-")[0]
    s = re.split(r"_(?:subsp|var|cv)\.?_", s)[0]
    parts = s.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else s


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", help="path to the unpacked FTP_download")
    ap.add_argument("-o", "--out", default="pmiren_db")
    ap.add_argument("--min-species", type=int, default=3,
                    help="minimum number of species a family must occur in (default 3)")
    ap.add_argument("--keep-mirn", action="store_true",
                    help="do not drop entries with miRN in the name")
    ap.add_argument("--legumes-only", action="store_true")
    ap.add_argument("--count-on-all", action="store_true",
                    help="count species per family over the WHOLE database but apply the "
                         "threshold to the selected subset. Safer for legume-"
                         "specific families with --legumes-only")
    ap.add_argument("--no-collapse-cultivars", action="store_true",
                    help="do not collapse cultivars and subspecies when counting")
    ap.add_argument("--min-len", type=int, default=18)
    ap.add_argument("--max-len", type=int, default=26)
    args = ap.parse_args()

    info = Path(args.src) / "miRNA_info"
    if not info.is_dir():
        sys.exit(f"directory not found: {info}")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    species = sorted(d.name for d in info.iterdir() if d.is_dir())
    print(f"species in the archive: {len(species)}")
    pea_found = [s for s in species if s in PEA]
    print(f"pea directories: {pea_found or 'none'}")

    def selected(s):
        if s in PEA:
            return False
        if args.legumes_only and not s.startswith(LEGUME_PREFIXES):
            return False
        return True

    ref_species = [s for s in species if selected(s)]
    print(f"species in the reference: {len(ref_species)}")

    # first pass: how many true species per family
    norm = (lambda d: d) if args.no_collapse_cultivars else base_species
    count_species = [s for s in species if s not in PEA] if args.count_on_all else ref_species
    print(f"threshold counted over {len(count_species)} directories "
          f"({len(set(map(norm, count_species)))} species after collapsing cultivars)")

    fam_species = collections.defaultdict(set)
    for s in count_species:
        f = info / s / f"{s}_mature.fa"
        if not f.exists():
            continue
        for name, _ in read_fasta(f):
            fam = family_of(name)
            if fam:
                fam_species[fam].add(norm(s))

    keep_fams = {f for f, sp in fam_species.items() if len(sp) >= args.min_species}
    print(f"families total: {len(fam_species)}; passing threshold "
          f"({args.min_species}+ species): {len(keep_fams)}")
    dropped = sorted(fam_species.keys() - keep_fams,
                     key=lambda f: -len(fam_species[f]))[:10]
    if dropped:
        print(f"  examples dropped: {', '.join(dropped)}")

    # second pass: writing
    def write(kind, species_list, path, apply_filters):
        n_in = n_out = 0
        with open(path, "w") as fh:
            for s in species_list:
                f = info / s / f"{s}_{kind}.fa"
                if not f.exists():
                    continue
                for name, seq in read_fasta(f):
                    n_in += 1
                    if not (args.min_len <= len(seq) <= args.max_len):
                        continue
                    if apply_filters:
                        if not args.keep_mirn and re.search(r"miRN", name):
                            continue
                        fam = family_of(name)
                        if fam is None or fam not in keep_fams:
                            continue
                    fh.write(f">{s}|{name.split()[0]}\n{seq}\n")
                    n_out += 1
        print(f"  {path.name}: {n_out} of {n_in}")

    tag = "legumes" if args.legumes_only else "all"
    print(f"\nreference ({tag}):")
    write("mature", ref_species, out / f"pmiren_{tag}_mature.fa", True)
    write("star", ref_species, out / f"pmiren_{tag}_star.fa", True)

    print("\npea, separately, NOT for the reference:")
    write("mature", pea_found, out / "pmiren_pea_mature.fa", False)
    write("star", pea_found, out / "pmiren_pea_star.fa", False)
    for s in pea_found:
        bi = info / s / f"{s}_basicinfo.tsv"
        if bi.exists():
            (out / bi.name).write_bytes(bi.read_bytes())
            print(f"  copied {bi.name}")

    with open(out / f"family_species_counts_{tag}.tsv", "w") as fh:
        fh.write("family\tn_species\tkept\n")
        for f, sp in sorted(fam_species.items(), key=lambda x: -len(x[1])):
            fh.write(f"{f}\t{len(sp)}\t{'yes' if f in keep_fams else 'no'}\n")
    print(f"\nfamily table: family_species_counts_{tag}.tsv")


if __name__ == "__main__":
    main()
