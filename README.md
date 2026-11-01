# Strand-resolved microRNA annotation of the pea (Pisum sativum) root

Scripts accompanying the manuscript on the strand-resolved microRNA annotation
of pea roots and its early response to rhizobial inoculation.

The central result is methodological: the arm selected by read abundance is not
a reliable guide to the functional microRNA strand. The pipeline therefore
assigns the functional arm from degradome evidence where available, and from
database homology otherwise, and quantifies the two arms separately.

## Scripts

Run in this order. Each script prints its own usage with `-h`.

| Script | Purpose |
| --- | --- |
| `build_pmiren_db.py` | Build reference FASTA from the PmiREN 3.1 archive, removing species-specific entries that would corrupt family assignment. |
| `annotate_mirna_arms.py` | Annotate both arms of every locus against miRBase and PmiREN; call the conserved arm by match strength, not abundance. |
| `curate_families.py` | Manual curation: restore conserved families overridden by species-specific entries, drop animal hits, set database priority. |
| `assign_functional_arm.py` | Assign the functional arm (degradome-validated, then homology, then none), using the stringent set of Zhai et al. 2011. |
| `assign_names.py` | Assign systematic psa-MIR names by chromosomal position, following miRBase conventions. |
| `fig_de_volcano.R` | Volcano overview of the differential-abundance analysis. |
| `mirna_de_pipeline.R` | Differential-abundance analysis of arm groups with DESeq2 across the three contrasts. |

## Inputs

- `mir.fasta` ShortStack FASTA with hairpin, mature and star records per locus
- `Results.txt` ShortStack output, for per-locus read counts
- `mature.fa` miRBase mature sequences (https://www.mirbase.org/download/)
- PmiREN 3.1 archive (https://pmiren.com)
- `arm_group_counts.tsv` per-arm read counts for the differential analysis
- Supplementary Table S4 of Zhai et al. 2011 (degradome-validated targets)

## Dependencies

- Python 3.9+: pandas, openpyxl, xlrd>=2.0.1 (for the .xls degradome table)
- R 4.4: DESeq2, tidyverse
- External tools (versions in Materials and Methods): ShortStack 4.0.2,
  Infernal 1.1.5, bowtie, samtools, bedtools

## Example

```bash
python3 build_pmiren_db.py FTP_download -o mirna_db --legumes-only
python3 annotate_mirna_arms.py --loci mir.fasta --db mature.fa \
    mirna_db/pmiren_legumes_mature.fa --results Results.txt --out annotation_final.tsv
python3 curate_families.py annotation_final.tsv annotation_curated.tsv
python3 assign_functional_arm.py --annotation annotation_curated.tsv \
    --pare Table_S4_PARE_validation.xls --out annotation_functional.tsv
python3 assign_names.py --functional annotation_functional.tsv \
    --coordinates annotation_curated.tsv --out mirna_names.tsv
```


## Figure scripts

| Script | Figure |
| --- | --- |
| `fig1_annotation_overview.R` | Figure 1b (UpSet); panel 1c is a template, panel 1a is a schematic |
| `fig2_size_classes.R` | Figure 2a, 2b (size-class bars); panel 2c printed to console |
| `fig3_arm_asymmetry.R` | Figure 3 (arm-usage asymmetry) |
| `figS_de_volcano.R` | Supplementary volcano of the differential analysis |

Figure scripts additionally need: `ComplexUpset`, `ggrepel`, `patchwork`.
Inputs for the figures: `size_by_arm_long.tsv`, `arm_group_counts.tsv`,
`arm_group_members.tsv`, `byloc.tsv`, `DE_arm_groups_summary.tsv`,
`annotation_functional.tsv`, `mirna_names.tsv`.

## License

MIT
