# AI-Assisted Staphylococcus Invasome Mining Pipeline

Reusable computational workflow for prioritizing candidate invasome-associated genes from assembled *Staphylococcus* genomes.

The pipeline combines genome annotation, biologically grouped keyword screening, optional AMRFinderPlus evidence, optional ABRicate/VFDB virulence screening, protein localization and surface-association heuristics, candidate scoring, and a final summary report.


## Purpose

The pipeline helps identify candidate genes associated with:

- host interaction
- adhesion
- invasion-related annotation signals
- toxin and hemolysis activity
- immune evasion
- secretion and surface localization
- iron/heme acquisition
- stress survival
- biofilm formation
- tissue-damaging enzymes
- antimicrobial resistance context
- optional virulence database hits from ABRicate
- heuristic protein localization and surface-associated clues

The output is a prioritized candidate list for review. It is not experimental proof of invasion or virulence.

## Input Requirement

The pipeline requires an assembled *Staphylococcus* genome FASTA file.

Accepted examples:

```text
sample.fasta
sample.fa
sample.fna
sample.fas
```

The FASTA must be located inside:

```text
PROJECT_DIR/01_Genome_FASTA/
```

Raw FASTQ reads are not accepted. Read quality control, trimming, assembly, and assembly assessment should be performed before using this pipeline.

## What the Pipeline Does

1. Loads `config.env`.
2. Validates required settings and input safety rules.
3. Runs Bakta genome annotation.
4. Extracts Bakta TSV features into a clean feature table.
5. Screens gene/product/database/annotation text for candidate invasome-related keywords.
6. Runs AMRFinderPlus if the `amrfinder` command is available and Bakta protein output exists.
7. Runs ABRicate if `RUN_ABRICATE="true"`, `ABRICATE_DB` is set, and the `abricate` command is available.
8. Summarizes optional AMRFinderPlus and ABRicate evidence if present.
9. Screens proteins for lightweight localization/surface-associated clues.
10. Scores candidate genes using simple, explainable rules.
11. Generates a final candidate report.

AMRFinderPlus and ABRicate are optional evidence layers. Missing tools, failed optional runs, or missing ABRicate database configuration do not stop the full pipeline. The localization module is heuristic and does not replace SignalP, TMHMM, Phobius, InterProScan, or experimental localization.

## Environment Setup

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate StaphInvasomeMining
```

The environment file includes:

- Python 3.11
- Bakta
- ncbi-amrfinderplus
- ABRicate
- Biopython
- BLAST
- HMMER
- DIAMOND
- SeqKit

Bakta requires a local Bakta database. AMRFinderPlus requires a valid AMRFinderPlus database. ABRicate requires installed/configured ABRicate databases before database screening can work.

## Prepare the MRSA252 Test Project

This repository is meant to be reusable as a template, while real analyses run in a separate project folder. A typical layout is:

```text
~/bioinformatics/00_Pipeline_Templates/ai_staphylococcus_invasome_mining_pipeline
~/bioinformatics/09_Projects/staph_invasome_test_project
```

For the MRSA252 test run, copy the template files into the project folder while preserving the input FASTA:

```bash
mkdir -p ~/bioinformatics/09_Projects/staph_invasome_test_project/01_Genome_FASTA
cp -r run_pipeline.sh config.example.env environment.yml README.md PIPELINE_EXPLANATION.md scripts ~/bioinformatics/09_Projects/staph_invasome_test_project/
```

Place the assembled genome at:

```text
~/bioinformatics/09_Projects/staph_invasome_test_project/01_Genome_FASTA/MRSA252.fasta
```

## Prepare `config.env`

The pipeline requires a private local `config.env`. Create it from the public example, then edit it for your sample:

```bash
cp config.example.env config.env
nano config.env
```

Example MRSA252 test settings:

```bash
PROJECT_DIR="$HOME/bioinformatics/09_Projects/staph_invasome_test_project"
SAMPLE_ID="MRSA252"
GENOME_FASTA="${PROJECT_DIR}/01_Genome_FASTA/MRSA252.fasta"
THREADS="4"
BAKTA_DB="$HOME/Bioinformatics/06_Tools/bakta_db/db-light"
AMRFINDER_DB="$HOME/Bioinformatics/06_Tools/bakta_db/db-light/amrfinderplus-db/2026-05-15.1"

RUN_ABRICATE="true"
ABRICATE_DB="vfdb"
```

Required validation rules:

- `PROJECT_DIR` must exist.
- `SAMPLE_ID` may contain only letters, numbers, underscores, and hyphens.
- `THREADS` must be a positive integer.
- `GENOME_FASTA` must exist.
- `GENOME_FASTA` must be inside `PROJECT_DIR/01_Genome_FASTA` after real path resolution.
- `BAKTA_DB` must exist.

Optional ABRicate settings:

- `AMRFINDER_DB` points to an optional AMRFinderPlus database directory. Leave it empty to use AMRFinderPlus defaults.
- `RUN_ABRICATE="true"` attempts ABRicate screening.
- `RUN_ABRICATE="false"` skips ABRicate.
- `ABRICATE_DB="vfdb"` requests the ABRicate VFDB database, assuming it has been installed and configured.

If ABRicate is missing, skipped, fails, or the requested database is unavailable, the pipeline records a clear status and continues.

The pipeline never deletes input FASTA files. Rerun cleanup removes generated output folders only.

## Run

From the repository folder:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

The main entrypoint is always:

```text
run_pipeline.sh
```

## Output Folders

The current folder layout is:

```text
01_Genome_FASTA/                     input FASTA files, ignored by Git
02_Bakta_Annotation/                 Bakta annotation output
03_Invasome_Feature_Extraction/      cleaned Bakta feature table
04_Invasome_Keyword_Screening/       keyword hit table
05_AMRFinderPlus/                    optional AMRFinderPlus output
06_ABRicate/                         optional ABRicate output
07_Protein_Localization/             heuristic protein localization/surface clues
08_Invasome_Summary/                 summaries and candidate scores
09_Final_Report/                     final report
logs/                                run logs
```

Important output files:

```text
03_Invasome_Feature_Extraction/bakta_features_extracted.tsv
04_Invasome_Keyword_Screening/invasome_keyword_hits.tsv
05_AMRFinderPlus/amrfinder_<SAMPLE_ID>.tsv
06_ABRicate/abricate_<SAMPLE_ID>.tsv
08_Invasome_Summary/invasome_category_counts.tsv
08_Invasome_Summary/amrfinder_status.tsv
08_Invasome_Summary/amrfinder_summary.tsv
08_Invasome_Summary/abricate_status.tsv
08_Invasome_Summary/abricate_virulence_summary.tsv
07_Protein_Localization/protein_localization_candidates.tsv
08_Invasome_Summary/invasome_candidate_scores.tsv
09_Final_Report/invasome_mining_final_report.txt
```

## Protein Localization Stage

The `07_Protein_Localization` stage uses simple reproducible heuristics to identify proteins with candidate surface, secreted, membrane-associated, cell-wall-anchored, adhesin, capsule, biofilm, toxin, or immune-interaction clues.

It checks annotation terms and lightweight sequence features such as:

- LPXTG-like C-terminal motifs
- N-terminal signal-like hydrophobic stretches
- possible lipoprotein signal motifs
- approximate transmembrane-like hydrophobic stretches

These are support clues only. They do not replace dedicated localization tools such as SignalP, TMHMM, Phobius, InterProScan, or experimental localization.

## Candidate Scoring

Scoring is simple and explainable:

- `+2` for a keyword/category match
- `+2` extra if the category is toxin/hemolysis, adhesion, immune evasion, invasion/host interaction, or secretion/surface
- `+3` if optional AMRFinderPlus provides valid AMR, stress, virulence, or pathogen-associated context
- `+3` if optional ABRicate produces a valid database hit connected to the candidate
- `+2` for high heuristic localization priority, or `+1` for medium heuristic localization priority
- `+1` for clear LPXTG/surface/secreted localization evidence
- `+1` if the gene/product text suggests surface, secreted, membrane, LPXTG, adhesin, hemolysin, toxin, immune, capsule, biofilm, iron, heme, or stress biology

Priority levels:

- High: score `>= 6`
- Medium: score `3-5`
- Low: score `1-2`

Scores prioritize manual review. They are not measurements of invasiveness.

## Keyword Categories

The keyword screen uses these categories:

- adhesion
- invasion_host_interaction
- toxin_hemolysis
- immune_evasion
- secretion_surface
- iron_heme_acquisition
- stress_survival
- biofilm
- enzyme_tissue_damage
- antimicrobial_resistance_context

## Limitations

This pipeline identifies candidates, not confirmed mechanisms.

Important limitations:

- Keyword screening can produce false positives.
- Keyword screening can miss genes with vague, novel, or sparse annotations.
- Bakta, AMRFinderPlus, and ABRicate results depend on local database versions.
- ABRicate database hits are supporting evidence, not experimental proof.
- Heuristic localization clues are candidate support, not confirmed localization.
- Optional evidence layers may be unavailable on some systems.
- The final report does not prove host invasion, pathogenicity, or disease causation.
- Avian erythrocyte invasion, or any other specific invasion phenotype, requires experimental validation or stronger specialized evidence.

## Future Milestones

Future optional modules may include:

- SignalP-style secretion prediction
- TMHMM/Phobius-style transmembrane and topology prediction
- HMMER searches against curated protein/domain profiles
- VFDB workflows beyond ABRicate defaults
- InterProScan or eggNOG functional layers
- improved scoring with curated evidence weights
- benchmark datasets for evaluation
- protein-language-model or deep-learning scoring only after reproducible evidence tables and validation data exist
- comparative evidence across multiple *Staphylococcus* genomes

## Faculty and Portfolio Positioning

This repository is intended as a transparent teaching and portfolio pipeline. The Bash entrypoint, configuration file, TSV evidence tables, and plain-text report are designed to be readable by students and reviewable by faculty. The pipeline prioritizes reproducibility and honest candidate interpretation over unsupported claims.
