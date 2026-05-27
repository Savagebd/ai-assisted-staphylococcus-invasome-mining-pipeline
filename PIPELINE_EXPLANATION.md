# Pipeline Explanation

## What Is an Invasome?

In this project, "invasome" is a practical umbrella term. It refers to genes and proteins that may help a bacterium interact with host tissue, attach to cells or extracellular matrix, damage tissue, evade immune responses, survive stress, acquire nutrients, form biofilms, or express virulence-associated traits.

It is not one single pathway. For *Staphylococcus*, invasion and host interaction can involve many different systems working together.

## Why *Staphylococcus*?

*Staphylococcus* species are important in clinical, veterinary, food, and environmental microbiology. Some strains are harmless colonizers, while others can cause skin infections, bloodstream infections, device-associated biofilms, pneumonia, toxin-mediated disease, or other severe outcomes.

Candidate invasome-related genes can include:

- adhesins that bind host molecules
- surface proteins that contact host tissue
- toxins and hemolysins
- immune evasion proteins
- iron and heme acquisition systems
- stress survival systems
- biofilm-associated genes
- enzymes that damage host tissue
- antimicrobial resistance context genes that may support survival during treatment

## Why Use an Assembled Genome FASTA?

The pipeline starts from an assembled genome FASTA because Bakta and ABRicate screen assembled contigs or complete genomes.

Raw FASTQ files need separate preprocessing:

- quality control
- trimming
- genome assembly
- assembly quality assessment

Those steps are important, but they are outside the current pipeline scope.

## Why Bakta?

Bakta is used because it provides reproducible bacterial genome annotation. It predicts genomic features and assigns gene names, product descriptions, locus tags, and database cross-references where possible.

For this pipeline, the most important Bakta output is the TSV annotation file. It gives structured text that can be searched and summarized.

The extracted feature table keeps fields such as:

- sample ID
- feature ID
- contig
- start and stop positions
- strand
- feature type
- gene
- product
- database cross-references
- estimated nucleotide and amino acid lengths
- raw annotation text

## Why Keyword Screening?

Keyword screening is a transparent first-pass method. It is easy for students and reviewers to understand because each hit is connected to visible annotation text.

For example:

- `fibronectin`, `fibrinogen`, or `adhesin` may suggest adhesion-related candidates.
- `hemolysin`, `leukocidin`, or `toxin` may suggest toxin or hemolysis candidates.
- `sortase`, `LPXTG`, `surface`, or `secreted` may suggest secretion or surface localization context.
- `iron`, `heme`, or `isd` may suggest nutrient acquisition during host association.
- `biofilm`, `ica`, or `autolysin` may suggest biofilm-related biology.

Keyword screening is useful, but it is imperfect. It can miss genes with unusual names and can flag genes that are not truly involved in invasion.

## Why AMRFinderPlus Is Optional

AMRFinderPlus can add useful context, especially for antimicrobial resistance genes and some virulence or pathogen-associated elements. However, not every user will have AMRFinderPlus installed or configured.

The pipeline therefore treats AMRFinderPlus as optional:

- If `amrfinder` is available and Bakta protein output exists, the pipeline runs it.
- If `AMRFINDER_DB` is configured, the pipeline passes that database path to AMRFinderPlus.
- If it is missing, the pipeline records that it was skipped and continues.
- If it fails, the pipeline records that failure and continues with Bakta keyword evidence.

This keeps the template usable for beginners while still supporting richer evidence when available.

## What ABRicate Does

ABRicate screens assembled genome sequences against nucleotide databases. Depending on which databases are installed, it can be used for antimicrobial resistance genes, plasmid markers, or virulence-associated databases such as VFDB-style resources.

In this pipeline, ABRicate is used as an optional database evidence layer for candidate invasome or virulence screening.

## Why ABRicate Helps for Virulence and Invasome Mining

Bakta annotation and keyword screening ask:

```text
What does the genome annotation say this feature might be?
```

ABRicate asks a different question:

```text
Does this assembled genome contain sequence similarity to entries in a selected database?
```

That can help find known virulence-associated genes even when ordinary annotation text is incomplete. For example, screening with an installed VFDB-style database may highlight genes with similarity to known virulence factors.

## Why ABRicate Is Optional

ABRicate depends on local installation and local database setup. Public GitHub users may not have the same databases installed, and some databases need separate download or setup steps.

The pipeline therefore uses these rules:

- `RUN_ABRICATE="true"` asks the pipeline to try ABRicate.
- `RUN_ABRICATE="false"` skips ABRicate.
- `ABRICATE_DB` names the requested ABRicate database, such as `vfdb`.
- Missing ABRicate, missing databases, skipped runs, and failed optional runs do not stop the pipeline.

This design keeps the baseline reproducible while allowing stronger evidence when a user has ABRicate configured.

## Why Database Screening Is Not Proof

An ABRicate hit is useful evidence, but it is not experimental proof of invasion.

Reasons include:

- database entries can vary in quality and scope
- sequence similarity does not guarantee expression or function
- a gene may be fragmented, pseudogenized, or context-dependent
- virulence depends on host environment, regulation, and strain background
- database hits can require manual confirmation of identity and coverage

ABRicate hits should be interpreted as supporting evidence for prioritization.

## Why Protein Localization Matters for Invasome Mining

Many invasome-relevant proteins act at the bacterial surface or outside the cell. Examples include adhesins, secreted toxins, immune interaction proteins, capsule-associated proteins, biofilm proteins, cell-wall anchored proteins, and membrane-associated transport or stress-survival systems.

For Gram-positive bacteria such as *Staphylococcus*, some surface proteins are attached to the cell wall by sortase enzymes. A common clue is an LPXTG-like motif near the C-terminal end of a protein. This motif does not prove surface localization by itself, but it is a useful screening clue when combined with annotation and database evidence.

The current protein localization module uses simple reproducible heuristics:

- annotation terms such as surface, secreted, cell wall, membrane, adhesin, capsule, biofilm, toxin, or immune interaction
- LPXTG-like motifs near the C-terminal region
- N-terminal signal-like hydrophobic stretches
- possible lipoprotein signal motifs
- approximate transmembrane-like hydrophobic stretches

These are localization/support clues only. They do not replace SignalP, TMHMM, Phobius, InterProScan, or experimental localization assays.

## Why Candidate Scoring?

The scoring system helps prioritize manual review. It gives higher scores to candidates with stronger or more relevant evidence categories.

The current scoring is intentionally simple:

- keyword/category evidence gives a baseline score
- categories closely related to host interaction receive extra weight
- AMRFinderPlus AMR, stress, virulence, or pathogen-associated context adds support when available
- ABRicate database hits add support when they can be connected to a candidate
- heuristic protein localization and surface-associated evidence adds support when available
- product text suggesting surface, secreted, membrane, toxin, immune, capsule, biofilm, iron/heme, or stress biology adds one small supporting point

This is not a predictive model. It is an evidence checklist converted into a sortable table.

## Why Scoring Does Not Prove Invasion

A score is not experimental proof.

A high-priority candidate means:

```text
This gene has annotation and optional database evidence that makes it worth reviewing first.
```

It does not mean:

```text
This gene is proven to cause host invasion.
```

Reasons include:

- annotations can be wrong or incomplete
- genes may be fragmented or not expressed
- sequence presence does not guarantee phenotype
- broad keywords can create false positives
- true invasion depends on host, strain background, regulation, and experimental context

Manual curation and laboratory validation are still required.

## What the Final Report Provides

The final report summarizes:

- sample ID
- input FASTA path
- number of extracted Bakta features
- number of keyword hits
- number of scored candidates
- high/medium/low priority counts
- keyword category counts
- AMRFinderPlus status
- ABRicate status
- total ABRicate hit rows
- protein localization candidate counts
- top candidate genes
- interpretation notes
- limitations
- future module ideas

The report is designed for quick review and teaching, while the TSV tables remain the primary reproducible evidence files.

## What Future Deep Learning Could Add

Deep learning could eventually help with:

- ranking candidate virulence genes from combined evidence
- learning patterns from curated invasive and non-invasive strain datasets
- integrating protein sequence, domains, localization, and annotation text
- generating hypothesis summaries from structured evidence tables

However, deep learning should not be added before the baseline evidence workflow is reproducible and testable. A model is only useful if the input labels, benchmark data, and evaluation rules are clear.

## Future Modules

Future optional modules may include:

- SignalP-style secretion prediction
- TMHMM-style transmembrane prediction
- HMMER domain/profile searches
- VFDB workflows beyond ABRicate defaults
- InterProScan or eggNOG annotation layers
- curated virulence gene review tables
- improved scoring with transparent weights
- optional AI-assisted interpretation after evidence tables are generated
