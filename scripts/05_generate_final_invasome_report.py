#!/usr/bin/env python3

"""Generate the Milestone 2 final candidate invasome report."""

import argparse
import csv
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the final Milestone 2 invasome mining report."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--genome-fasta", required=True, help="Input genome FASTA path.")
    parser.add_argument("--features", required=True, help="Extracted Bakta feature TSV.")
    parser.add_argument("--keyword-hits", required=True, help="Keyword hits TSV.")
    parser.add_argument("--category-counts", required=True, help="Category counts TSV.")
    parser.add_argument("--amrfinder-status", required=True, help="AMRFinderPlus status TSV.")
    parser.add_argument("--amrfinder-summary", required=True, help="AMRFinderPlus summary TSV.")
    parser.add_argument("--abricate-status", required=True, help="ABRicate status TSV.")
    parser.add_argument("--abricate-summary", required=True, help="ABRicate summary TSV.")
    parser.add_argument("--localization-summary", required=True, help="Heuristic protein localization candidate TSV.")
    parser.add_argument("--candidate-scores", required=True, help="Candidate scores TSV.")
    parser.add_argument("--output", required=True, help="Final report output path.")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Required report input not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header row: {path}")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def read_category_counts(path: Path) -> list[tuple[str, str]]:
    rows = read_tsv(path)
    output: list[tuple[str, str]] = []
    for row in rows:
        category = row.get("category", "")
        count = row.get("candidate_hit_count", row.get("hit_count", "0"))
        output.append((category, count))
    return output


def read_status(path: Path, tool_name: str) -> tuple[str, str]:
    rows = read_tsv(path)
    if not rows:
        return "unknown", f"{tool_name} status file was empty."
    return rows[0].get("status", "unknown"), rows[0].get("note", "")


def priority_counts(candidate_rows: list[dict[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in candidate_rows:
        counts[row.get("priority_level", "Unprioritized")] += 1
    return counts


def safe_int(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return 0


def top_candidates(candidate_rows: list[dict[str, str]], limit: int = 15) -> list[dict[str, str]]:
    return sorted(
        candidate_rows,
        key=lambda row: (
            safe_int(row.get("candidate_score", "0")),
            row.get("priority_level", ""),
            row.get("gene", ""),
        ),
        reverse=True,
    )[:limit]


def count_valid_abricate_hits(rows: list[dict[str, str]]) -> int:
    return sum(
        1
        for row in rows
        if row.get("evidence_source", "") == "optional_abricate_database_hit"
    )


def count_valid_amrfinder_hits(rows: list[dict[str, str]]) -> int:
    return sum(
        1
        for row in rows
        if row.get("evidence_class", "") != "no_amrfinder_evidence"
    )


def count_candidates_with_amrfinder(rows: list[dict[str, str]]) -> int:
    return sum(1 for row in rows if row.get("amrfinder_evidence", "none") != "none")


def count_candidates_with_abricate(rows: list[dict[str, str]]) -> int:
    return sum(1 for row in rows if row.get("abricate_evidence", "none") != "none")


def localization_priority_counts(rows: list[dict[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in rows:
        counts[row.get("localization_priority", "Unprioritized")] += 1
    return counts


def count_candidates_with_localization(rows: list[dict[str, str]]) -> int:
    return sum(1 for row in rows if row.get("localization_evidence", "none") != "none")


def clipped(value: str, width: int) -> str:
    if len(value) <= width:
        return value
    return value[: max(width - 3, 0)] + "..."


def write_candidate_table(handle, rows: list[dict[str, str]]) -> None:
    if not rows:
        handle.write("No candidates were scored in this run.\n")
        return

    header = f"{'Priority':<8}  {'Score':<5}  {'Gene':<16}  {'Feature ID':<22}  {'Categories':<42}  Product\n"
    handle.write(header)
    handle.write("-" * 120 + "\n")
    for row in rows:
        handle.write(
            f"{clipped(row.get('priority_level', ''), 8):<8}  "
            f"{clipped(row.get('candidate_score', ''), 5):<5}  "
            f"{clipped(row.get('gene', '') or 'NA', 16):<16}  "
            f"{clipped(row.get('feature_id', '') or 'NA', 22):<22}  "
            f"{clipped(row.get('matched_categories', ''), 42):<42}  "
            f"{clipped(row.get('product', '') or 'NA', 80)}\n"
        )


def write_report(args: argparse.Namespace) -> None:
    features = read_tsv(Path(args.features))
    keyword_hits = read_tsv(Path(args.keyword_hits))
    category_counts = read_category_counts(Path(args.category_counts))
    amrfinder_status, amrfinder_note = read_status(
        Path(args.amrfinder_status), "AMRFinderPlus"
    )
    amrfinder_rows = read_tsv(Path(args.amrfinder_summary))
    abricate_status, abricate_note = read_status(Path(args.abricate_status), "ABRicate")
    abricate_rows = read_tsv(Path(args.abricate_summary))
    localization_rows = read_tsv(Path(args.localization_summary))
    candidates = read_tsv(Path(args.candidate_scores))
    counts = priority_counts(candidates)
    loc_counts = localization_priority_counts(localization_rows)
    total_amrfinder_hits = count_valid_amrfinder_hits(amrfinder_rows)
    candidates_with_amrfinder = count_candidates_with_amrfinder(candidates)
    total_abricate_hits = count_valid_abricate_hits(abricate_rows)
    candidates_with_abricate = count_candidates_with_abricate(candidates)
    candidates_with_localization = count_candidates_with_localization(candidates)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("AI-Assisted Staphylococcus Invasome Mining Pipeline\n")
        handle.write("Candidate Invasome Gene Report\n")
        handle.write("=" * 64 + "\n\n")
        handle.write(f"Sample ID: {args.sample_id}\n")
        handle.write(f"Input genome FASTA: {args.genome_fasta}\n")
        handle.write(f"Report generated: {datetime.now().isoformat(timespec='seconds')}\n\n")

        handle.write("Run summary\n")
        handle.write(f"- Total extracted Bakta features: {len(features)}\n")
        handle.write(f"- Total keyword hit records: {len(keyword_hits)}\n")
        handle.write(f"- Total scored candidates: {len(candidates)}\n")
        handle.write(f"- High priority candidates: {counts.get('High', 0)}\n")
        handle.write(f"- Medium priority candidates: {counts.get('Medium', 0)}\n")
        handle.write(f"- Low priority candidates: {counts.get('Low', 0)}\n")
        handle.write(f"- AMRFinderPlus status: {amrfinder_status}\n")
        handle.write(f"- AMRFinderPlus note: {amrfinder_note}\n")
        handle.write(f"- Total AMRFinderPlus summary hit rows: {total_amrfinder_hits}\n")
        handle.write(
            f"- Candidates with AMRFinderPlus evidence contributing to scoring: {candidates_with_amrfinder}\n\n"
        )
        handle.write(f"- ABRicate status: {abricate_status}\n")
        handle.write(f"- ABRicate note: {abricate_note}\n")
        handle.write(f"- Total ABRicate database hit rows: {total_abricate_hits}\n")
        handle.write(
            f"- Candidates with ABRicate evidence contributing to scoring: {candidates_with_abricate}\n\n"
        )
        handle.write("- Protein localization status: completed heuristic screening\n")
        handle.write(f"- Total protein localization candidate rows: {len(localization_rows)}\n")
        handle.write(f"- High localization candidates: {loc_counts.get('High', 0)}\n")
        handle.write(f"- Medium localization candidates: {loc_counts.get('Medium', 0)}\n")
        handle.write(f"- Low localization candidates: {loc_counts.get('Low', 0)}\n")
        handle.write(
            f"- Candidates with localization evidence contributing to scoring: {candidates_with_localization}\n\n"
        )

        handle.write("Candidate category counts\n")
        for category, count in category_counts:
            handle.write(f"- {category}: {count}\n")
        handle.write("\n")

        handle.write("Top candidate genes\n")
        write_candidate_table(handle, top_candidates(candidates))
        handle.write("\n")

        handle.write("Interpretation\n")
        handle.write(
            "The candidates listed above are prioritized from one or more evidence layers: "
            "Bakta annotation keyword/category matches, optional AMRFinderPlus context, optional "
            "ABRicate database hits, and heuristic protein localization/surface clues. These evidence layers can highlight genes associated with adhesion, "
            "host interaction, toxin or hemolysis, immune evasion, secretion or surface localization, "
            "iron/heme acquisition, stress survival, biofilm biology, tissue-damaging enzymes, "
            "antimicrobial resistance context, or curated virulence database similarity.\n\n"
        )
        handle.write(
            "Scores are intended to help decide which genes deserve manual review first. "
            "They are not measurements of invasiveness and do not prove that a strain invades host cells.\n\n"
        )

        handle.write("Limitations\n")
        handle.write("- Candidate wording is deliberate: this is evidence prioritization, not experimental proof.\n")
        handle.write("- Keyword screening can produce false positives when broad terms appear in unrelated annotations.\n")
        handle.write("- Keyword screening can miss genes with sparse, novel, or unusual annotation text.\n")
        handle.write("- Bakta, AMRFinderPlus, and ABRicate results depend on genome assembly quality and local database versions.\n")
        handle.write("- AMRFinderPlus is optional; skipped or failed runs can still produce candidate scores from other evidence layers.\n")
        handle.write("- ABRicate is optional database screening; database hits are supporting evidence, not experimental proof of invasion.\n")
        handle.write("- Protein localization screening is heuristic and does not replace SignalP, TMHMM, Phobius, or experimental localization.\n")
        handle.write("- Avian erythrocyte invasion would require experimental validation or stronger specialized evidence.\n\n")

        handle.write("Future modules\n")
        handle.write("- SignalP-style secretion prediction.\n")
        handle.write("- TMHMM-style transmembrane prediction.\n")
        handle.write("- HMMER profile searches for curated adhesin, secretion, and virulence-associated domains.\n")
        handle.write("- VFDB, InterProScan, eggNOG, or other curated annotation layers.\n")
        handle.write("- Protein-language-model or deep-learning scoring after reproducible evidence tables and benchmark data are available.\n")
        handle.write("- Comparative evidence across multiple Staphylococcus genomes.\n\n")

        handle.write("Main output files\n")
        handle.write(f"- Extracted Bakta features: {args.features}\n")
        handle.write(f"- Keyword hits: {args.keyword_hits}\n")
        handle.write(f"- Category counts: {args.category_counts}\n")
        handle.write(f"- AMRFinderPlus summary: {args.amrfinder_summary}\n")
        handle.write(f"- ABRicate summary: {args.abricate_summary}\n")
        handle.write(f"- Protein localization candidates: {args.localization_summary}\n")
        handle.write(f"- Candidate scores: {args.candidate_scores}\n")


def main() -> int:
    args = parse_args()
    try:
        write_report(args)
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Final report written to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
