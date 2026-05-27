#!/usr/bin/env python3

"""Score candidate invasome-related genes using simple explainable evidence."""

import argparse
import csv
import sys
from pathlib import Path


OUTPUT_COLUMNS = [
    "sample_id",
    "feature_id",
    "gene",
    "product",
    "matched_categories",
    "matched_keywords",
    "amrfinder_evidence",
    "abricate_evidence",
    "localization_evidence",
    "candidate_score",
    "priority_level",
    "interpretation_note",
]


HIGH_WEIGHT_CATEGORIES = {
    "toxin_hemolysis",
    "adhesion",
    "immune_evasion",
    "invasion_host_interaction",
    "secretion_surface",
}


SUPPORTING_TERMS = [
    "surface",
    "secreted",
    "secretion",
    "membrane",
    "lpxtg",
    "adhesin",
    "adhesion",
    "hemolysin",
    "haemolysin",
    "toxin",
    "immune",
    "capsule",
    "capsular",
    "biofilm",
    "iron",
    "heme",
    "haem",
    "stress",
    "virulence",
    "vfdb",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score candidate invasome-related genes from keyword, AMRFinderPlus, ABRicate, and localization evidence."
    )
    parser.add_argument("--keyword-hits", required=True, help="Keyword hit TSV.")
    parser.add_argument("--amrfinder-summary", required=False, help="Optional summarized AMRFinderPlus TSV.")
    parser.add_argument("--abricate-summary", required=False, help="Optional summarized ABRicate virulence/database TSV.")
    parser.add_argument("--localization-summary", required=False, help="Optional heuristic protein localization TSV.")
    parser.add_argument("--output", required=True, help="Candidate score TSV output.")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Required input table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header row: {path}")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def optional_tsv(path_text: str | None) -> list[dict[str, str]]:
    if not path_text:
        return []
    path = Path(path_text)
    if not path.exists():
        return []
    return read_tsv(path)


def norm(value: str) -> str:
    return value.strip().lower()


def split_semicolon(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def merge_semicolon(existing: str, additions: list[str]) -> str:
    merged = split_semicolon(existing)
    merged.extend(item for item in additions if item and item != "none")
    return ";".join(sorted(set(merged)))


def priority_level(score: int) -> str:
    if score >= 6:
        return "High"
    if score >= 3:
        return "Medium"
    if score >= 1:
        return "Low"
    return "Unprioritized"


def candidate_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    feature_id = norm(row.get("feature_id", ""))
    gene = norm(row.get("gene", ""))
    product = norm(row.get("product", ""))
    accession = norm(row.get("accession", ""))
    location = "|".join([norm(row.get("sequence", "")), norm(row.get("start", "")), norm(row.get("end", ""))])

    if feature_id:
        return ("feature", feature_id, gene, product, "")
    if gene:
        return ("gene", gene, accession, location, product)
    if accession:
        return ("accession", accession, location, product, "")
    if product:
        return ("product", product, location, "", "")
    return ("row", location, "", "", "")


def empty_candidate() -> dict[str, str]:
    return {
        "sample_id": "",
        "feature_id": "",
        "gene": "",
        "product": "",
        "matched_categories": "",
        "matched_keywords": "",
        "amrfinder_evidence": "none",
        "abricate_evidence": "none",
        "localization_evidence": "none",
    }


def ensure_candidate(candidates: dict[tuple[str, str, str, str, str], dict[str, str]], key: tuple[str, str, str, str, str]) -> dict[str, str]:
    if key not in candidates:
        candidates[key] = empty_candidate()
    return candidates[key]


def fill_if_blank(candidate: dict[str, str], field: str, value: str) -> None:
    if value and not candidate.get(field, ""):
        candidate[field] = value


def add_keyword_candidates(candidates: dict[tuple[str, str, str, str, str], dict[str, str]], rows: list[dict[str, str]]) -> None:
    for row in rows:
        candidate = ensure_candidate(candidates, candidate_key(row))
        fill_if_blank(candidate, "sample_id", row.get("sample_id", ""))
        fill_if_blank(candidate, "feature_id", row.get("feature_id", ""))
        fill_if_blank(candidate, "gene", row.get("gene", ""))
        fill_if_blank(candidate, "product", row.get("product", ""))
        candidate["matched_categories"] = merge_semicolon(candidate["matched_categories"], split_semicolon(row.get("matched_categories", "")))
        candidate["matched_keywords"] = merge_semicolon(candidate["matched_keywords"], split_semicolon(row.get("matched_keywords", "")))


def amrfinder_note(row: dict[str, str]) -> str:
    parts = [row.get("gene", ""), row.get("product", ""), row.get("amrfinder_element_type", ""), row.get("amrfinder_subtype", ""), row.get("evidence_class", "")]
    return " | ".join(part for part in parts if part)


def add_amrfinder_candidates(candidates: dict[tuple[str, str, str, str, str], dict[str, str]], rows: list[dict[str, str]]) -> None:
    for row in rows:
        if row.get("evidence_class", "") == "no_amrfinder_evidence":
            continue
        candidate = ensure_candidate(candidates, candidate_key(row))
        fill_if_blank(candidate, "feature_id", row.get("feature_id", ""))
        fill_if_blank(candidate, "gene", row.get("gene", ""))
        fill_if_blank(candidate, "product", row.get("product", ""))
        candidate["amrfinder_evidence"] = merge_semicolon(
            "" if candidate["amrfinder_evidence"] == "none" else candidate["amrfinder_evidence"],
            [amrfinder_note(row)],
        ) or "none"


def abricate_note(row: dict[str, str]) -> str:
    parts = [
        row.get("gene", ""),
        row.get("database", ""),
        row.get("accession", ""),
        row.get("product", ""),
        row.get("sequence", ""),
        row.get("start", ""),
        row.get("end", ""),
        row.get("percent_identity", ""),
        row.get("percent_coverage", ""),
    ]
    return " | ".join(part for part in parts if part)


def add_abricate_candidates(candidates: dict[tuple[str, str, str, str, str], dict[str, str]], rows: list[dict[str, str]]) -> None:
    for row in rows:
        if row.get("evidence_source", "") != "optional_abricate_database_hit":
            continue
        candidate = ensure_candidate(candidates, candidate_key(row))
        fill_if_blank(candidate, "sample_id", row.get("sample_id", ""))
        fill_if_blank(candidate, "gene", row.get("gene", ""))
        fill_if_blank(candidate, "product", row.get("product", ""))
        candidate["abricate_evidence"] = merge_semicolon(
            "" if candidate["abricate_evidence"] == "none" else candidate["abricate_evidence"],
            [abricate_note(row)],
        ) or "none"


def localization_note(row: dict[str, str]) -> str:
    parts = [
        row.get("localization_priority", ""),
        row.get("localization_score", ""),
        row.get("localization_clues", ""),
        row.get("matched_localization_terms", ""),
        row.get("sequence_motif_clues", ""),
    ]
    return " | ".join(part for part in parts if part)


def add_localization_candidates(candidates: dict[tuple[str, str, str, str, str], dict[str, str]], rows: list[dict[str, str]]) -> None:
    for row in rows:
        if not row.get("localization_score", ""):
            continue
        candidate = ensure_candidate(candidates, candidate_key(row))
        fill_if_blank(candidate, "sample_id", row.get("sample_id", ""))
        fill_if_blank(candidate, "feature_id", row.get("feature_id", ""))
        fill_if_blank(candidate, "gene", row.get("gene", ""))
        fill_if_blank(candidate, "product", row.get("product", ""))
        candidate["localization_evidence"] = merge_semicolon(
            "" if candidate["localization_evidence"] == "none" else candidate["localization_evidence"],
            [localization_note(row)],
        ) or "none"


def has_evidence(candidate: dict[str, str], field: str) -> bool:
    return candidate.get(field, "none") != "none"


def localization_score_bonus(candidate: dict[str, str]) -> int:
    evidence = candidate.get("localization_evidence", "")
    if "High" in evidence:
        bonus = 2
    elif "Medium" in evidence:
        bonus = 1
    else:
        bonus = 0
    if any(term in evidence.lower() for term in ["lpxtg", "surface", "secreted", "cell wall", "cell-wall"]):
        bonus += 1
    return bonus


def score_candidate(candidate: dict[str, str]) -> int:
    score = 0
    categories = split_semicolon(candidate.get("matched_categories", ""))
    text = " ".join(
        [
            candidate.get("gene", ""),
            candidate.get("product", ""),
            candidate.get("matched_keywords", ""),
            candidate.get("amrfinder_evidence", ""),
            candidate.get("abricate_evidence", ""),
            candidate.get("localization_evidence", ""),
        ]
    ).lower()

    if categories:
        score += 2
    if any(category in HIGH_WEIGHT_CATEGORIES for category in categories):
        score += 2
    if has_evidence(candidate, "amrfinder_evidence"):
        score += 3
    if has_evidence(candidate, "abricate_evidence"):
        score += 3
    score += localization_score_bonus(candidate)
    if any(term in text for term in SUPPORTING_TERMS):
        score += 1
    return score


def interpretation(candidate: dict[str, str], score: int) -> str:
    level = priority_level(score)
    layers: list[str] = []
    if candidate.get("matched_categories", ""):
        layers.append(f"keyword categories ({candidate['matched_categories']})")
    if has_evidence(candidate, "amrfinder_evidence"):
        layers.append("optional AMRFinderPlus context")
    if has_evidence(candidate, "abricate_evidence"):
        layers.append("optional ABRicate database evidence")
    if has_evidence(candidate, "localization_evidence"):
        layers.append("heuristic protein localization/surface clues")
    if not layers:
        layers.append("limited supporting text evidence")
    return f"{level} priority candidate based on {' and '.join(layers)}. This is evidence prioritization, not proof of invasion."


def score_candidates(keyword_rows: list[dict[str, str]], amrfinder_rows: list[dict[str, str]], abricate_rows: list[dict[str, str]], localization_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    candidates: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    add_keyword_candidates(candidates, keyword_rows)
    add_amrfinder_candidates(candidates, amrfinder_rows)
    add_abricate_candidates(candidates, abricate_rows)
    add_localization_candidates(candidates, localization_rows)

    scored: list[dict[str, str]] = []
    for candidate in candidates.values():
        score = score_candidate(candidate)
        scored.append(
            {
                "sample_id": candidate.get("sample_id", ""),
                "feature_id": candidate.get("feature_id", ""),
                "gene": candidate.get("gene", ""),
                "product": candidate.get("product", ""),
                "matched_categories": candidate.get("matched_categories", ""),
                "matched_keywords": candidate.get("matched_keywords", ""),
                "amrfinder_evidence": candidate.get("amrfinder_evidence", "none"),
                "abricate_evidence": candidate.get("abricate_evidence", "none"),
                "localization_evidence": candidate.get("localization_evidence", "none"),
                "candidate_score": str(score),
                "priority_level": priority_level(score),
                "interpretation_note": interpretation(candidate, score),
            }
        )
    scored.sort(key=lambda item: (int(item["candidate_score"]), item.get("priority_level", ""), item.get("gene", ""), item.get("product", "")), reverse=True)
    return scored


def write_scores(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        scored = score_candidates(
            read_tsv(Path(args.keyword_hits)),
            optional_tsv(args.amrfinder_summary),
            optional_tsv(args.abricate_summary),
            optional_tsv(args.localization_summary),
        )
        write_scores(scored, Path(args.output))
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Scored candidate records: {len(scored)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
