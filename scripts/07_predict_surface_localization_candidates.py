#!/usr/bin/env python3

"""Predict heuristic surface/localization support clues for Bakta proteins.

This module deliberately uses simple, explainable rules. It does not replace
SignalP, TMHMM, Phobius, InterProScan, or curated experimental localization.
"""

import argparse
import csv
import re
import sys
from pathlib import Path


OUTPUT_COLUMNS = [
    "sample_id",
    "feature_id",
    "gene",
    "product",
    "protein_length",
    "localization_clues",
    "matched_localization_terms",
    "sequence_motif_clues",
    "localization_score",
    "localization_priority",
    "evidence_note",
]


LOCALIZATION_TERMS = [
    "surface",
    "cell wall",
    "cell-wall",
    "membrane",
    "lipoprotein",
    "adhesin",
    "adhesion",
    "mscramm",
    "sortase",
    "lpxtg",
    "anchored",
    "extracellular",
    "secreted",
    "secretion",
    "capsule",
    "capsular",
    "biofilm",
    "fibrinogen",
    "fibronectin",
    "collagen",
    "immunoglobulin",
    "hemolysin",
    "haemolysin",
    "toxin",
    "leukocidin",
    "nuclease",
    "protease",
]


HYDROPHOBIC = set("AILMFWVY")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find heuristic protein localization and surface-associated candidate clues."
    )
    parser.add_argument("--features", required=True, help="Extracted Bakta feature TSV.")
    parser.add_argument(
        "--protein-fasta",
        required=False,
        default="",
        help="Optional Bakta protein FASTA file.",
    )
    parser.add_argument("--output", required=True, help="Localization candidate TSV.")
    return parser.parse_args()


def read_features(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Feature table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Feature table has no header row: {path}")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def fasta_header_ids(header: str) -> list[str]:
    primary = header.split()[0]
    ids = [primary]
    for pattern in [r"locus_tag=([^\]\s]+)", r"ID=([^;\s]+)", r"Name=([^;\s]+)"]:
        match = re.search(pattern, header)
        if match:
            ids.append(match.group(1))
    return sorted(set(ids))


def read_fasta(path_text: str) -> dict[str, str]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists() or not path.is_file():
        return {}

    sequences: dict[str, str] = {}
    current_header = ""
    chunks: list[str] = []

    def store_record() -> None:
        if not current_header:
            return
        sequence = "".join(chunks).upper()
        for identifier in fasta_header_ids(current_header):
            sequences[identifier] = sequence

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store_record()
                current_header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        store_record()

    return sequences


def feature_sequence(feature: dict[str, str], sequences: dict[str, str]) -> str:
    for key in [
        feature.get("feature_id", ""),
        feature.get("gene", ""),
        feature.get("locus_tag", ""),
    ]:
        if key and key in sequences:
            return sequences[key]
    return ""


def hydrophobic_count(window: str) -> int:
    return sum(1 for aa in window if aa in HYDROPHOBIC)


def has_signal_like_n_terminus(sequence: str) -> bool:
    n_term = sequence[:35]
    if len(n_term) < 18:
        return False
    for index in range(0, max(len(n_term) - 8, 1)):
        window = n_term[index : index + 9]
        if len(window) == 9 and hydrophobic_count(window) >= 7:
            return True
    return False


def has_lipoprotein_like_motif(sequence: str) -> bool:
    n_term = sequence[:35]
    return bool(re.search(r"[A-Z]{1,20}[LVI][ASTVI][GAS]C", n_term))


def has_lpxtg_like_motif(sequence: str) -> bool:
    c_term = sequence[-80:]
    return bool(re.search(r"LP.TG", c_term))


def count_tm_like_stretches(sequence: str) -> int:
    count = 0
    index = 0
    while index <= len(sequence) - 18:
        window = sequence[index : index + 18]
        if hydrophobic_count(window) >= 14:
            count += 1
            index += 17
        else:
            index += 1
    return min(count, 10)


def matched_terms(feature: dict[str, str]) -> list[str]:
    text = " ".join(
        [
            feature.get("gene", ""),
            feature.get("product", ""),
            feature.get("db_xrefs", ""),
            feature.get("raw_annotation_text", ""),
        ]
    ).lower()
    return sorted({term for term in LOCALIZATION_TERMS if term in text})


def priority(score: int) -> str:
    if score >= 4:
        return "High"
    if score >= 2:
        return "Medium"
    if score >= 1:
        return "Low"
    return "Unprioritized"


def score_feature(feature: dict[str, str], sequence: str) -> dict[str, str] | None:
    terms = matched_terms(feature)
    motif_clues: list[str] = []
    text_clues: list[str] = []
    score = 0

    if terms:
        text_clues.append("annotation_localization_terms")
        score += 2

    if sequence:
        if has_lpxtg_like_motif(sequence):
            motif_clues.append("c_terminal_lpxtg_like_motif")
            score += 2
        if has_signal_like_n_terminus(sequence):
            motif_clues.append("n_terminal_signal_like_hydrophobic_stretch")
            score += 1
        if has_lipoprotein_like_motif(sequence):
            motif_clues.append("possible_lipoprotein_signal_motif")
            score += 1
        tm_count = count_tm_like_stretches(sequence)
        if tm_count >= 2:
            motif_clues.append(f"multiple_tm_like_hydrophobic_stretches:{tm_count}")
            score += 1
        elif tm_count == 1:
            motif_clues.append("single_tm_like_hydrophobic_stretch")

    if score <= 0:
        return None

    protein_length = str(len(sequence)) if sequence else feature.get("aa_sequence_length", "")
    clues = sorted(set(text_clues + motif_clues))
    return {
        "sample_id": feature.get("sample_id", ""),
        "feature_id": feature.get("feature_id", ""),
        "gene": feature.get("gene", ""),
        "product": feature.get("product", ""),
        "protein_length": protein_length,
        "localization_clues": ";".join(clues),
        "matched_localization_terms": ";".join(terms),
        "sequence_motif_clues": ";".join(motif_clues),
        "localization_score": str(score),
        "localization_priority": priority(score),
        "evidence_note": "Heuristic localization/support clues only; not confirmed localization.",
    }


def write_rows(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        features = read_features(Path(args.features))
        sequences = read_fasta(args.protein_fasta)
        rows = []
        for feature in features:
            row = score_feature(feature, feature_sequence(feature, sequences))
            if row is not None:
                rows.append(row)
        write_rows(rows, Path(args.output))
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Protein localization candidate rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
