#!/usr/bin/env python3

"""Summarize optional AMRFinderPlus results for Milestone 2 evidence tables."""

import argparse
import csv
import sys
from pathlib import Path


OUTPUT_COLUMNS = [
    "feature_id",
    "gene",
    "product",
    "amrfinder_element_type",
    "amrfinder_subtype",
    "evidence_class",
    "evidence_note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize AMRFinderPlus TSV results if present."
    )
    parser.add_argument(
        "--amrfinder-tsv",
        required=True,
        help="AMRFinderPlus TSV output. It may be missing if AMRFinderPlus was skipped.",
    )
    parser.add_argument(
        "--status-file",
        required=True,
        help="TSV file describing whether AMRFinderPlus completed, failed, or was skipped.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path for summarized AMRFinderPlus evidence.",
    )
    return parser.parse_args()


def normalize_header(value: str) -> str:
    return (
        value.strip()
        .lstrip("#")
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
    )


def first_available(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value:
            return value.strip()
    return ""


def read_status(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "unknown", "AMRFinderPlus status file was not found."

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return row.get("status", "unknown"), row.get("note", "")

    return "unknown", "AMRFinderPlus status file was empty."


def classify_evidence(row: dict[str, str]) -> str:
    text = " ".join(row.values()).lower()
    if any(word in text for word in ["virulence", "toxin", "pathogen", "adhesin"]):
        return "virulence_or_pathogen_associated"
    if any(word in text for word in ["stress", "biocide", "metal", "acid", "oxidative"]):
        return "stress_or_survival_associated"
    if any(word in text for word in ["resistance", "antimicrobial", "antibiotic", "amr"]):
        return "amr_context"
    return "amrfinder_associated"


def note_row(status: str, note: str) -> dict[str, str]:
    return {
        "feature_id": "",
        "gene": "",
        "product": "",
        "amrfinder_element_type": status,
        "amrfinder_subtype": "",
        "evidence_class": "no_amrfinder_evidence",
        "evidence_note": note,
    }


def summarize_amrfinder(path: Path, status_file: Path) -> list[dict[str, str]]:
    status, status_note = read_status(status_file)
    if status != "completed":
        return [
            note_row(
                status,
                status_note
                or "AMRFinderPlus was not completed; no AMRFinderPlus evidence was used.",
            )
        ]

    if not path.exists() or path.stat().st_size == 0:
        return [note_row(status, status_note or "AMRFinderPlus output was missing or empty.")]

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [note_row(status, "AMRFinderPlus output had no header row.")]

        reader.fieldnames = [normalize_header(name) for name in reader.fieldnames]
        rows: list[dict[str, str]] = []

        for row in reader:
            normalized = {
                normalize_header(key): (value or "").strip()
                for key, value in row.items()
                if key is not None
            }
            gene = first_available(
                normalized,
                [
                    "gene_symbol",
                    "gene",
                    "element_symbol",
                    "element_name",
                    "protein_identifier",
                ],
            )
            product = first_available(
                normalized,
                [
                    "sequence_name",
                    "product",
                    "element_name",
                    "protein_name",
                    "class",
                ],
            )
            feature_id = first_available(
                normalized,
                [
                    "protein_identifier",
                    "contig_id",
                    "sequence_identifier",
                    "accession_of_closest_sequence",
                    "gene_symbol",
                ],
            )
            element_type = first_available(
                normalized,
                ["element_type", "type", "class", "subclass", "method"],
            )
            subtype = first_available(
                normalized,
                ["element_subtype", "subtype", "subclass", "scope"],
            )
            evidence_class = classify_evidence(normalized)
            rows.append(
                {
                    "feature_id": feature_id,
                    "gene": gene,
                    "product": product,
                    "amrfinder_element_type": element_type,
                    "amrfinder_subtype": subtype,
                    "evidence_class": evidence_class,
                    "evidence_note": "Optional AMRFinderPlus evidence; interpret as supporting context, not proof of invasion.",
                }
            )

    if not rows:
        return [note_row(status, "AMRFinderPlus output contained no result rows.")]
    return rows


def write_summary(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        rows = summarize_amrfinder(Path(args.amrfinder_tsv), Path(args.status_file))
        write_summary(rows, Path(args.output))
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"AMRFinderPlus summary rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
