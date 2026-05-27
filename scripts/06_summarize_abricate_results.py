#!/usr/bin/env python3

"""Summarize optional ABRicate virulence/database screening results."""

import argparse
import csv
import sys
from pathlib import Path


OUTPUT_COLUMNS = [
    "sample_id",
    "gene",
    "database",
    "accession",
    "percent_identity",
    "percent_coverage",
    "sequence",
    "start",
    "end",
    "product",
    "evidence_source",
    "evidence_note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize ABRicate TSV results if optional ABRicate screening completed."
    )
    parser.add_argument(
        "--abricate-tsv",
        required=True,
        help="Raw ABRicate TSV output. It may be missing if ABRicate was skipped.",
    )
    parser.add_argument(
        "--status-file",
        required=True,
        help="TSV file describing whether ABRicate completed, failed, or was skipped.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path for summarized ABRicate virulence/database evidence.",
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
        .replace("%", "percent")
    )


def first_available(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value:
            return value.strip()
    return ""


def sample_from_file(value: str) -> str:
    if not value:
        return ""
    path = Path(value)
    for suffix in [".fasta", ".fa", ".fna", ".fas"]:
        if path.name.lower().endswith(suffix):
            return path.name[: -len(suffix)]
    return path.stem or value


def read_status(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "unknown", "ABRicate status file was not found."

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return row.get("status", "unknown"), row.get("note", "")

    return "unknown", "ABRicate status file was empty."


def note_row(status: str, note: str) -> dict[str, str]:
    return {
        "sample_id": "",
        "gene": "",
        "database": "",
        "accession": "",
        "percent_identity": "",
        "percent_coverage": "",
        "sequence": "",
        "start": "",
        "end": "",
        "product": "",
        "evidence_source": f"abricate_{status}",
        "evidence_note": note,
    }


def summarize_abricate(path: Path, status_file: Path) -> list[dict[str, str]]:
    status, status_note = read_status(status_file)
    if status != "completed":
        return [
            note_row(
                status,
                status_note
                or "ABRicate was not completed; no ABRicate evidence was used.",
            )
        ]

    if not path.exists() or path.stat().st_size == 0:
        return [note_row(status, "ABRicate completed but output was missing or empty.")]

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [note_row(status, "ABRicate output had no header row.")]

        reader.fieldnames = [normalize_header(name) for name in reader.fieldnames]
        rows: list[dict[str, str]] = []

        for raw_row in reader:
            row = {
                normalize_header(key): (value or "").strip()
                for key, value in raw_row.items()
                if key is not None
            }
            gene = first_available(row, ["gene", "resistance_gene", "element"])
            database = first_available(row, ["database", "db"])
            product = first_available(
                row,
                ["product", "resistance", "resistance_gene", "annotation"],
            )
            sample_id = sample_from_file(
                first_available(row, ["file", "filename", "sample"])
            )
            rows.append(
                {
                    "sample_id": sample_id,
                    "gene": gene,
                    "database": database,
                    "accession": first_available(row, ["accession", "acc"]),
                    "percent_identity": first_available(
                        row, ["percentidentity", "identity", "percent_id", "pident"]
                    ),
                    "percent_coverage": first_available(
                        row,
                        [
                            "percentcoverage",
                            "coverage",
                            "percent_cov",
                            "pcov",
                            "percent_coverage",
                        ],
                    ),
                    "sequence": first_available(row, ["sequence", "seq", "contig"]),
                    "start": first_available(row, ["start"]),
                    "end": first_available(row, ["end", "stop"]),
                    "product": product,
                    "evidence_source": "optional_abricate_database_hit",
                    "evidence_note": (
                        "Optional ABRicate database hit; useful supporting evidence "
                        "for candidate prioritization, not proof of invasion."
                    ),
                }
            )

    if not rows:
        return [note_row(status, "ABRicate completed but contained no hit rows.")]
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
        rows = summarize_abricate(Path(args.abricate_tsv), Path(args.status_file))
        write_summary(rows, Path(args.output))
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"ABRicate summary rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
