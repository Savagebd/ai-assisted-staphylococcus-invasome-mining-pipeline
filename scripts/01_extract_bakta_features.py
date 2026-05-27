#!/usr/bin/env python3

"""Extract a robust, clean feature table from a Bakta TSV annotation file."""

import argparse
import csv
import io
import sys
from pathlib import Path


OUTPUT_COLUMNS = [
    "sample_id",
    "feature_id",
    "contig",
    "start",
    "stop",
    "strand",
    "feature_type",
    "gene",
    "product",
    "db_xrefs",
    "nt_sequence_length",
    "aa_sequence_length",
    "raw_annotation_text",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract selected fields from a Bakta TSV file."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument(
        "--bakta-tsv",
        required=True,
        help="Input Bakta TSV file, usually 02_Bakta_Annotation/<sample>.tsv.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path for extracted feature records.",
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


def sequence_length(start: str, stop: str) -> str:
    try:
        left = int(start)
        right = int(stop)
    except ValueError:
        return ""

    return str(abs(right - left) + 1)


def aa_length(nt_length: str, feature_type: str) -> str:
    if not nt_length:
        return ""
    if feature_type.lower() not in {"cds", "coding_sequence", "protein"}:
        return ""
    try:
        return str(max((int(nt_length) // 3) - 1, 0))
    except ValueError:
        return ""


def build_feature_id(row: dict[str, str], fallback_index: int) -> str:
    explicit_id = first_available(
        row,
        [
            "locus_tag",
            "locus",
            "id",
            "protein_id",
            "gene",
            "db_xrefs",
            "db_xref",
        ],
    )
    if explicit_id:
        return explicit_id
    return f"feature_{fallback_index:06d}"


def raw_annotation_text(row: dict[str, str]) -> str:
    useful_values = [value.strip() for value in row.values() if value and value.strip()]
    return " | ".join(useful_values)


def read_bakta_records(path: Path, sample_id: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Bakta TSV file not found: {path}")
    if not path.is_file():
        raise ValueError(f"Bakta TSV path is not a file: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        data_lines: list[str] = []
        header_found = False
        for line in handle:
            if not header_found:
                if line.startswith("#") and "\t" in line and "Sequence Id" in line:
                    data_lines.append(line.lstrip("#"))
                    header_found = True
                elif not line.startswith("#") and "\t" in line:
                    data_lines.append(line)
                    header_found = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            data_lines.append(line)

        if not data_lines:
            raise ValueError(f"Bakta TSV file has no parseable table: {path}")

        reader = csv.DictReader(io.StringIO("".join(data_lines)), delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Bakta TSV file has no header row: {path}")

        reader.fieldnames = [normalize_header(name) for name in reader.fieldnames]
        records: list[dict[str, str]] = []

        for index, row in enumerate(reader, start=1):
            normalized = {
                normalize_header(key): (value or "").strip()
                for key, value in row.items()
                if key is not None
            }
            start = first_available(normalized, ["start", "begin"])
            stop = first_available(normalized, ["stop", "end"])
            feature_type = first_available(
                normalized, ["type", "feature_type", "feature"]
            )
            nt_len = first_available(
                normalized,
                [
                    "nt_sequence_length",
                    "nucleotide_sequence_length",
                    "nucleotide_length",
                    "length",
                ],
            )
            if not nt_len:
                nt_len = sequence_length(start, stop)

            aa_len = first_available(
                normalized,
                ["aa_sequence_length", "protein_sequence_length", "protein_length"],
            )
            if not aa_len:
                aa_len = aa_length(nt_len, feature_type)

            records.append(
                {
                    "sample_id": sample_id,
                    "feature_id": build_feature_id(normalized, index),
                    "contig": first_available(
                        normalized,
                        ["contig", "sequence_id", "sequence", "seq_id", "replicon"],
                    ),
                    "start": start,
                    "stop": stop,
                    "strand": first_available(normalized, ["strand"]),
                    "feature_type": feature_type,
                    "gene": first_available(normalized, ["gene", "gene_name"]),
                    "product": first_available(
                        normalized, ["product", "function", "description"]
                    ),
                    "db_xrefs": first_available(
                        normalized,
                        ["db_xrefs", "dbxref", "db_xref", "database", "inference"],
                    ),
                    "nt_sequence_length": nt_len,
                    "aa_sequence_length": aa_len,
                    "raw_annotation_text": raw_annotation_text(normalized),
                }
            )

    return records


def write_records(records: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    args = parse_args()
    try:
        records = read_bakta_records(Path(args.bakta_tsv), args.sample_id)
        write_records(records, Path(args.output))
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Extracted {len(records)} Bakta feature records.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
