#!/usr/bin/env python3

"""Screen extracted Bakta features for candidate invasome-related keywords."""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


KEYWORD_CATEGORIES: dict[str, list[str]] = {
    "adhesion": [
        "adhesin",
        "adhesion",
        "fibronectin",
        "fibrinogen",
        "collagen",
        "elastin",
        "clf",
        "clfa",
        "clfb",
        "fnb",
        "fnba",
        "fnbb",
        "sdr",
        "sdrd",
        "sdre",
        "ebp",
        "map",
        "spa",
        "surface protein",
    ],
    "invasion_host_interaction": [
        "invasion",
        "internalization",
        "host",
        "host cell",
        "eap",
        "extracellular adherence protein",
        "von willebrand",
        "coagulase",
        "clumping factor",
    ],
    "toxin_hemolysis": [
        "toxin",
        "hemolysin",
        "haemolysin",
        "leukocidin",
        "panton-valentine",
        "pvl",
        "alpha-hemolysin",
        "beta-hemolysin",
        "gamma-hemolysin",
        "hly",
        "hla",
        "hlb",
        "hlg",
        "luk",
        "exfoliative toxin",
        "enterotoxin",
    ],
    "immune_evasion": [
        "capsule",
        "cap",
        "complement",
        "immune",
        "immunoglobulin",
        "protein a",
        "staphylokinase",
        "sak",
        "scin",
        "chp",
        "eap",
        "chemotaxis inhibitory",
    ],
    "secretion_surface": [
        "secretion",
        "secreted",
        "signal peptide",
        "sortase",
        "srt",
        "lpxtg",
        "cell wall",
        "membrane",
        "surface",
        "exported",
        "transporter",
        "lipoprotein",
    ],
    "iron_heme_acquisition": [
        "iron",
        "heme",
        "haem",
        "siderophore",
        "isd",
        "fhu",
        "hpu",
        "ferric",
        "ferrous",
        "transferrin",
        "hemoglobin",
        "haemoglobin",
    ],
    "stress_survival": [
        "stress",
        "oxidative",
        "peroxide",
        "catalase",
        "superoxide",
        "dismutase",
        "heat shock",
        "cold shock",
        "chaperone",
        "alkaline shock",
        "acid shock",
    ],
    "biofilm": [
        "biofilm",
        "ica",
        "polysaccharide intercellular adhesin",
        "slime",
        "matrix",
        "autolysin",
        "atl",
        "bap",
        "accumulation-associated protein",
    ],
    "enzyme_tissue_damage": [
        "protease",
        "lipase",
        "hyaluronidase",
        "nuclease",
        "thermonuclease",
        "staphylocoagulase",
        "phospholipase",
        "tissue",
        "degradation",
        "metalloprotease",
    ],
    "antimicrobial_resistance_context": [
        "resistance",
        "beta-lactamase",
        "betalactamase",
        "mec",
        "meca",
        "mecc",
        "pbp2a",
        "penicillin-binding protein 2a",
        "efflux",
        "multidrug",
        "vancomycin",
        "methicillin",
    ],
}


HIT_COLUMNS = [
    "sample_id",
    "feature_id",
    "gene",
    "product",
    "matched_categories",
    "matched_keywords",
    "evidence_source",
    "raw_annotation_text",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find candidate invasome-related keyword hits in extracted Bakta features."
    )
    parser.add_argument("--features", required=True, help="Input extracted feature TSV.")
    parser.add_argument(
        "--hits-output", required=True, help="Output TSV for candidate keyword hits."
    )
    parser.add_argument(
        "--counts-output", required=True, help="Output TSV for category hit counts."
    )
    return parser.parse_args()


def load_features(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Feature table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Feature table has no header row: {path}")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def keyword_pattern(keyword: str) -> re.Pattern[str]:
    escaped = re.escape(keyword.lower())
    if re.match(r"^[a-z0-9_-]+$", keyword.lower()):
        return re.compile(rf"(?<![a-z0-9]){escaped}(?![a-z0-9])")
    return re.compile(escaped)


COMPILED_KEYWORDS = {
    category: [(keyword, keyword_pattern(keyword)) for keyword in keywords]
    for category, keywords in KEYWORD_CATEGORIES.items()
}


def searchable_text(feature: dict[str, str]) -> str:
    fields = [
        feature.get("gene", ""),
        feature.get("product", ""),
        feature.get("db_xrefs", ""),
        feature.get("raw_annotation_text", ""),
    ]
    return " ".join(fields).lower()


def find_matches(text: str) -> tuple[list[str], list[str]]:
    categories: list[str] = []
    keywords: list[str] = []

    for category, keyword_patterns in COMPILED_KEYWORDS.items():
        category_keywords = [
            keyword for keyword, pattern in keyword_patterns if pattern.search(text)
        ]
        if category_keywords:
            categories.append(category)
            keywords.extend(category_keywords)

    return categories, sorted(set(keywords))


def write_hits(features: list[dict[str, str]], output_path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_COLUMNS, delimiter="\t")
        writer.writeheader()

        for feature in features:
            categories, keywords = find_matches(searchable_text(feature))
            if not categories:
                continue

            for category in categories:
                counts[category] += 1

            writer.writerow(
                {
                    "sample_id": feature.get("sample_id", ""),
                    "feature_id": feature.get("feature_id", ""),
                    "gene": feature.get("gene", ""),
                    "product": feature.get("product", ""),
                    "matched_categories": ";".join(categories),
                    "matched_keywords": ";".join(keywords),
                    "evidence_source": "candidate_keyword_match_from_bakta_annotation",
                    "raw_annotation_text": feature.get("raw_annotation_text", ""),
                }
            )

    return counts


def write_counts(counts: Counter[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "candidate_hit_count"])
        for category in KEYWORD_CATEGORIES:
            writer.writerow([category, counts.get(category, 0)])


def main() -> int:
    args = parse_args()
    try:
        features = load_features(Path(args.features))
        counts = write_hits(features, Path(args.hits_output))
        write_counts(counts, Path(args.counts_output))
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Screened {len(features)} features and found {sum(counts.values())} candidate category matches."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
