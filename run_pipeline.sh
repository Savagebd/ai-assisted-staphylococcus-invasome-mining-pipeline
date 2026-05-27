#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# AI-Assisted Staphylococcus Invasome Mining Pipeline
# Milestone 2+: Bakta annotation + keyword screening + optional AMRFinderPlus
#               + optional ABRicate + explainable candidate scoring + report
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.env"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

info() {
    echo "[INFO] $*"
}

warn() {
    echo "[WARN] $*" >&2
}

require_command() {
    local command_name="$1"
    command -v "${command_name}" >/dev/null 2>&1 || die "Required command not found: ${command_name}"
}

realpath_safe() {
    local path_value="$1"
    python -c 'import os, sys; print(os.path.realpath(os.path.expanduser(sys.argv[1])))' "${path_value}"
}

validate_config_value() {
    local name="$1"
    local value="${!name:-}"
    [[ -n "${value}" ]] || die "${name} is missing or empty in config.env"
}

validate_positive_integer() {
    local value="$1"
    [[ "${value}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer. Found: ${value}"
}

validate_sample_id() {
    local value="$1"
    [[ "${value}" =~ ^[A-Za-z0-9_-]+$ ]] || die "SAMPLE_ID may only contain letters, numbers, underscores, and hyphens. Found: ${value}"
}

is_true() {
    local value="${1:-false}"
    case "${value,,}" in
        true|yes|y|1) return 0 ;;
        *) return 1 ;;
    esac
}

prepare_output_directories() {
    mkdir -p \
        "${PROJECT_DIR}/01_Genome_FASTA" \
        "${PROJECT_DIR}/02_Bakta_Annotation" \
        "${PROJECT_DIR}/03_Invasome_Feature_Extraction" \
        "${PROJECT_DIR}/04_Invasome_Keyword_Screening" \
        "${PROJECT_DIR}/05_AMRFinderPlus" \
        "${PROJECT_DIR}/06_ABRicate" \
        "${PROJECT_DIR}/07_Protein_Localization" \
        "${PROJECT_DIR}/08_Invasome_Summary" \
        "${PROJECT_DIR}/09_Final_Report" \
        "${PROJECT_DIR}/logs"
}

clean_generated_outputs_for_rerun() {
    info "Cleaning generated outputs from previous runs. Input FASTA files are not touched."
    rm -rf \
        "${PROJECT_DIR}/02_Bakta_Annotation" \
        "${PROJECT_DIR}/03_Invasome_Feature_Extraction" \
        "${PROJECT_DIR}/04_Invasome_Keyword_Screening" \
        "${PROJECT_DIR}/05_AMRFinderPlus" \
        "${PROJECT_DIR}/06_ABRicate" \
        "${PROJECT_DIR}/07_Protein_Localization" \
        "${PROJECT_DIR}/08_Invasome_Summary" \
        "${PROJECT_DIR}/09_Final_Report" \
        "${PROJECT_DIR}/logs"
}

write_status_tsv() {
    local output_file="$1"
    local status="$2"
    local note="$3"
    {
        printf "status\tnote\n"
        printf "%s\t%s\n" "${status}" "${note}"
    } > "${output_file}"
}

[[ -f "${CONFIG_FILE}" ]] || die "Missing config.env. Create one with: cp config.example.env config.env, then edit config.env."

info "Loading configuration from config.env."

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

RUN_ABRICATE="${RUN_ABRICATE:-false}"
ABRICATE_DB="${ABRICATE_DB:-}"
AMRFINDER_DB="${AMRFINDER_DB:-}"

validate_config_value "PROJECT_DIR"
validate_config_value "SAMPLE_ID"
validate_config_value "GENOME_FASTA"
validate_config_value "THREADS"
validate_config_value "BAKTA_DB"

validate_sample_id "${SAMPLE_ID}"
validate_positive_integer "${THREADS}"

require_command "python"
require_command "bakta"

PROJECT_DIR="$(realpath_safe "${PROJECT_DIR}")"
GENOME_FASTA="$(realpath_safe "${GENOME_FASTA}")"
BAKTA_DB="$(realpath_safe "${BAKTA_DB}")"
GENOME_INPUT_DIR="$(realpath_safe "${PROJECT_DIR}/01_Genome_FASTA")"
if [[ -n "${AMRFINDER_DB}" ]]; then
    AMRFINDER_DB="$(realpath_safe "${AMRFINDER_DB}")"
fi

[[ -d "${PROJECT_DIR}" ]] || die "PROJECT_DIR does not exist: ${PROJECT_DIR}"
[[ "${PROJECT_DIR}" != "/" ]] || die "PROJECT_DIR cannot be the filesystem root."
[[ -f "${GENOME_FASTA}" ]] || die "GENOME_FASTA does not exist: ${GENOME_FASTA}"
[[ -d "${BAKTA_DB}" ]] || die "BAKTA_DB does not exist: ${BAKTA_DB}"

case "${GENOME_FASTA}" in
    "${GENOME_INPUT_DIR}"/*) ;;
    *) die "GENOME_FASTA must be inside PROJECT_DIR/01_Genome_FASTA after resolving real paths. Found: ${GENOME_FASTA}" ;;
esac

case "${GENOME_FASTA}" in
    *.fa|*.fna|*.fasta|*.fas|*.FA|*.FNA|*.FASTA|*.FAS) ;;
    *) die "GENOME_FASTA must be an assembled genome FASTA file with a FASTA-like extension, not raw FASTQ." ;;
esac

clean_generated_outputs_for_rerun
prepare_output_directories

BAKTA_OUT_DIR="${PROJECT_DIR}/02_Bakta_Annotation"
FEATURE_TABLE="${PROJECT_DIR}/03_Invasome_Feature_Extraction/bakta_features_extracted.tsv"
KEYWORD_HITS="${PROJECT_DIR}/04_Invasome_Keyword_Screening/invasome_keyword_hits.tsv"
CATEGORY_COUNTS="${PROJECT_DIR}/08_Invasome_Summary/invasome_category_counts.tsv"
AMRFINDER_OUT="${PROJECT_DIR}/05_AMRFinderPlus/amrfinder_${SAMPLE_ID}.tsv"
AMRFINDER_SUMMARY="${PROJECT_DIR}/08_Invasome_Summary/amrfinder_summary.tsv"
AMRFINDER_STATUS="${PROJECT_DIR}/08_Invasome_Summary/amrfinder_status.tsv"
ABRICATE_OUT="${PROJECT_DIR}/06_ABRicate/abricate_${SAMPLE_ID}.tsv"
ABRICATE_SUMMARY="${PROJECT_DIR}/08_Invasome_Summary/abricate_virulence_summary.tsv"
ABRICATE_STATUS="${PROJECT_DIR}/08_Invasome_Summary/abricate_status.tsv"
LOCALIZATION_CANDIDATES="${PROJECT_DIR}/07_Protein_Localization/protein_localization_candidates.tsv"
CANDIDATE_SCORES="${PROJECT_DIR}/08_Invasome_Summary/invasome_candidate_scores.tsv"
FINAL_REPORT="${PROJECT_DIR}/09_Final_Report/invasome_mining_final_report.txt"
BAKTA_LOG="${PROJECT_DIR}/logs/bakta.log"
AMRFINDER_LOG="${PROJECT_DIR}/logs/amrfinderplus.log"
ABRICATE_LOG="${PROJECT_DIR}/logs/abricate.log"

info "Running Bakta annotation for sample: ${SAMPLE_ID}"
bakta \
    --db "${BAKTA_DB}" \
    --threads "${THREADS}" \
    --output "${BAKTA_OUT_DIR}" \
    --prefix "${SAMPLE_ID}" \
    --force \
    "${GENOME_FASTA}" \
    > "${BAKTA_LOG}" 2>&1

BAKTA_TSV="${BAKTA_OUT_DIR}/${SAMPLE_ID}.tsv"
BAKTA_PROTEINS="${BAKTA_OUT_DIR}/${SAMPLE_ID}.faa"
[[ -f "${BAKTA_TSV}" ]] || die "Bakta completed but expected TSV was not found: ${BAKTA_TSV}"

info "Extracting Bakta features."
python "${SCRIPT_DIR}/scripts/01_extract_bakta_features.py" \
    --sample-id "${SAMPLE_ID}" \
    --bakta-tsv "${BAKTA_TSV}" \
    --output "${FEATURE_TABLE}"

info "Screening annotation text for candidate invasome-related keyword evidence."
python "${SCRIPT_DIR}/scripts/02_screen_invasome_keywords.py" \
    --features "${FEATURE_TABLE}" \
    --hits-output "${KEYWORD_HITS}" \
    --counts-output "${CATEGORY_COUNTS}"

if command -v amrfinder >/dev/null 2>&1; then
    if [[ -f "${BAKTA_PROTEINS}" ]]; then
        AMRFINDER_DATABASE_ARGS=()
        AMRFINDER_DB_NOTE="default AMRFinderPlus database location"
        if [[ -n "${AMRFINDER_DB}" ]]; then
            if [[ -d "${AMRFINDER_DB}" ]]; then
                AMRFINDER_DATABASE_ARGS=(--database "${AMRFINDER_DB}")
                AMRFINDER_DB_NOTE="${AMRFINDER_DB}"
            else
                warn "Configured AMRFINDER_DB does not exist. Skipping AMRFinderPlus."
                write_status_tsv "${AMRFINDER_STATUS}" "skipped" "Configured AMRFINDER_DB does not exist: ${AMRFINDER_DB}"
            fi
        fi

        if [[ ! -f "${AMRFINDER_STATUS}" ]]; then
            info "Running optional AMRFinderPlus screening on Bakta protein output."
            if amrfinder \
                -p "${BAKTA_PROTEINS}" \
                -o "${AMRFINDER_OUT}" \
                --threads "${THREADS}" \
                "${AMRFINDER_DATABASE_ARGS[@]}" \
                > "${AMRFINDER_LOG}" 2>&1; then
                write_status_tsv "${AMRFINDER_STATUS}" "completed" "AMRFinderPlus completed successfully using ${AMRFINDER_DB_NOTE}."
            else
                warn "AMRFinderPlus was found but did not complete successfully. Continuing without AMRFinderPlus evidence."
                write_status_tsv "${AMRFINDER_STATUS}" "failed" "AMRFinderPlus was installed but failed. See logs/amrfinderplus.log."
            fi
        fi
    else
        warn "Bakta protein FASTA was not found. Skipping AMRFinderPlus."
        write_status_tsv "${AMRFINDER_STATUS}" "skipped" "Bakta protein FASTA was not found: ${BAKTA_PROTEINS}"
    fi
else
    warn "AMRFinderPlus command 'amrfinder' was not found. Skipping optional AMRFinderPlus screening."
    write_status_tsv "${AMRFINDER_STATUS}" "skipped" "AMRFinderPlus command 'amrfinder' was not available in the active environment."
fi

info "Summarizing optional AMRFinderPlus evidence."
python "${SCRIPT_DIR}/scripts/03_summarize_amrfinder_results.py" \
    --amrfinder-tsv "${AMRFINDER_OUT}" \
    --status-file "${AMRFINDER_STATUS}" \
    --output "${AMRFINDER_SUMMARY}"

if is_true "${RUN_ABRICATE}"; then
    if [[ -z "${ABRICATE_DB}" ]]; then
        warn "RUN_ABRICATE is true but ABRICATE_DB is empty. Skipping ABRicate."
        write_status_tsv "${ABRICATE_STATUS}" "skipped" "RUN_ABRICATE=true but ABRICATE_DB was empty."
    elif command -v abricate >/dev/null 2>&1; then
        info "Running optional ABRicate screening with database: ${ABRICATE_DB}"
        if abricate \
            --db "${ABRICATE_DB}" \
            "${GENOME_FASTA}" \
            > "${ABRICATE_OUT}" 2> "${ABRICATE_LOG}"; then
            write_status_tsv "${ABRICATE_STATUS}" "completed" "ABRicate completed successfully with database: ${ABRICATE_DB}"
        else
            warn "ABRicate was found but did not complete successfully. Continuing without ABRicate evidence."
            write_status_tsv "${ABRICATE_STATUS}" "failed" "ABRicate failed or the requested database was unavailable: ${ABRICATE_DB}. See logs/abricate.log."
        fi
    else
        warn "ABRicate command 'abricate' was not found. Skipping optional ABRicate screening."
        write_status_tsv "${ABRICATE_STATUS}" "skipped" "ABRicate command 'abricate' was not available in the active environment."
    fi
else
    info "RUN_ABRICATE is not true. Skipping optional ABRicate screening."
    write_status_tsv "${ABRICATE_STATUS}" "skipped" "RUN_ABRICATE was not set to true."
fi

info "Summarizing optional ABRicate virulence evidence."
python "${SCRIPT_DIR}/scripts/06_summarize_abricate_results.py" \
    --abricate-tsv "${ABRICATE_OUT}" \
    --status-file "${ABRICATE_STATUS}" \
    --output "${ABRICATE_SUMMARY}"

info "Predicting heuristic protein localization and surface-associated candidates."
python "${SCRIPT_DIR}/scripts/07_predict_surface_localization_candidates.py" \
    --features "${FEATURE_TABLE}" \
    --protein-fasta "${BAKTA_PROTEINS}" \
    --output "${LOCALIZATION_CANDIDATES}"

info "Scoring candidate invasome-related genes."
python "${SCRIPT_DIR}/scripts/04_score_invasome_candidates.py" \
    --keyword-hits "${KEYWORD_HITS}" \
    --amrfinder-summary "${AMRFINDER_SUMMARY}" \
    --abricate-summary "${ABRICATE_SUMMARY}" \
    --localization-summary "${LOCALIZATION_CANDIDATES}" \
    --output "${CANDIDATE_SCORES}"

info "Generating final report."
python "${SCRIPT_DIR}/scripts/05_generate_final_invasome_report.py" \
    --sample-id "${SAMPLE_ID}" \
    --genome-fasta "${GENOME_FASTA}" \
    --features "${FEATURE_TABLE}" \
    --keyword-hits "${KEYWORD_HITS}" \
    --category-counts "${CATEGORY_COUNTS}" \
    --amrfinder-status "${AMRFINDER_STATUS}" \
    --amrfinder-summary "${AMRFINDER_SUMMARY}" \
    --abricate-status "${ABRICATE_STATUS}" \
    --abricate-summary "${ABRICATE_SUMMARY}" \
    --localization-summary "${LOCALIZATION_CANDIDATES}" \
    --candidate-scores "${CANDIDATE_SCORES}" \
    --output "${FINAL_REPORT}"

echo
echo "Pipeline completed successfully."
echo "Key outputs:"
echo "  Bakta annotation:        ${BAKTA_OUT_DIR}"
echo "  Extracted features:      ${FEATURE_TABLE}"
echo "  Keyword hits:            ${KEYWORD_HITS}"
echo "  Category counts:         ${CATEGORY_COUNTS}"
echo "  AMRFinderPlus summary:   ${AMRFINDER_SUMMARY}"
echo "  ABRicate summary:        ${ABRICATE_SUMMARY}"
echo "  Localization candidates: ${LOCALIZATION_CANDIDATES}"
echo "  Candidate scores:        ${CANDIDATE_SCORES}"
echo "  Final report:            ${FINAL_REPORT}"
echo "  Logs:                    ${PROJECT_DIR}/logs"
