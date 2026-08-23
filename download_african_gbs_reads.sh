#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# GBS-Sentinel — African GBS RAW-READ cohort (FASTQ)
#
# Nigeria      : 21  (PRJNA632041, BioSamples SAMN18352391–SAMN18352411)
# South Africa : 19  (PRJNA479604, human S. agalactiae paired WGS)
# Total        : 40
#
# Uses SRA Toolkit (prefetch + fasterq-dump), NOT datasets for SRA.
# Does NOT modify raw_data/, curated_data/, or previous results.
# Output: african_raw_reads/
# ============================================================

PROJECT_DIR="${HOME}/projects/GBS-Sentinel"
OUT_DIR="${PROJECT_DIR}/african_raw_reads"
WORK_DIR="${PROJECT_DIR}/african_raw_work"

NIGERIA_PROJECT="PRJNA632041"
SOUTH_AFRICA_PROJECT="PRJNA479604"
NIGERIA_TARGET=10
SOUTH_AFRICA_TARGET=10
TOTAL_TARGET=20

# Published Port Harcourt maternal GBS BioSamples (Bob-Manuel et al.)
NIGERIA_BIOSAMPLES_START=18352391
NIGERIA_BIOSAMPLES_END=18352411

METADATA="${OUT_DIR}/african_gbs_reads_metadata.tsv"
ACCESSIONS="${OUT_DIR}/african_gbs_run_accessions.txt"
MANIFEST="${OUT_DIR}/african_gbs_fastq_manifest.tsv"

NIGERIA_RUNINFO="${WORK_DIR}/nigeria_runinfo.csv"
SOUTH_AFRICA_RUNINFO="${WORK_DIR}/south_africa_runinfo.csv"
NIGERIA_SELECTED="${WORK_DIR}/nigeria_selected.tsv"
SOUTH_AFRICA_SELECTED="${WORK_DIR}/south_africa_selected.tsv"

mkdir -p "${OUT_DIR}" "${WORK_DIR}"

echo
echo "============================================================"
echo "GBS-Sentinel African GBS RAW-READ Cohort"
echo "============================================================"
echo
echo "Project : ${PROJECT_DIR}"
echo "FASTQ   : ${OUT_DIR}"
echo
echo "Target:"
echo "  Nigeria      : ${NIGERIA_TARGET}  (${NIGERIA_PROJECT})"
echo "  South Africa : ${SOUTH_AFRICA_TARGET}  (${SOUTH_AFRICA_PROJECT})"
echo "  Total        : ${TOTAL_TARGET}"
echo

# ============================================================
# 1. SRA Toolkit
# ============================================================
echo "------------------------------------------------------------"
echo "STEP 1 — Checking SRA Toolkit"
echo "------------------------------------------------------------"

if ! command -v prefetch >/dev/null 2>&1 || ! command -v fasterq-dump >/dev/null 2>&1; then
    echo
    echo "ERROR: SRA Toolkit not found (need prefetch and fasterq-dump)."
    echo "Install, for example:"
    echo "  conda install -c bioconda sra-tools -y"
    exit 1
fi

echo "prefetch:     $(command -v prefetch)"
echo "fasterq-dump: $(command -v fasterq-dump)"
prefetch --version 2>&1 | head -n 3 || true

# ============================================================
# 2. Fetch RunInfo (SRA Run Selector CSV via ENA/NCBI)
# ============================================================
echo
echo "------------------------------------------------------------"
echo "STEP 2 — Downloading RunInfo tables"
echo "------------------------------------------------------------"

# NCBI SRA RunInfo CSV for a BioProject
download_runinfo() {
    local bioproject="$1"
    local outfile="$2"
    local url="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-info/runinfo?acc=${bioproject}"

    echo "Fetching RunInfo for ${bioproject} ..."
    if ! curl -fsSL --retry 3 --retry-delay 2 "${url}" -o "${outfile}"; then
        # Fallback: ENA filereport
        local ena="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${bioproject}&result=read_run&fields=run_accession,sample_accession,scientific_name,library_layout,fastq_ftp,base_count,read_count,country,host,instrument_model&format=tsv&limit=0"
        echo "NCBI RunInfo failed; trying ENA filereport..."
        curl -fsSL --retry 3 --retry-delay 2 "${ena}" -o "${outfile}"
    fi

    if [ ! -s "${outfile}" ]; then
        echo "ERROR: Empty RunInfo for ${bioproject}"
        exit 1
    fi
    echo "  Saved: ${outfile} ($(wc -l < "${outfile}" | tr -d ' ') lines)"
}

download_runinfo "${NIGERIA_PROJECT}" "${NIGERIA_RUNINFO}"
download_runinfo "${SOUTH_AFRICA_PROJECT}" "${SOUTH_AFRICA_RUNINFO}"

# ============================================================
# 3. Select Nigeria: SAMN18352391–SAMN18352411 when present
# ============================================================
echo
echo "------------------------------------------------------------"
echo "STEP 3 — Selecting Nigerian runs (Port Harcourt BioSamples)"
echo "------------------------------------------------------------"

python3 - \
    "${NIGERIA_RUNINFO}" \
    "${NIGERIA_SELECTED}" \
    "${NIGERIA_TARGET}" \
    "${NIGERIA_BIOSAMPLES_START}" \
    "${NIGERIA_BIOSAMPLES_END}" \
    <<'PY'
import csv
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
out = Path(sys.argv[2])
target = int(sys.argv[3])
bs_start = int(sys.argv[4])
bs_end = int(sys.argv[5])

wanted = {f"SAMN{i}" for i in range(bs_start, bs_end + 1)}

def norm_header(h):
    return re.sub(r"[^a-z0-9]+", "_", (h or "").strip().lower()).strip("_")

text = path.read_text(encoding="utf-8", errors="replace")
# Detect delimiter
dialect = csv.Sniffer().sniff(text.splitlines()[0] if text else "a,b", delimiters=",\t")
reader = csv.DictReader(text.splitlines(), delimiter=dialect.delimiter)
if not reader.fieldnames:
    raise SystemExit("ERROR: Nigeria RunInfo has no header")

# Map normalized -> original
keymap = {norm_header(k): k for k in reader.fieldnames}

def get(row, *cands):
    for c in cands:
        k = keymap.get(norm_header(c))
        if k and row.get(k):
            return str(row[k]).strip()
    # fuzzy
    for nk, ok in keymap.items():
        for c in cands:
            if norm_header(c) in nk or nk in norm_header(c):
                if row.get(ok):
                    return str(row[ok]).strip()
    return ""

rows = []
for row in reader:
    run = get(row, "Run", "run", "run_accession", "Run Accession")
    biosample = get(row, "BioSample", "biosample", "sample_accession", "Sample Name", "Sample")
    org = get(row, "ScientificName", "scientific_name", "Organism", "organism")
    layout = get(row, "LibraryLayout", "library_layout", "Library Layout")
    spots = get(row, "spots", "read_count", "Read Count")
    bases = get(row, "bases", "base_count", "Base Count")
    size = get(row, "size_MB", "size", "Bytes", "size_mb")
    country = get(row, "geo_loc_name_country_calc", "Country", "country", "geo_loc_name")
    host = get(row, "Host", "host")

    if not run:
        continue
    # Prefer published BioSamples
    if biosample and biosample not in wanted:
        # still allow if org is GBS and we'll fill later
        pass

    org_l = org.lower()
    if org and "agalactiae" not in org_l and "group b" not in org_l:
        continue

    layout_l = layout.upper()
    # Prefer PAIRED
    rows.append({
        "run": run,
        "biosample": biosample,
        "organism": org or "Streptococcus agalactiae",
        "layout": layout or "",
        "spots": spots,
        "bases": bases,
        "size": size,
        "country": country or "Nigeria",
        "host": host,
        "in_published_set": "yes" if biosample in wanted else "no",
    })

# Prefer published BioSamples, paired, unique biosample
published = [r for r in rows if r["in_published_set"] == "yes"]
others = [r for r in rows if r["in_published_set"] != "yes"]

def sort_key(r):
    paired = 0 if "PAIR" in r["layout"].upper() else 1
    return (paired, r["biosample"] or "", r["run"])

published.sort(key=sort_key)
others.sort(key=sort_key)

selected = []
seen_bs = set()
seen_run = set()

for pool in (published, others):
    for r in pool:
        if len(selected) >= target:
            break
        if r["run"] in seen_run:
            continue
        bs = r["biosample"] or r["run"]
        if bs in seen_bs:
            continue
        seen_bs.add(bs)
        seen_run.add(r["run"])
        selected.append(r)
    if len(selected) >= target:
        break

selected = selected[:target]

with out.open("w", encoding="utf-8", newline="") as fh:
    cols = ["run", "biosample", "organism", "layout", "spots", "bases", "size", "country", "host", "in_published_set"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in selected:
        w.writerow(r)

print(f"Nigeria candidate rows in RunInfo (GBS-like): {len(rows)}")
print(f"Published BioSample hits: {sum(1 for r in rows if r['in_published_set']=='yes')}")
print(f"Selected Nigeria runs: {len(selected)}")
for i, r in enumerate(selected, 1):
    print(f"  {i:2d}. {r['run']}  {r['biosample']}  {r['layout']}  published={r['in_published_set']}")
if len(selected) < target:
    print(f"WARNING: only {len(selected)} Nigeria runs selected (wanted {target})")
PY

NGA_N=$(awk 'NR>1{n++} END{print n+0}' "${NIGERIA_SELECTED}")
if [ "${NGA_N}" -lt 1 ]; then
    echo "ERROR: No Nigerian runs selected. Inspect ${NIGERIA_RUNINFO}"
    exit 1
fi
if [ "${NGA_N}" -lt "${NIGERIA_TARGET}" ]; then
    echo "WARNING: Nigeria selected ${NGA_N} < ${NIGERIA_TARGET}. Continuing with available runs."
fi

# ============================================================
# 4. Select South Africa: human GBS paired from PRJNA479604
# ============================================================
echo
echo "------------------------------------------------------------"
echo "STEP 4 — Selecting South African runs (PRJNA479604)"
echo "------------------------------------------------------------"

python3 - \
    "${SOUTH_AFRICA_RUNINFO}" \
    "${SOUTH_AFRICA_SELECTED}" \
    "${SOUTH_AFRICA_TARGET}" \
    <<'PY'
import csv
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
out = Path(sys.argv[2])
target = int(sys.argv[3])

def norm_header(h):
    return re.sub(r"[^a-z0-9]+", "_", (h or "").strip().lower()).strip("_")

text = path.read_text(encoding="utf-8", errors="replace")
dialect = csv.Sniffer().sniff(text.splitlines()[0] if text else "a,b", delimiters=",\t")
reader = csv.DictReader(text.splitlines(), delimiter=dialect.delimiter)
if not reader.fieldnames:
    raise SystemExit("ERROR: South Africa RunInfo has no header")

keymap = {norm_header(k): k for k in reader.fieldnames}

def get(row, *cands):
    for c in cands:
        k = keymap.get(norm_header(c))
        if k and row.get(k):
            return str(row[k]).strip()
    for nk, ok in keymap.items():
        for c in cands:
            if norm_header(c) in nk or nk in norm_header(c):
                if row.get(ok):
                    return str(row[ok]).strip()
    return ""

rows = []
for row in reader:
    run = get(row, "Run", "run", "run_accession")
    biosample = get(row, "BioSample", "biosample", "sample_accession", "Sample")
    org = get(row, "ScientificName", "scientific_name", "Organism")
    layout = get(row, "LibraryLayout", "library_layout")
    spots = get(row, "spots", "read_count")
    bases = get(row, "bases", "base_count")
    size = get(row, "size_MB", "size", "Bytes")
    country = get(row, "geo_loc_name_country_calc", "Country", "country", "geo_loc_name")
    host = get(row, "Host", "host")

    if not run:
        continue
    org_l = (org or "").lower()
    if org and "agalactiae" not in org_l and "group b" not in org_l:
        continue

    host_l = (host or "").lower()
    country_l = (country or "").lower()
    # Prefer human + South Africa when fields exist; do not drop if missing
    # (project is SA GBS study)
    rows.append({
        "run": run,
        "biosample": biosample,
        "organism": org or "Streptococcus agalactiae",
        "layout": layout or "",
        "spots": spots,
        "bases": bases,
        "size": size,
        "country": country or "South Africa",
        "host": host or "",
    })

def score(r):
    paired = 0 if "PAIR" in r["layout"].upper() else 1
    human = 0 if ("human" in r["host"].lower() or "homo" in r["host"].lower() or not r["host"]) else 1
    sa = 0 if ("south africa" in r["country"].lower() or not r["country"] or r["country"] == "South Africa") else 1
    return (paired, human, sa, r["biosample"] or "", r["run"])

rows.sort(key=score)

selected = []
seen_bs = set()
seen_run = set()
for r in rows:
    if len(selected) >= target:
        break
    if r["run"] in seen_run:
        continue
    bs = r["biosample"] or r["run"]
    if bs in seen_bs:
        continue
    # Prefer paired
    if r["layout"] and "PAIR" not in r["layout"].upper():
        continue
    seen_bs.add(bs)
    seen_run.add(r["run"])
    selected.append(r)

# If not enough paired, allow non-paired to fill
if len(selected) < target:
    for r in rows:
        if len(selected) >= target:
            break
        if r["run"] in seen_run:
            continue
        bs = r["biosample"] or r["run"]
        if bs in seen_bs:
            continue
        seen_bs.add(bs)
        seen_run.add(r["run"])
        selected.append(r)

selected = selected[:target]

with out.open("w", encoding="utf-8", newline="") as fh:
    cols = ["run", "biosample", "organism", "layout", "spots", "bases", "size", "country", "host"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in selected:
        w.writerow(r)

print(f"South Africa GBS-like rows: {len(rows)}")
print(f"Selected SA runs: {len(selected)}")
for i, r in enumerate(selected, 1):
    print(f"  {i:2d}. {r['run']}  {r['biosample']}  {r['layout']}  {r['country']}")
if len(selected) < target:
    print(f"WARNING: only {len(selected)} SA runs (wanted {target})")
PY

ZAF_N=$(awk 'NR>1{n++} END{print n+0}' "${SOUTH_AFRICA_SELECTED}")
if [ "${ZAF_N}" -lt 1 ]; then
    echo "ERROR: No South African runs selected. Inspect ${SOUTH_AFRICA_RUNINFO}"
    exit 1
fi

# ============================================================
# 5. Combined manifest + size estimate BEFORE download
# ============================================================
echo
echo "============================================================"
echo "SELECTED COHORT — BEFORE DOWNLOAD"
echo "============================================================"

{
    printf "country\tbioproject\trun\tbiosample\torganism\tlayout\tbases\tsize_field\n"
    awk -F '\t' -v p="${NIGERIA_PROJECT}" 'NR>1{
        printf "Nigeria\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", p,$1,$2,$3,$4,$6,$7
    }' "${NIGERIA_SELECTED}"
    awk -F '\t' -v p="${SOUTH_AFRICA_PROJECT}" 'NR>1{
        printf "South Africa\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", p,$1,$2,$3,$4,$6,$7
    }' "${SOUTH_AFRICA_SELECTED}"
} > "${METADATA}"

{
    awk -F '\t' 'NR>1{print $1}' "${NIGERIA_SELECTED}"
    awk -F '\t' 'NR>1{print $1}' "${SOUTH_AFRICA_SELECTED}"
} > "${ACCESSIONS}"

TOTAL_RUNS=$(wc -l < "${ACCESSIONS}" | tr -d ' ')
echo
column -t -s $'\t' "${METADATA}" || cat "${METADATA}"
echo
echo "Runs selected: ${TOTAL_RUNS}"
echo "Metadata: ${METADATA}"
echo "Run list: ${ACCESSIONS}"
echo
echo "IMPORTANT: Review the table above."
echo "Disk: 40 WGS pairs can be tens of GB compressed / more unzipped."
echo
echo "Press Enter to start download, or Ctrl+C to abort and trim the lists."
read -r _

# ============================================================
# 6. Download with prefetch + fasterq-dump
# ============================================================
echo
echo "------------------------------------------------------------"
echo "STEP 6 — Downloading and converting to FASTQ"
echo "------------------------------------------------------------"

download_one_run() {
    local run="$1"
    local prefix="$2"   # NGA or ZAF
    local sample="${prefix}_${run}"

    local r1="${OUT_DIR}/${sample}_1.fastq.gz"
    local r2="${OUT_DIR}/${sample}_2.fastq.gz"

    if [ -s "${r1}" ] && [ -s "${r2}" ]; then
        echo "  SKIP (exists): ${sample}"
        return 0
    fi

    echo "  prefetch ${run} ..."
    prefetch "${run}" --output-directory "${WORK_DIR}/sra" 2>&1 | tail -n 5

    echo "  fasterq-dump ${run} ..."
    (
        cd "${WORK_DIR}"
        fasterq-dump "${run}" \
            --outdir "${OUT_DIR}" \
            --split-files \
            --threads 4 \
            --skip-technical \
            --force
    )

    # fasterq-dump names: RUN_1.fastq RUN_2.fastq
    if [ -f "${OUT_DIR}/${run}_1.fastq" ]; then
        gzip -f "${OUT_DIR}/${run}_1.fastq"
        mv -f "${OUT_DIR}/${run}_1.fastq.gz" "${r1}"
    fi
    if [ -f "${OUT_DIR}/${run}_2.fastq" ]; then
        gzip -f "${OUT_DIR}/${run}_2.fastq"
        mv -f "${OUT_DIR}/${run}_2.fastq.gz" "${r2}"
    fi

    # Single-end fallback
    if [ ! -s "${r1}" ] && [ -f "${OUT_DIR}/${run}.fastq" ]; then
        gzip -f "${OUT_DIR}/${run}.fastq"
        mv -f "${OUT_DIR}/${run}.fastq.gz" "${OUT_DIR}/${sample}.fastq.gz"
        echo "  NOTE: single-end ${sample}.fastq.gz"
        return 0
    fi

    if [ ! -s "${r1}" ]; then
        echo "ERROR: Missing R1 for ${run}"
        exit 1
    fi
    echo "  OK ${r1}"
    if [ -s "${r2}" ]; then
        echo "  OK ${r2}"
    fi
}

echo "Nigeria..."
while IFS=$'\t' read -r run rest; do
    [ -z "${run}" ] && continue
    download_one_run "${run}" "NGA"
done < <(awk -F '\t' 'NR>1{print $1}' "${NIGERIA_SELECTED}")

echo "South Africa..."
while IFS=$'\t' read -r run rest; do
    [ -z "${run}" ] && continue
    download_one_run "${run}" "ZAF"
done < <(awk -F '\t' 'NR>1{print $1}' "${SOUTH_AFRICA_SELECTED}")

# ============================================================
# 7. Manifest for Nextflow
# ============================================================
echo
echo "------------------------------------------------------------"
echo "STEP 7 — Building FASTQ manifest"
echo "------------------------------------------------------------"

{
    printf "sample_id\tcountry\trun\tr1\tr2\n"
    for f in "${OUT_DIR}"/NGA_*_1.fastq.gz; do
        [ -e "$f" ] || continue
        base=$(basename "$f" _1.fastq.gz)
        run=${base#NGA_}
        r2="${OUT_DIR}/${base}_2.fastq.gz"
        printf "NGA_%s\tNigeria\t%s\t%s\t%s\n" "${run}" "${run}" "$f" "$r2"
    done
    for f in "${OUT_DIR}"/ZAF_*_1.fastq.gz; do
        [ -e "$f" ] || continue
        base=$(basename "$f" _1.fastq.gz)
        run=${base#ZAF_}
        r2="${OUT_DIR}/${base}_2.fastq.gz"
        printf "ZAF_%s\tSouth Africa\t%s\t%s\t%s\n" "${run}" "${run}" "$f" "$r2"
    done
} > "${MANIFEST}"

echo "Manifest: ${MANIFEST}"
wc -l "${MANIFEST}"

# ============================================================
# 8. Verify pairs
# ============================================================
echo
echo "------------------------------------------------------------"
echo "STEP 8 — FASTQ verification"
echo "------------------------------------------------------------"

python3 - "${OUT_DIR}" "${NIGERIA_TARGET}" "${SOUTH_AFRICA_TARGET}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
nga_t = int(sys.argv[2])
zaf_t = int(sys.argv[3])

nga = sorted(out.glob("NGA_*_1.fastq.gz"))
zaf = sorted(out.glob("ZAF_*_1.fastq.gz"))
errors = []

print(f"Nigeria R1 files : {len(nga)}")
print(f"South Africa R1  : {len(zaf)}")

for r1 in nga + zaf:
    r2 = Path(str(r1).replace("_1.fastq.gz", "_2.fastq.gz"))
    if not r2.exists():
        errors.append(f"Missing R2 for {r1.name}")
    elif r1.stat().st_size < 1000 or r2.stat().st_size < 1000:
        errors.append(f"Very small files: {r1.name}")

if len(nga) < 1 or len(zaf) < 1:
    errors.append("Need at least one Nigeria and one South Africa pair")

if errors:
    print("FASTQ VERIFICATION FAILED")
    for e in errors:
        print(f"- {e}")
    raise SystemExit(1)

print("FASTQ VERIFICATION PASSED")
print(f"Nigeria pairs: {len(nga)} (target was {nga_t})")
print(f"South Africa pairs: {len(zaf)} (target was {zaf_t})")
print("Ready for GBS-Sentinel raw-read analysis.")
PY

echo
echo "============================================================"
echo "FINAL AFRICAN RAW-READ COHORT"
echo "============================================================"
echo
echo "FASTQ directory : ${OUT_DIR}"
echo "Metadata        : ${METADATA}"
echo "Manifest        : ${MANIFEST}"
echo "Run accessions  : ${ACCESSIONS}"
echo
echo "Run pipeline (example):"
echo "  cd ${PROJECT_DIR}"
echo "  nextflow run pipeline.nf \\"
echo "    -profile docker \\"
echo "    --raw_dir african_raw_reads \\"
echo "    --outdir results_africa_raw \\"
echo "    --max_cpus 8 \\"
echo "    --max_memory 18g"
echo
echo "STATUS: READY (after verification passed)"

echo "============================================================"
