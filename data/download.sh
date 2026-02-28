#!/usr/bin/env bash
# =============================================================================
# Terra Incognita — Phase 1 Reference Data Downloader
# P1.1: data acquisition for 5 terrain class reference regions
#
# Usage:
#   ./data/download.sh [merit|geomorpho90m|koppen|verify|all]
#
# Prerequisites:
#   - sha256sum
#
# Both MERIT-DEM and Geomorpho90m must be downloaded MANUALLY (automated
# download is no longer supported for either dataset — see data/sources.md).
# Place extracted TIF files in data/raw/merit/ and data/raw/geomorpho90m/
# before running this script. Köppen-Geiger still downloads automatically.
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/raw"

mkdir -p "${RAW_DIR}/merit"
mkdir -p "${RAW_DIR}/geomorpho90m"
mkdir -p "${RAW_DIR}/koppen"

# ─── Helpers ──────────────────────────────────────────────────────────────────

log()  { echo "[download.sh] $*"; }
warn() { echo "[download.sh] WARNING: $*" >&2; }
die()  { echo "[download.sh] ERROR: $*" >&2; exit 1; }

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1. Install it and retry."
}

# Download a URL to a destination path. Skips if file already exists.
fetch() {
    local url="$1" dest="$2"
    if [[ -f "${dest}" ]]; then
        log "Already exists, skipping: $(basename "${dest}")"
        return 0
    fi
    log "Downloading: $(basename "${dest}")"
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 --retry-delay 5 -o "${dest}.tmp" "${url}" && mv "${dest}.tmp" "${dest}"
    else
        wget -q --tries=3 -O "${dest}.tmp" "${url}" && mv "${dest}.tmp" "${dest}"
    fi
}

# Download a GDAL COG sub-window for a given bounding box.
# Uses /vsicurl/ for on-the-fly streaming — no need to download full global file.
fetch_cog_clip() {
    local url="$1" dest="$2"
    local min_lon="$3" max_lon="$4" min_lat="$5" max_lat="$6"
    if [[ -f "${dest}" ]]; then
        log "Already exists, skipping: $(basename "${dest}")"
        return 0
    fi
    require_cmd gdal_translate
    log "Clipping COG region (${min_lat}..${max_lat}°N ${min_lon}..${max_lon}°E): $(basename "${dest}")"
    # gdal_translate -projwin: <ulx uly lrx lry> = <left top right bottom>
    gdal_translate -q \
        -projwin "${min_lon}" "${max_lat}" "${max_lon}" "${min_lat}" \
        "/vsicurl/${url}" "${dest}"
}

# ─── Köppen-Geiger ────────────────────────────────────────────────────────────
download_koppen() {
    log "=== Downloading Köppen-Geiger classification ==="
    # Beck et al. (2023) v3 — 1991-2020 present-day 1km global map
    # Available from Figshare: https://figshare.com/articles/dataset/21789074
    # The zip contains multiple period subfolders and resolutions; extract the
    # 1991_2020/koppen_geiger_0p00833333.tif and place it directly in raw/koppen/:
    #   unzip -j koppen_geiger_tif.zip "1991_2020/koppen_geiger_0p00833333.tif" -d raw/koppen/
    local KOPPEN_ZIP="${RAW_DIR}/koppen/koppen_geiger_tif.zip"
    local KOPPEN_TIF="${RAW_DIR}/koppen/koppen_geiger_0p00833333.tif"

    # NOTE: Figshare direct download URLs are versioned. Verify the file ID below
    # at https://figshare.com/articles/dataset/21789074 before running.
    local KOPPEN_URL="https://figshare.com/ndownloader/articles/21789074"

    if [[ ! -f "${KOPPEN_TIF}" ]]; then
        fetch "${KOPPEN_URL}" "${KOPPEN_ZIP}"
        log "Extracting Köppen archive..."
        cd "${RAW_DIR}/koppen" && unzip -q "$(basename "${KOPPEN_ZIP}")" && cd "${SCRIPT_DIR}"
    else
        log "Köppen TIF already present."
    fi
}

# ─── Geomorpho90m ─────────────────────────────────────────────────────────────
download_geomorpho90m() {
    log "=== Geomorpho90m: checking for manually downloaded tiles ==="
    # Automated download is no longer viable (OpenTopography S3 returns 403).
    # Download the 6 archives below manually — see data/sources.md §1.2 for instructions.
    # Each archive is .tar.gz; sampler extracts internally. Do NOT pre-extract.
    # Place archives in: ${RAW_DIR}/geomorpho90m/

    local all_present=true
    # Archive → region:  n30e060,n30e090 → Himalaya
    #                    n00e000         → Congo + Ahaggar (non-overlapping bboxes)
    #                    n30w120         → Colorado Plateau
    #                    n30w090         → Atlantic Coastal
    #                    s30w090         → Ucayali (Amazon headwaters, second FluvialHumid)
    for tile_id in n30e060 n30e090 n00e000 n30w120 n30w090 s30w090; do
        local archive="${RAW_DIR}/geomorpho90m/geom_90M_${tile_id}.tar.gz"
        if [[ -f "${archive}" ]]; then
            log "Found: geom_90M_${tile_id}.tar.gz"
        else
            warn "Missing: ${archive}"
            all_present=false
        fi
    done

    if ! "${all_present}"; then
        die "One or more Geomorpho90m archives missing. See data/sources.md §1.2 for download instructions."
    fi
    log "All Geomorpho90m archives present."
}

# ─── MERIT-DEM ────────────────────────────────────────────────────────────────
download_merit() {
    log "=== MERIT-DEM: checking for manually downloaded tiles ==="
    # MERIT-DEM automated download is no longer supported in this script.
    # Download the 6 archives below manually — see data/sources.md §1.1 for instructions.
    # Each archive is .tar (not .tar.gz); sampler extracts internally. Do NOT pre-extract.
    # Place archives in: ${RAW_DIR}/merit/

    local all_present=true
    # Archive → region:  n30e060,n30e090 → Himalaya
    #                    n00e000         → Congo + Ahaggar (non-overlapping bboxes)
    #                    n30w120         → Colorado Plateau
    #                    n30w090         → Atlantic Coastal
    #                    s30w090         → Ucayali (Amazon headwaters, second FluvialHumid)
    for tile_id in n30e060 n30e090 n00e000 n30w120 n30w090 s30w090; do
        local archive="${RAW_DIR}/merit/dem_tif_${tile_id}.tar"
        if [[ -f "${archive}" ]]; then
            log "Found: dem_tif_${tile_id}.tar"
        else
            warn "Missing: ${archive}"
            all_present=false
        fi
    done

    if ! "${all_present}"; then
        die "One or more MERIT-DEM archives missing. See data/sources.md §1.1 for download instructions."
    fi
    log "All MERIT-DEM archives present."
}

# ─── Checksum verification ────────────────────────────────────────────────────
verify_checksums() {
    log "=== Verifying checksums ==="
    require_cmd sha256sum

    local checksum_file="${SCRIPT_DIR}/checksums.txt"
    local all_ok=true

    while IFS='  ' read -r expected_hash filepath; do
        # Skip comments and PENDING entries
        [[ "${expected_hash}" == \#* ]] && continue
        [[ "${expected_hash}" == "PENDING" ]] && continue
        [[ -z "${expected_hash}" ]] && continue

        local full_path="${SCRIPT_DIR}/${filepath}"
        if [[ ! -f "${full_path}" ]]; then
            warn "File not found: ${filepath}"
            all_ok=false
            continue
        fi

        local actual_hash
        actual_hash=$(sha256sum "${full_path}" | awk '{print $1}')
        if [[ "${actual_hash}" != "${expected_hash}" ]]; then
            warn "Checksum MISMATCH for ${filepath}"
            warn "  Expected: ${expected_hash}"
            warn "  Actual:   ${actual_hash}"
            all_ok=false
        else
            log "OK: ${filepath}"
        fi
    done < "${checksum_file}"

    if "${all_ok}"; then
        log "All checksums verified."
    else
        die "Checksum verification failed. Re-run downloads for failed files."
    fi
}

# Update checksums.txt with actual hashes after download
update_checksums() {
    log "=== Updating checksums.txt with computed hashes ==="
    require_cmd sha256sum
    local checksum_file="${SCRIPT_DIR}/checksums.txt"

    # Replace PENDING entries with actual hashes
    while read -r line; do
        if [[ "${line}" == PENDING* ]]; then
            local filepath="${line#PENDING  }"
            local full_path="${SCRIPT_DIR}/${filepath}"
            if [[ -f "${full_path}" ]]; then
                local hash
                hash=$(sha256sum "${full_path}" | awk '{print $1}')
                echo "${hash}  ${filepath}"
            else
                echo "${line}  # FILE MISSING"
            fi
        else
            echo "${line}"
        fi
    done < "${checksum_file}" > "${checksum_file}.tmp"
    mv "${checksum_file}.tmp" "${checksum_file}"
    log "checksums.txt updated."
}

# ─── Main ─────────────────────────────────────────────────────────────────────
CMD="${1:-all}"

case "${CMD}" in
    koppen)        download_koppen ;;
    geomorpho90m)  download_geomorpho90m ;;
    merit)         download_merit ;;
    verify)        verify_checksums ;;
    update-checksums) update_checksums ;;
    all)
        download_koppen
        download_geomorpho90m
        download_merit
        update_checksums
        verify_checksums
        log "=== All Phase 1 data acquisition complete ==="
        ;;
    *)
        echo "Usage: $0 [merit|geomorpho90m|koppen|verify|update-checksums|all]"
        exit 1
        ;;
esac
