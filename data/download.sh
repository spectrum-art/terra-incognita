#!/usr/bin/env bash
# =============================================================================
# Terra Incognita — Phase 1 Reference Data Downloader
# P1.1: data acquisition for 5 terrain class reference regions
#
# Usage:
#   ./data/download.sh [merit|geomorpho90m|koppen|verify|all]
#
# Prerequisites:
#   - curl or wget
#   - gdal_translate  (for Geomorpho90m COG clips)
#   - sha256sum
#
# MERIT-DEM requires registration at:
#   https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/
# Set credentials via environment variables BEFORE running the merit step:
#   export MERIT_USER="your@email.com"
#   export MERIT_PASS="yourpassword"
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
    # Beck et al. (2018) — present-day 1km global map
    # Available from Figshare: https://figshare.com/articles/6396959
    # The zip contains Beck_KG_V1_present_0p0083.tif (~350MB)
    local KOPPEN_ZIP="${RAW_DIR}/koppen/Beck_KG_V1_present.zip"
    local KOPPEN_TIF="${RAW_DIR}/koppen/Beck_KG_V1_present_0p0083.tif"

    # NOTE: Figshare direct download URLs are versioned. Verify the URL below
    # at https://figshare.com/articles/dataset/6396959 before running.
    local KOPPEN_URL="https://figshare.com/ndownloader/files/12407516"

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
    log "=== Downloading Geomorpho90m geomorphon tiles ==="
    require_cmd gdal_translate

    # Geomorpho90m is served as COG tiles on OpenTopography.
    # Base URL — verify at https://portal.opentopography.org/dataspace/dataset?opentopoID=OTDS.012020.4326.1
    # before running. The path below is the expected S3/CDN path; it may need updating.
    local BASE_URL="https://opentopography.s3.sdsc.edu/dataspace/OTDS.012020.4326.1/raster/geom"

    # Tile naming: geomorpho90m_geom_merit_dem_egm96_wgs84_{TILE}_c_20181019.tif
    # 15°×15° tiles named by SW corner (N/S + lat, E/W + lon, 2-digit lat, 3-digit lon)
    # Example: N15E075 = 15-30°N, 75-90°E

    # Function to fetch one geomorphon tile by its 15° tile ID
    fetch_geom_tile() {
        local tile_id="$1"
        local filename="geomorpho90m_geom_merit_dem_egm96_wgs84_${tile_id}_c_20181019.tif"
        local url="${BASE_URL}/${filename}"
        local dest="${RAW_DIR}/geomorpho90m/${filename}"
        fetch "${url}" "${dest}"
    }

    # Tiles needed to cover all 5 regions (see data/regions.json)
    # Alpine (Himalayas 25-35°N 78-92°E) → 15° tiles: N15E075, N15E090
    fetch_geom_tile "N15E075"
    fetch_geom_tile "N15E090"
    # FluvialHumid (Congo 8°S-2°N 15-28°E) → S15E015, S15E030
    fetch_geom_tile "S15E015"
    fetch_geom_tile "S15E030"
    # Cratonic (Ahaggar 20-27°N 4-12°E) → N15E000 (covers 15-30°N 0-15°E)
    fetch_geom_tile "N15E000"
    # FluvialArid (Colorado Plateau 35-38°N 107-113°W) → N30W120, N30W105
    fetch_geom_tile "N30W120"
    fetch_geom_tile "N30W105"
    # Coastal (Atlantic 34-38°N 75-82°W) → N30W090, N30W075
    fetch_geom_tile "N30W090"
    fetch_geom_tile "N30W075"

    log "Geomorpho90m download complete."
    warn "If URLs returned 404, verify current tile paths at the OpenTopography portal"
    warn "and update BASE_URL in this script accordingly."
}

# ─── MERIT-DEM ────────────────────────────────────────────────────────────────
download_merit() {
    log "=== Downloading MERIT-DEM tiles ==="

    if [[ -z "${MERIT_USER:-}" ]] || [[ -z "${MERIT_PASS:-}" ]]; then
        die "MERIT-DEM requires credentials. Set MERIT_USER and MERIT_PASS environment variables.
     Register at: https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/"
    fi

    # MERIT-DEM download URL pattern (post-registration).
    # Each tile is a tar archive containing [tile_id]_dem.tif
    # URL format verified at time of writing; re-check if downloads fail.
    local BASE_URL="https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/data"

    fetch_merit_tile() {
        local tile_id="$1"
        local tar_file="${RAW_DIR}/merit/${tile_id}_dem.tar"
        local tif_file="${RAW_DIR}/merit/${tile_id}_dem.tif"

        if [[ -f "${tif_file}" ]]; then
            log "Already extracted: ${tile_id}_dem.tif"
            return 0
        fi

        local url="${BASE_URL}/${tile_id}/${tile_id}_dem.tar"
        log "Downloading MERIT tile: ${tile_id}"
        if command -v curl >/dev/null 2>&1; then
            curl -fL --retry 3 --user "${MERIT_USER}:${MERIT_PASS}" \
                -o "${tar_file}.tmp" "${url}" && mv "${tar_file}.tmp" "${tar_file}"
        else
            wget -q --tries=3 --user="${MERIT_USER}" --password="${MERIT_PASS}" \
                -O "${tar_file}.tmp" "${url}" && mv "${tar_file}.tmp" "${tar_file}"
        fi

        log "Extracting ${tile_id}..."
        tar -xf "${tar_file}" -C "${RAW_DIR}/merit/" && rm "${tar_file}"
    }

    # ── Alpine: Central Himalayas ──────────────────────────────────────────────
    fetch_merit_tile "n25e080"
    fetch_merit_tile "n25e085"
    fetch_merit_tile "n30e080"
    fetch_merit_tile "n30e085"

    # ── FluvialHumid: Congo Basin margins ─────────────────────────────────────
    fetch_merit_tile "s00e015"
    fetch_merit_tile "s00e020"
    fetch_merit_tile "s05e015"
    fetch_merit_tile "s05e020"
    fetch_merit_tile "s05e025"

    # ── Cratonic: Ahaggar Massif ───────────────────────────────────────────────
    fetch_merit_tile "n20e005"
    fetch_merit_tile "n20e010"
    fetch_merit_tile "n25e005"

    # ── FluvialArid: Colorado Plateau ─────────────────────────────────────────
    fetch_merit_tile "n35w115"
    fetch_merit_tile "n35w110"

    # ── Coastal: US Atlantic coastal plain ────────────────────────────────────
    fetch_merit_tile "n35w080"
    fetch_merit_tile "n30w080"
    fetch_merit_tile "n35w085"

    log "MERIT-DEM download complete."
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
