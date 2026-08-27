#!/usr/bin/env bash
# Fetch the published Liu et al. (2022) artifacts used as validation reference.
#
# Both archives are DOI-pinned and immutable, so they are fetched on demand
# rather than vendored into this repository:
#
#   Zenodo  10.5281/zenodo.5823021   software: exe, inputs, MATLAB post-processing
#   Pangaea 10.1594/PANGAEA.940262   results:  Model A/B/C simulation output + mesh/
#   Paper   10.1029/2021JB023420     (paywalled; not fetched)
#
# Layout produced (matches what scripts/paleo_site_stats.py expects):
#   <dest>/zenodo_software/   <dest>/pangaea_results/
#
# Usage:  bash scripts/fetch_published_reference.sh [dest_dir]   (default archive/published)
set -euo pipefail

DEST="${1:-archive/published}"

ZENODO_DOI="10.5281/zenodo.5823021"
ZENODO_URL="https://zenodo.org/records/5823021/files/dunyuliu/Multicycle_dynamic_SSAF_NSJF-v1.0.0.zip?download=1"
ZENODO_ZIP="Multicycle_dynamic_SSAF_NSJF-v1.0.0.zip"
ZENODO_MD5="6ef08516580b86d0005b5f0490f485cb"

PANGAEA_DOI="10.1594/PANGAEA.940262"
PANGAEA_URL="https://download.pangaea.de/dataset/940262/files/mesh_simulation_result.zip"
PANGAEA_ZIP="mesh_simulation_result.zip"
PANGAEA_MD5="894584931b3916a84b3c7e311d7b7678"

have() { command -v "$1" >/dev/null 2>&1; }
have unzip || { echo "error: unzip not found" >&2; exit 1; }
if have curl; then GET() { curl -fSL --retry 3 -o "$1" "$2"; }
elif have wget; then GET() { wget -q --show-progress -O "$1" "$2"; }
else echo "error: need curl or wget" >&2; exit 1; fi

check_md5() {  # file expected_md5 -> 0 if match
  [ -f "$1" ] || return 1
  [ "$(md5sum "$1" | awk '{print $1}')" = "$2" ]
}

fetch() {  # name url zip md5 subdir doi
  local name="$1" url="$2" zip="$3" md5="$4" sub="$5" doi="$6"
  local dir="$DEST/$sub"
  mkdir -p "$dir"
  if check_md5 "$dir/$zip" "$md5"; then
    echo "$name: $zip already present and matches md5, skipping download"
  else
    echo "$name: fetching $doi"
    GET "$dir/$zip" "$url"
    if ! check_md5 "$dir/$zip" "$md5"; then
      echo "error: $name md5 mismatch." >&2
      echo "  expected $md5" >&2
      echo "  got      $(md5sum "$dir/$zip" | awk '{print $1}')" >&2
      echo "  These DOIs are immutable, so a mismatch means a corrupt or" >&2
      echo "  redirected download, not a new version. Not unpacking." >&2
      exit 1
    fi
  fi
  echo "$name: unpacking"
  unzip -q -o "$dir/$zip" -d "$dir"
}

fetch "zenodo"  "$ZENODO_URL"  "$ZENODO_ZIP"  "$ZENODO_MD5"  zenodo_software/dunyuliu "$ZENODO_DOI"
fetch "pangaea" "$PANGAEA_URL" "$PANGAEA_ZIP" "$PANGAEA_MD5" pangaea_results "$PANGAEA_DOI"

cat <<MSG

Done. Reference tree under $DEST:
  zenodo_software/  software     $ZENODO_DOI
  pangaea_results/  results      $PANGAEA_DOI

Validate the post-processing port against published Table 2:
  python3 scripts/paleo_site_stats.py $DEST/pangaea_results/work_vis7_fs0.5

Note: neither archive contains a Rate_direction.txt in the rad/s convention the
compiled binary consumes at runtime; see $DEST/README.md if present.
MSG
