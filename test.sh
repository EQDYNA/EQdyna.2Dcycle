#! /bin/bash
# Usage: bash test.sh <safpub|saf|subei|gulang>
#
# safpub  -> compset/test.saf.published  (C_mesh=2, validates against 2.0.2 baseline)
# saf     -> compset/test.saf            (C_mesh=3, gmsh-based; should match 2.0.2)
# subei   -> compset/test.subei          (C_mesh=3, gmsh-based)
# gulang  -> compset/test.gulang         (C_mesh=3, gmsh-based)

set -e

if [ $# -ne 1 ]; then
    echo "Usage: bash test.sh <safpub|saf|subei|gulang>"
    exit 1
fi

case "$1" in
    safpub)
        COMPSET="test.saf.published"
        ;;
    saf)
        COMPSET="test.saf"
        ;;
    subei)
        COMPSET="test.subei"
        ;;
    gulang)
        COMPSET="test.gulang"
        ;;
    *)
        echo "Error: unknown case '$1'. Must be safpub, saf, subei, or gulang."
        exit 1
        ;;
esac

CASEDIR="$(pwd)/work/test.$1"

# Set up environment
source "$(dirname "$0")/install.eqdyna.2dcycle.sh"

echo "=== Creating case: $COMPSET -> $CASEDIR ==="
create.newcase "$CASEDIR" "$COMPSET"

cd "$CASEDIR"
echo "=== Running case.setup ==="
python3 case.setup

echo "=== Running simulation ==="
bash run.sh
