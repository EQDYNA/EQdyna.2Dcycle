#! /bin/bash
# install.sh — build the binary and export env vars.
#
# Usage:
#   ./install.sh                  # autodetect OS, build + install python deps
#   ./install.sh -m <machine>     # build only (skip python deps)
#   ./install.sh -e <machine>     # build + install python deps
#   ./install.sh -c <machine>     # configure only (skip build)
#   source ./install.sh           # just export EQDYNA2DCYCLEROOT + PATH
#
# Supported machines: ubuntu, macos, ls6 (TACC), grace (TAMU)

set -e

MACH=""
ENV=""
CONFIG=""

while getopts "m:e:c:h" OPTION; do
    case $OPTION in
        m) MACH=$OPTARG ;;
        e) MACH=$OPTARG; ENV="True" ;;
        c) MACH=$OPTARG; CONFIG="True" ;;
        h)
            sed -n '2,12p' "$0"; exit 0 ;;
    esac
done

# Autodetect OS if no flag was given
if [ -z "$MACH" ] && [ "$0" = "$BASH_SOURCE" ]; then
    case "$(uname -s)" in
        Linux)  MACH=ubuntu; ENV="True" ;;
        Darwin) MACH=macos;  ENV="True" ;;
        *)      echo "Cannot autodetect OS. Run with -m <ubuntu|macos|ls6|grace>." ; exit 1 ;;
    esac
    echo "Autodetected: -e $MACH"
fi

if [ -n "$MACH" ]; then
    export MACHINE=$MACH

    if [ -n "$ENV" ]; then
        case "$MACHINE" in
            ubuntu)
                echo "Installing python deps (Ubuntu)..."
                pip3 install --user --quiet numpy matplotlib xarray gmsh meshio nbconvert \
                    || echo "WARNING: pip3 install failed; install numpy/matplotlib/xarray manually"
                ;;
            macos)
                echo "Installing python deps (macOS via Homebrew)..."
                brew install gcc python
                pip3 install --break-system-packages numpy matplotlib xarray gmsh meshio nbconvert
                ;;
            ls6|grace)
                echo "Skipping python deps on $MACHINE (use module load)."
                ;;
        esac
    fi

    if [ -z "$CONFIG" ]; then
        echo "Building EQdyna.2Dcycle on $MACHINE..."
        ( cd src && make )
        mkdir -p bin
        # rm first so a busy binary (running sim) can be replaced safely:
        # unlink keeps the in-use inode alive for running procs; mv lays a fresh file.
        rm -f bin/run_eqdyna2d_* 2>/dev/null || true
        mv src/run_eqdyna2d_* bin/
        # Clean .o/.mod artifacts; src/ holds source only.
        rm -f src/*.o src/*.mod
    fi

    chmod -R 755 scripts
fi

# Always export env (works whether sourced or executed)
export EQDYNA2DCYCLEROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export PATH="$EQDYNA2DCYCLEROOT/bin:$EQDYNA2DCYCLEROOT/scripts:$PATH"

echo
echo "EQDYNA2DCYCLEROOT=$EQDYNA2DCYCLEROOT"
echo "Binary: $(ls "$EQDYNA2DCYCLEROOT/bin/" 2>/dev/null | grep run_eqdyna2d_ | tail -1)"
echo
echo "Next: bash example_workflow.sh   # one-step demo"
echo "Or:   create.newcase --work_dir work/my_case --compset paper.saf.A"
