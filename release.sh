#!/usr/bin/env bash
# release.sh — gated release process for EQdyna.2Dcycle.
#
# What it does:
#   1. Refuses to run on a dirty tree.
#   2. Reads VERSION; refuses if v$VERSION tag already exists.
#   3. Confirms CHANGELOG.md has an entry for this version (not just [Unreleased]).
#   4. Runs test_system/smoke.py.
#   5. Creates annotated tag v$VERSION.
#   6. Prints next steps.
#
# Usage:
#   ./release.sh           # do all gates and tag
#   ./release.sh --dry-run # show what would happen, don't tag

set -euo pipefail
cd "$(dirname "$0")"

DRY_RUN=""
SKIP_SMOKE=""
for arg in "$@"; do
    case "$arg" in
        --dry-run)    DRY_RUN="(dry-run) " ;;
        --skip-smoke) SKIP_SMOKE=1 ;;
        *) echo "Unknown flag: $arg"; exit 1 ;;
    esac
done

VERSION="$(cat VERSION)"
TAG="v${VERSION}"

echo "${DRY_RUN}Releasing ${TAG}"

# 1. clean tree
if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "ERROR: working tree is dirty. Commit or stash first." >&2
    git status -s
    exit 1
fi

# 2. tag must not exist
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "ERROR: tag $TAG already exists. Bump VERSION before releasing." >&2
    exit 1
fi

# 3. CHANGELOG must mention this version
if ! grep -qE "^## \[${VERSION}\]" CHANGELOG.md; then
    echo "ERROR: CHANGELOG.md has no '## [${VERSION}] - YYYY-MM-DD' header." >&2
    echo "       Rename '## [Unreleased]' and add today's date." >&2
    exit 1
fi

# 4. docs checklist: README.md / CLAUDE.md must have been touched since last tag
PREV_TAG="$(git describe --tags --abbrev=0 2>/dev/null || true)"
if [ -n "$PREV_TAG" ]; then
    DOC_CHANGES="$(git diff --name-only "${PREV_TAG}..HEAD" -- README.md CLAUDE.md test_system/README.md work/README.md 2>/dev/null || true)"
    if [ -z "$DOC_CHANGES" ]; then
        echo "WARNING: no doc files (README.md / CLAUDE.md / sub-READMEs) changed since ${PREV_TAG}."
        echo "         Releases should keep docs in lockstep with code changes."
        if [ -z "$DRY_RUN" ]; then
            read -p "Continue anyway? [y/N] " ans
            [ "$ans" = "y" ] || [ "$ans" = "Y" ] || { echo "Aborted."; exit 1; }
        else
            echo "(dry-run) would prompt for confirmation"
        fi
    else
        echo "Docs touched since ${PREV_TAG}: $(echo "$DOC_CHANGES" | tr '\n' ' ')"
    fi
fi

# 5. smoke test
if [ -n "$SKIP_SMOKE" ]; then
    echo "WARNING: --skip-smoke given; smoke test bypassed."
elif [ -z "$DRY_RUN" ]; then
    echo "Running smoke test..."
    python3 test_system/smoke.py
else
    echo "(dry-run) skipping smoke test"
fi

# 6. tag
NOTES="$(awk "/^## \[${VERSION}\]/,/^## \[/" CHANGELOG.md | sed '$d')"
if [ -z "$DRY_RUN" ]; then
    git tag -a "$TAG" -m "Release ${TAG}

${NOTES}"
    echo "Created tag ${TAG}"
else
    echo "(dry-run) would create annotated tag ${TAG} with notes:"
    echo "---"
    echo "$NOTES"
    echo "---"
fi

echo
echo "Next steps:"
echo "  git push origin main"
echo "  git push origin ${TAG}"
echo "  gh release create ${TAG} --notes-from-tag   # optional GitHub release"
