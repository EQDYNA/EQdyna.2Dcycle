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
SKIP_DOC_CHECK=""
for arg in "$@"; do
    case "$arg" in
        --dry-run)        DRY_RUN="(dry-run) " ;;
        --skip-smoke)     SKIP_SMOKE=1 ;;
        --skip-doc-check) SKIP_DOC_CHECK=1 ;;
        *) echo "Unknown flag: $arg"; exit 1 ;;
    esac
done

VERSION="$(cat VERSION)"
TAG="v${VERSION}"

echo "${DRY_RUN}Releasing ${TAG}"

# 1. docs lockstep — README.md / CLAUDE.md / sub-READMEs must be touched
#    (modified or newly added) since the last tag. Hard-fail (with
#    --skip-doc-check to bypass), so docs never silently fall behind code.
PREV_TAG="$(git describe --tags --abbrev=0 2>/dev/null || true)"
if [ -n "$PREV_TAG" ] && [ -z "$SKIP_DOC_CHECK" ]; then
    DOC_PATHS="README.md CLAUDE.md test_system/README.md work/README.md"
    DOC_PATHS="$DOC_PATHS $(find compset -maxdepth 2 -name 'README.md' 2>/dev/null)"
    DOC_CHANGES="$(git log --name-only --pretty=format: "${PREV_TAG}..HEAD" -- $DOC_PATHS 2>/dev/null \
                   | sort -u | grep -v '^$' || true)"
    if [ -z "$DOC_CHANGES" ]; then
        echo "ERROR: no doc files touched since ${PREV_TAG}." >&2
        echo "       Update README.md / CLAUDE.md / sub-READMEs to reflect what" >&2
        echo "       changed in this release, or pass --skip-doc-check to bypass." >&2
        exit 1
    fi
    echo "Docs touched since ${PREV_TAG}:"
    echo "$DOC_CHANGES" | sed 's/^/  /'
fi

# 2. clean tree
if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "ERROR: working tree is dirty. Commit or stash first." >&2
    git status -s
    exit 1
fi

# 3. tag must not exist
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "ERROR: tag $TAG already exists. Bump VERSION before releasing." >&2
    exit 1
fi

# 4. CHANGELOG must mention this version
if ! grep -qE "^## \[${VERSION}\]" CHANGELOG.md; then
    echo "ERROR: CHANGELOG.md has no '## [${VERSION}] - YYYY-MM-DD' header." >&2
    echo "       Rename '## [Unreleased]' and add today's date." >&2
    exit 1
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
# Extract this version's CHANGELOG section for the tag message.
# The obvious awk range /^## \[$VERSION\]/,/^## \[/ does NOT work: the
# start line itself matches the end pattern, so awk closes the range on
# that same line and the notes come out empty. Every tag through
# v2.0.7-rc7 carries an empty message for exactly this reason. Skip the
# heading, then stop at the next one.
NOTES="$(awk -v hdr="## [${VERSION}]" '
    index($0, hdr) == 1 { found = 1; next }
    found && /^## \[/ { exit }
    found { print }
' CHANGELOG.md)"
if [ -z "$(echo "$NOTES" | tr -d '[:space:]')" ]; then
    echo "ERROR: CHANGELOG.md has no body under '## [${VERSION}]'." >&2
    exit 1
fi
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
