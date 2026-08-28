#!/usr/bin/env bash
# Resume the EQdyna.2Dcycle Claude Code session.
#   ./resume_claude.sh              # resume the pinned session below
#   ./resume_claude.sh <session-id> # resume a different session
set -euo pipefail

SESSION_ID="${1:-cbbebfd4-f0c9-4986-9423-2de74831f444}"

cd "$(dirname "$0")"
exec claude --resume "$SESSION_ID" --dangerously-skip-permissions
