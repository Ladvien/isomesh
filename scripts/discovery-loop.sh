#!/usr/bin/env bash
# ~/isomesh/scripts/discovery-loop.sh
set -euo pipefail
cd ~/isomesh
PROMPT=docs/research/2026-09-12-discovery-loop-prompt.md
LEDGER=docs/research/questions.md
MAX_ITER=${1:-100}
for i in $(seq 1 "$MAX_ITER"); do
  claude -p "$(sed -n '/^## The prompt/,/^## Wrapper/p' "$PROMPT")" \
    --allowedTools "Bash,Read,Edit,Write,mcp__home-still__*" \
    --max-turns 200 || true
  if [ -f "$LEDGER" ] && grep -q '^ASK_USER' "$LEDGER"; then
    echo "iteration $i: agent asked for the user"; grep '^ASK_USER' "$LEDGER"; exit 2
  fi
  echo "iteration $i: $(git log -1 --pretty=%s)"
  sleep 5
done
echo "reached MAX_ITER=$MAX_ITER"
