#!/bin/bash
# Build a styled, self-contained HTML page (version 2) from an exercise Markdown file.
# Usage, from the repo folder:
#   html_build/build_v2.sh unixExs_EFNB3.md "Unix Exercises: Getting to Know the Dog Genome" "NGS course 2026 · Day 1"
set -euo pipefail

md=$1
title=$2
subtitle=${3:-"NGS course 2026"}
here=$(cd "$(dirname "$0")" && pwd)
out="${md%.md}_v2.html"

script=$(mktemp)
trap 'rm -f "$script"' EXIT
{ echo "<script>"; cat "$here/course-v2.js"; echo "</script>"; } > "$script"

pandoc "$md" -f gfm -t html5 -o "$out" \
  --standalone --embed-resources \
  --lua-filter "$here/callouts.lua" \
  --css "$here/course-v2.css" \
  --include-after-body "$script" \
  --toc --toc-depth=2 \
  --metadata title="$title" \
  --metadata subtitle="$subtitle" \
  --metadata pagetitle="$title"

echo "Wrote $out"
