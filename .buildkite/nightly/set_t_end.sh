#!/bin/bash
# Update the `t_end`s of each config to T_END.
# Usage: set_t_end.sh T_END CONFIG...
set -euo pipefail

t_end=$1
shift
if [[ ! $t_end =~ ^[0-9]+(secs|mins|hours|days)$ ]]; then
    echo "Bad t_end '$t_end'. Examples: 600secs, 96mins, 3hours, 2days" >&2
    exit 1
fi

for config in "$@"; do
    sed -i -E "s/^t_end:.*/t_end: \"$t_end\"/" "$config"
    grep -H "^t_end:" "$config"
done
