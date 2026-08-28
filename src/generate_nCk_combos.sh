#!/bin/bash
# For use in sex chromosome assembly CBS pipeline
# Shell script version of generate_nCk_combos.py

# Usage: generate_nCk_combos.sh {COMBO_SIZE} {ID1} {ID2} ...

# Recursive function to generate k-choose combos
combinations() {
  local k="$1"
  local arr=("${@:2}")

  if (( k == 0 )); then
    echo ""
    return
  fi

  local n="${#arr[@]}"
  if (( k > n )); then
    return
  fi

  for (( i=0; i<=n-k; i++ )); do
    local start="${arr[i]}"
    local rest=("${arr[@]:i+1}")

    combinations $((k-1)) "${rest[@]}" | while IFS= read -r end; do
     if [[ -z "$end" ]]; then
       echo "$start"
     else
       echo "$start $end"
     fi
    done
  done
}

# Run the function
combinations "$1" "${@:2}"
