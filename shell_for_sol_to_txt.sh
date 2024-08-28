#!/bin/bash

# Check if a directory is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 directory_path"
  exit 1
fi

DIRECTORY=$1

# Check if the provided argument is a directory
if [ -d "$DIRECTORY" ]; then
  # Loop through the files in the directory
  for FILE in "$DIRECTORY"/*; doc
    # Check if it's a file
    if [ -f "$FILE" ]; then
      python3 OA_sol_to_txt.py $(basename "$FILE")
    fi
  done
else
  echo "$DIRECTORY is not a valid directory"
  exit 1
fi  
