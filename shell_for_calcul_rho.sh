#!/bin/bash

# Check if the file path is provided as an argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <file_path>"
  exit 1
fi

FILE=$1

# Check if the file exists
if [ ! -f "$FILE" ]; then
  echo "File not found!"
  exit 1
fi

# Read and print lines that start with a number
while IFS= read -r line; do
  if [[ $line =~ ^[0-9] ]]; then
    echo "$line"
  fi
done < "$FILE"
