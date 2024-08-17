#!/bin/bash

# Check if a file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <json-file>"
    exit 1
fi

# Pretty print the JSON file using Python's json.tool
python3 -m json.tool "$1" | less
