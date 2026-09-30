#!/bin/bash
# Directory to clear
TARGET_DIR="utils/starlink_tles"
# Delete all files and subdirectories, including hidden ones
find "$TARGET_DIR" -maxdepth 1 -type f -exec rm -f {} \;
echo "All TLE files in $TARGET_DIR have been deleted."