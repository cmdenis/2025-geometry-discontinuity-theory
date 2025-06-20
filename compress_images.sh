#!/bin/bash

# Usage check
if [ -z "$1" ]; then
  echo "Usage: $0 /path/to/folder [scale percentage (default 50)]"
  exit 1
fi

FOLDER="$1"
SCALE=${2:-50}  # Default to 50% if no scale is given

# Check if folder exists
if [ ! -d "$FOLDER" ]; then
  echo "Error: Folder '$FOLDER' does not exist."
  exit 1
fi

echo "Resizing images in $FOLDER to $SCALE%..."

shopt -s nullglob
for img in "$FOLDER"/*.{jpg,jpeg,png,JPG,JPEG,PNG}; do
  echo "Resizing $img"
  mogrify -resize "${SCALE}%" "$img"
done

echo "Resizing complete."