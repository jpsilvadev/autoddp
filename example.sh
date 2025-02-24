#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 {short|long}"
    exit 1
fi

# set the source directory based on the argument
if [ "$1" == "short" ]; then
    SRC_DIR="examples/short"
elif [ "$1" == "long" ]; then
    SRC_DIR="examples/long"
else
    echo "Invalid argument: $1"
    echo "Usage: $0 {short|long}"
    exit 1
fi


TARGET_DIR="./tmp"
mkdir -p "$TARGET_DIR"
cp -r "$SRC_DIR"/* "$TARGET_DIR"/

cp main.py "$TARGET_DIR"/
python3 "$TARGET_DIR/main.py"
