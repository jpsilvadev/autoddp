#!/bin/bash

SRC_DIR="examples/standard"
TARGET_DIR="./run_example"

mkdir -p "$TARGET_DIR"
cp -r "$SRC_DIR"/* "$TARGET_DIR"/

cp autoddp.py "$TARGET_DIR"/
cd "$TARGET_DIR"
python3 autoddp.py --receptor receptor.pdbqt --config conf.txt --ligands ligands.sdf --pH 7.4 --complexes 10

