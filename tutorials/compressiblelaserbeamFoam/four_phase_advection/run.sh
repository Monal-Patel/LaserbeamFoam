#!/bin/bash
set -e

cd "$(dirname "$0")"

rm -rf 0 processor* log.* constant/polyMesh

cp -r 0.orig 0

echo "=== Running blockMesh ==="
blockMesh > log.blockMesh 2>&1

echo "=== Initializing fields ==="
python3 init_fields.py

echo "=== Running compressibleLaserbeamFoam ==="
compressibleLaserbeamFoam > log.compressibleLaserbeamFoam 2>&1

echo "=== Done ==="
tail -5 log.compressibleLaserbeamFoam
