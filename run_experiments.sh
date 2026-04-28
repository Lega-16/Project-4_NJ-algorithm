#!/bin/bash

echo "archivo,tiempo_nj" > tiempos_nj.csv

for f in unique_distance_matrices/*.phy; do
    echo "Procesando $f"
    tiempo=$(/usr/bin/time -f "%e" python3 run_nj.py "$f" 2>&1 >/dev/null)
    echo "$(basename $f),$tiempo" >> tiempos_nj.csv
done
