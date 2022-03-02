#! /bin/bash

gnuplot gnu_plot

AIMDIR="../OUTPUT/KMC_temperature_change_OUTPUT/N150/"

#change name of output file
for file in *.png
do
mv "$file" "N150_10d-6_$file"
mv "N150_10d-6_$file" "$AIMDIR"
done;

mv "stat_moves.DATA" "N150_10d-6_stat_moves.DATA"
cp "N150_10d-6_stat_moves.DATA" "$AIMDIR"
