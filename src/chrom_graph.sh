#!/bin/bash

# Assign command-line arguments to variables
input_bam=$1
output_graph=$2

samtools idxstats $input_bam | cut -f 1,3 | grep -v "\*" > output.txt

echo "set terminal pdf size 8,6" > plot.gp
echo "set output 'Stats/output.pdf'" >> plot.gp
echo "set style data histograms" >> plot.gp
echo "set style fill solid 1.0 border -1" >> plot.gp
echo "set xtics rotate by -45" >> plot.gp
echo "plot 'output.txt' using 2:xtic(1) title 'Number of reads aligned to each chromosome'" >> plot.gp
gnuplot plot.gp

mv Stats/output.pdf $output_graph
