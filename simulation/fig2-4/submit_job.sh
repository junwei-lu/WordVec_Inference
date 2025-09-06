#!/bin/bash

# Calculate total combinations
d_values=(50 100 200)
ord_values=(1)
c_len_values=(1000)
n_values=(400 800 1600)

# Multiply the lengths of all arrays to get total combinations
total_combinations=$((${#d_values[@]} * ${#ord_values[@]} * ${#c_len_values[@]} * ${#n_values[@]}))

# Submit the job with the calculated array size
sbatch --array=1-${total_combinations} PMITable1_241230.sh
