#!/bin/bash
#SBATCH --job-name=PMI_analysis
#SBATCH --output=output_d%a.log
#SBATCH --error=error_d%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=144:00:00


d_values=(50 100 200)
ord_values=(1)
c_len_values=(1000)
n_values=(400 800 1600)

get_parameters() {
    local task_id=$1
    local id=$((task_id - 1))

    local n_size=${#n_values[@]}
    local c_len_size=${#c_len_values[@]}
    local ord_size=${#ord_values[@]}
    local d_size=${#d_values[@]}
    
    local d_idx=$(( id / (n_size * c_len_size * ord_size) ))
    local remaining=$(( id % (n_size * c_len_size * ord_size) ))
    local ord_idx=$(( remaining / (n_size * c_len_size) ))
    remaining=$(( remaining % (n_size * c_len_size) ))
    local c_len_idx=$(( remaining / n_size ))
    local n_idx=$(( remaining % n_size ))
    d=${d_values[$d_idx]}
    ord=${ord_values[$ord_idx]}
    c_len=${c_len_values[$c_len_idx]}
    n=${n_values[$n_idx]}
}

# Get parameters for current task
get_parameters $SLURM_ARRAY_TASK_ID

cd /home/zhiweixu/pmi

module load R/4.4.0

export R_LIBS_USER=/home/zhiweixu/R/x86_64-pc-linux-gnu-library/4.4

echo "Running job array ${SLURM_ARRAY_TASK_ID}"
echo "Parameters: d=${d}, ord=${ord}, c_len=${c_len}, n=${n}"
echo "R_LIBS_USER is set to $R_LIBS_USER"
module list

Rscript PMITable1.R ${d} ${ord} ${c_len} ${n}
