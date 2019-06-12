#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --partition=euan,mrivas,normal,owners
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --err process_cnv_input.err
#SBATCH -x sh-113-15

module load py-pandas/0.23.0_py27
module load py-numpy/1.14.3_py27
module load py-scikit-image/0.15.0_py27

python process_cnv_input.py > process_cnv_input.log 

