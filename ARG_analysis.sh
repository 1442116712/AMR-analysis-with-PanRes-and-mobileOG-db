#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --job-name=mag.blast
#SBATCH --error=stderr.txt
#SBATCH --output=stdout.txt
#SBATCH --partition=k2-medpri
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=zwu18@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
roary -f ./Achromobacter_xylosoxidans_output_i80 -e -n -v -i 80 -cd 95 ./Achromobacter_xylosoxidans_gff/*.gff
