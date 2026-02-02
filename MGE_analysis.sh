#download and activate mobileOG-db
module load mamba/1.3.1 
mamba create -n mobileOG-db python=3.6.15 (then restart session)
mamba activate mobileOG-db
mamba install -c conda-forge biopython
mamba install -c bioconda prodigal
mamba install -c bioconda diamond
mamba install -c anaconda pandas
mkdir mobileOG-db_workdir
cd mobileOG-db_workdir
# create the job
nano mobileOG.sh
# copy and paste the code below
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
module load mamba/1.3.1
mamba activate mobileOG-db
./mobileOGs-pl-kyanite.sh -i PCG.fna -d mobileOG-db-beatrix-1.X.dmnd -m mobileOG-db-beatrix-1.X-All.csv -k 15 -e 1e-20 -p 90 -q 90 # 'X' is based on the version you use
# run the job
sbatch mobileOG.sh
