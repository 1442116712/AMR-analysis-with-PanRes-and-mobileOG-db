#### download PanRes(v1.0.2) database 
wget https://zenodo.org/records/13885013/files/panres_genes.fa?download=1

#### collect rumen (cow-rumen-v1-0-1 and sheep-rumen-v1-0) MAGs from Mgnify
wget -r -H -nd -A fna https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/cow-rumen/v1.0.1/species_catalogue/
wget -r -H -nd -A fna https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/sheep-rumen/v1.0/species_catalogue/

# remove the pan-genome file and save the original fna
rm *pan* 

# combine all the sequences
cat *.fna > MAG.fna 

# count the sequence number of the MAG.fna and split into 2 for the task 
grep -c "^>" MAG.fna
awk '/^>/ {n++} {if (n <= X) print > "spilt1.ffn"; else print > "spilt2.ffn"}' MAG.fna

# create the job
nano MAG_blast.sh

# copy the code below
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
module load apps/ncbiblast/2.15.0/gcc-14.1.0
blastn -db [PanRes_db path] -query spilt1.ffn -out spilt1_output.txt -outfmt 6
cat spilt1_output.txt | awk '{ if ($3 > 70 && $11 < 0.000001) { print $0 } }' > filter_output1.txt

# run the batch
sbatch MAG_blast.sh

# grep the MAGs accession ID by searching the blast query contig id
grep -o -w -f ID.txt *.fna | sed 's/:/ /' | awk '{print $2 "\t" $1}' | sort -u > all_MAG_matches.txt


#### collect rumen pure culture whole genomes from NCBI (should have your target genome accession ID - ID.txt) 
# download and activate the NCBI dataset tool
wget https://anaconda.org/conda-forge/ncbi-datasets-cli
module load mamba/1.3.1
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

# download the target genome from the NCBI dataset
datasets download genome accession --inputfile ID.txt --filename rumen_WGS_genomes.zip
unzip rumen_WGS_genomes.zip
cd rumen_WGS_genomes

# cp all fna files to a new directory and merge
mkdir rumen_WGS
find . -type f -name "*.fna" -exec cp -t rumen_WGS {} +
cd rumen_WGS
cat *.fna > PCG.fna

# create the job
nano PCG_blast.sh

# copy the code below
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
module load apps/ncbiblast/2.15.0/gcc-14.1.0
blastn -db [PanRes_db path] -query PCG.fna -out PCG_output.txt -outfmt 6
cat PCG_output.txt | awk '{ if ($3 > 70 && $11 < 0.000001) { print $0 } }' > filter_PCG_output.txt

# run the batch
sbatch PCG_blast.sh

# grep the MAGs accession ID by searching the blast query contig id
grep -o -w -f ID.txt *.fna | sed 's/:/ /' | awk '{print $2 "\t" $1}' | sort -u > all_PCG_matches.txt
