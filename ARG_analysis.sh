# download PanRes(v1.0.2) database 
wget https://zenodo.org/records/13885013/files/panres_genes.fa?download=1
# collect rumen (cow-rumen-v1-0-1 and sheep-rumen-v1-0) MAGs from Mgnify
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
grep -o -w -f ID.txt *.fna | sed 's/:/ /' | awk '{print $2 "\t" $1}' | sort -u > all_matches.txt
