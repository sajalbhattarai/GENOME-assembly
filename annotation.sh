#!/bin/bash -l

#SBATCH --time=4:00:00
#SBATCH -A scholar
#SBATCH -N 1 -n 20
#SBATCH --ntasks 20
#SBATCH --output=/depot/lindems-class/peptide/annotations.out
#SBATCH --error=/depot/lindems-class/peptide/annotations.error.out
#SBATCH --job-name <FILL IN>

#---- MAKE SURE TO CHECK DIRECTORIES FOR OUT/ERROR FILES AS WELL AS INPUT FILES ----

# I highly recommend using Video Code Studio (free download) for terminal commands, cluster access, and slurm creation
# Upload this slurm to the same directory as your input file (our amino_acid file from RAST)
# Sbatch this slurm from that directory 

#load bioinfo to run modules below BUT bioinfo may not be available in the future so biocontainers may need used (ml biocontainers)
ml bioinfo


#STEP 1: Mummer (if needed)

ml MUMmer
ml gnuplot

cd /depot/lindems-class/peptide/gln/HW5

mummer -mum -b -c /depot/lindems-class/data/Final.assembly.fasta \
 /depot/lindems-class/peptide/gln/HW#4/scaffolds.fasta > mummer.mums

mummerplot --postscript --prefix MUMmer mummer.mums

gnuplot MUMmer.gp

# Check out file (annotations.out) to make sure this ran properly. Should print the below statement when it runs properly.
echo "Finished MUMmer" 
#------------------------------------------------------------------------------------------------------------------------------------


#STEP 2: TIGRfam and Pfam

ml HMMER

#run tigrfam 
hmmsearch -o my_filtered_rast.TIGR --tblout my_rast.TIGR.tsv --cut_tc ../../../data/databases/TIGRfam/TIGRFAMs_15.0_HMM.LIB amino_acids.faa

#run pfam
hmmsearch -o my_filtered_rast.Pfam --tblout my_rast.Pfam.tsv --cut_tc ../../../data/databases/Pfam_36/Pfam-A.hmm amino_acids.faa 

# Check out file (annotations.out) to make sure this ran properly. Should print the below statement when it runs properly.
echo "Finished TIGRfam and Pfam" 
#------------------------------------------------------------------------------------------------------------------------------------


#STEP 3: MEROPS and TCDB

ml bioinfo diamond

diamond blastp --sensitive -d /depot/lindems-class/data/databases/MEROPES/meropscan -q amino_acids.faa -o MEROPS.txt

diamond blastp --sensitive -d /depot/lindems-class/data/databases/TCDB/tcdb -q amino_acids.faa -o TCDB.txt

# Check out file (annotations.out) to make sure this ran properly. Should print the below statement when it runs properly.
echo "Finished MEROPS and TCDB" 
#------------------------------------------------------------------------------------------------------------------------------------


#STEP 4: DBCAN and CAZy

ml anaconda
conda create -n dbcan4 python=3.8 dbcan -c conda-forge -c bioconda
conda activate dbcan4

#steps above: load anaconda (a package depository for python programs, similar to CRAN)
#create a new env with python 3.8, and download dbcan from either conda-forge or bioconda channel
#activate new env

run_dbcan amino_acids.faa protein --out_dir /depot/lindems-class/peptide/gln/HW7/dbcan4_result --db_dir /depot/lindems-class/data/databases/CAZyme_V12 --tools all --dia_cpu 10 --hmm_cpu 10 --dbcan_thread 10

# Check out file (annotations.out) to make sure this ran properly. Should print the below statement when it runs properly.
echo "Finished dbCAN and CAZyme" 
#------------------------------------------------------------------------------------------------------------------------------------


#STEP 5: TIGR Descriptions

# After TIGR is ran online, use this to loop through the database and tag your TIGR.txt with the descriptions provided by the professors

cd /depot/lindems-class/data/databases/TIGRfam/TIGRFAMs_15.0_INFO

# Simple command to read lines in a file called "tigr" (whatever you call your output). 
# The function then goes through all of the .INFO files online and prints the descriptions from line "DE" to OUTFILE.txt.
# OUTFILE.txt is then joined to your TIGR output file. Mine is called "TIGR_Outputs_clean.txt". Needed to "clean" the file of some Windows OS nonsense that Unix doesn't like.

while read tigr
do
awk '$1=="DE" {print}' $(echo $tigr".INFO") >> OUTFILE.txt
done < TIGR_Outputs_clean.txt 
paste TIGR_Outputs_clean.txt OUTFILE.txt

# Check out file (annotations.out) to make sure this ran properly. Should print the below statement when it runs properly.
echo "TIGR Loop"
#------------------------------------------------------------------------------------------------------------------------------------

# Message me with any questions or if anything is missing!