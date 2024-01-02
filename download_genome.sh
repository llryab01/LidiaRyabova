wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
wget http://ftp.ensembl.org/pub/release104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Unzip the downloaded files
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.104.gtf.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

# Make the script executable
#chmod +x download_genome.sh

# Run the script to download and unzip the files
#./download_genome.sh
