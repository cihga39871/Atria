
# Genome information
# Arabidopsis thaliana (thale cress) reference genome TAIR10.1
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

# Simulation program information
# ART (Skewer modified version)
# ART: a next-generation sequencing read simulator. Bioinformatics. 2012 Feb 15; 28(4): 593–594.
# Simulated from a real public sequence: SRR7243169.1


working_dir=~/analysis/atria-benchmark/simulate

mkdir -p $working_dir
cd $working_dir

# download genome
mkdir genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz -O genomes/TAIR10.1.fasta.gz
pigz -d genomes/TAIR10.1.fasta.gz

### simulate data
# Download from https://sourceforge.net/projects/skewer/files/Simulator/
cd ~/analysis/atria-benchmark/simulate/ART/art_profiler_illumina

# download real data to generate profiles
fastq-dump --split-files --origfmt SRR7243169.1
# generate profiles
./Illumina_readprofile_art profile_SRR7243169 . fastq

for i in `seq 1 3`
do
    echo "Replicate $i ----------------------------------"
    bash $atria/benchmark/simulate-run-bench.bash
done


julia $atria/benchmark/replicates-stats.jl replicate_*/summary.AdapterRemoval
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.Atria
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.Skewer
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.TrimGalore
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.Trimmomatic
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.Ktrim
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.Fastp
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.Atropos
julia $atria/benchmark/replicates-stats.jl replicate_*/summary.SeqPurge

awk 'BEGIN {print "Trimmer\tTP\tFP_ft\tFP_ot\tFN_fr\tFN_ut\tTN\tPPV\tSen.\tSpec.\tmCC" }; NR%3==2{FNR="\t"; sub("stats.summary.", "", FILENAME); print FILENAME"\t"$0}' stats.summary.* > performance_stats.df.txt

julia $atria/benchmark/replicates-stats.jl replicate_*/time_benchmark.df.txt
julia $atria/benchmark/replicates-stats.jl replicate_*/time_benchmark_gz.df.txt

julia $atria/benchmark/replicates-stats.jl replicate_*/time_benchmark.new.df.txt
julia $atria/benchmark/replicates-stats.jl replicate_*/time_benchmark.new_gz.df.txt
