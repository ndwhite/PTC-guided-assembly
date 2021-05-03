# PTC-guided-assembly
Guided assembly pipeline for phototransduction cascade probe set work. IN REVIEW. Decoy sequences available here: 10.5281/zenodo.4734556. In order to create hte species-specific references, need to concatenate the decoy sequences with the species references included here. "Swarm" commands are for a slurm hpc manager, but give you an idea of resources needed.

## Trimming:	Trimmomatic v0.39
I toss unpaired reads, but you could use them.

    #swarm -g 12 -t 16 -f Trimmo.sh --module trimmomatic/0.39
    java -jar $TRIMMOJAR PE -threads 16 -trimlog trimlog_SAMPLE.log -summary sumstats_SAMPLE.txt -phred33 SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz SAMPLE_R1-trimmed.fastq.gz SAMPLE_R1-unpaired.fastq.gz SAMPLE_R2-trimmed.fastq.gz SAMPLE_R2-unpaired.fastq.gz ILLUMINACLIP:All_adapters.fa:2:30:10:1:TRUE


## Assembly:	Trinity v2.8.4
May need to increase ram as high as 128G for very large files.

    #swarm -g 30 -t 12 -f Trinity.sh --module trinity --time 3-00:00:00
    Trinity --seqType fq --left SAMPLE_R1-trimmed.fastq.gz --right SAMPLE_R2-trimmed.fastq.gz --CPU 12 --max_memory 28G --time=3-00:00:00 --output SAMPLE_PTC_temp_trinity

Then, remove whitespace from headers.


## Read mapping:	BWA-MEM v0.7.17

First create index:

    #swarm --module bwa/0.7.17 --partition=quick -f index.sh
    bwa index -p Chicken Chicken_exome_with_decoys.fa

Map single-end reads:

    #swarm --module bwa/0.7.17 --partition=quick,norm -f BWA.sh
    bwa mem Chicken SAMPLE-trimmed.fastq.gz -B 8 > SAMPLE_Chicken.sam

Map paired-end reads:

    #swarm --module bwa/0.7.17 --partition=quick,norm -f BWA.sh
    bwa mem Chicken SAMPLE_R1-trimmed.fastq.gz SAMPLE_R2-trimmed.fastq.gz -B 8 > SAMPLE_Chicken.sam


## Masking and creating contigs:	Samtools v1.11; Bedtools v2.29.2
Altered from Ryan Schott's script (Schott RK, et al. 2017. Targeted capture of complete coding regions across divergent species. Genome Biology and Evolution 9: 398–414).

    module load samtools/1.11
    module load bedtools/2.29.2
    samtools sort -O bam -o Sorted_SAMPLE_Chicken.bam SAMPLE_Chicken.sam
    samtools view -b -F 260 Sorted_SAMPLE_Chicken.bam > Mapped_SAMPLE_Chicken.bam
    samtools depth -aa -d 1000000 Sorted_SAMPLE_Chicken.bam > depth_SAMPLE_Chicken.txt
    cat depth_SAMPLE_Chicken.txt | awk '$3==0 {print}' | awk -v OFS='\t' '{print $1,int($2)-1,$2}' > zero_coverage_SAMPLE_Chicken.bed
    bedtools maskfasta -fi Chicken_exome_with_decoys.fa -bed zero_coverage_SAMPLE_Chicken.bed -fo zero_masked_ref_SAMPLE_Chicken.fas
    bcftools mpileup -Ou -d 1000000 -f zero_masked_ref_SAMPLE_Chicken.fas Mapped_SAMPLE_Chicken.bam | bcftools call -Ou -mv | bcftools norm -Oz -f zero_masked_ref_SAMPLE_Chicken.fas > normalized_calls_zero_masked_SAMPLE_Chicken.zcf.gz
    bcftools index normalized_calls_zero_masked_SAMPLE_Chicken.zcf.gz
    bcftools consensus -f zero_masked_ref_SAMPLE_Chicken.fas normalized_calls_zero_masked_SAMPLE_Chicken.zcf.gz >   bcftools_consensus_zero_masked_SAMPLE_Chicken.fas
    module load python/3.7
    python depth_of_coverage_fast.py -t depth_SAMPLE_Chicken.txt -o  Depth_stats_SAMPLE_Chicken_FAST.csv 



## Alignment:	MASCE v2.03
If wanted, concatenate exons into CDS (with something like Cat_exons_to_CDS.pl). Then, de-interleave contigs (with something like DE-interleave_this_fasta.pl). Then do the following several steps:

Align SAG seqs --seq = reference seq (more reliable), --seq-lr = our sequences

    java -jar macse_v2.03.jar -prog alignSequences -out_AA 3.aligned_profiles/AA/SAG_all.fa -out_NT 3.aligned_profiles/NT/SAG_all.fa -seq 0.raw_input/ref_seqs/SAG.fa -seq_lr 1.all_seqs/SAG.fa

Remove the reference seq by subsetting to a list that excludes the ref seq

    java -jar macse_v2.03.jar -prog splitAlignment -align 3.aligned_profiles/NT/SAG_all.fa -out_subset 3.aligned_profiles/NT/SAG_all_filter_1.fa -subset 3.aligned_profiles/filter_lists/SAG_filter_1.txt

Trim to remove bases that do not have alignments now that ref seq is removed

    java -jar macse_v2.03.jar -prog trimAlignment -align 3.aligned_profiles/NT/SAG_all_filter_1.fa -out_NT 3.aligned_profiles/NT/SAG_all_filter_1_trimmed.fa

Export alignment

    java -jar macse_v2.03.jar -prog exportAlignment -align 3.aligned_profiles/NT/SAG_all_filter_1_trimmed.fa -codonForExternalFS NNN -codonForFinalStop NNN -codonForInternalFS NNN -codonForInternalStop NNN -out_AA 4.seqs_prepared/all_seqs/AA/SAG.fa -out_AA_consensus 4.seqs_prepared/all_seqs/AA/Consensus/SAG.fa -out_NT 4.seqs_prepared/all_seqs/NT/SAG.fa -out_NT_consensus 4.seqs_prepared/all_seqs/NT/Consensus/SAG.fa -out_stat_per_seq 4.seqs_prepared/all_seqs/Stats/SAGseq_stat.txt -out_stat_per_site 4.seqs_prepared/all_seqs/Stats/SAGsite_stat.txt


