###### RESISTANCE AMPLICONS ######

#1 Rename the files to be the correct sample names

#2 Create a sample_file.csv
ls *.fastq.gz | sed 's/.fastq.gz//' > sample_file.csv
### Make sure to add 'sample' to head of column now ###

# 3 Check bed file matches annotation in the gff file:
bedtools intersect -a /mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf3d7_drugresamplicons.bed -b /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.modified.new.gff3 -wa | sort | uniq > matched.bed
comm -23 <(cut -f1-4 /mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf3d7_drugresamplicons.bed | sort) matched.bed
# this prints all regions that match between the bed file and the gff file to matched.bed and then prints any 
# that do not match in the terminal for you to see. You want the output to terminal to be blank, which indicates that all of the bed file regions are also in the gff.
# Check that regions in the bedfile match the reference fastq file
awk 'NR==FNR {chrlen[$1]=$2; next} $1 in chrlen && $3 <= chrlen[$1] {next} {print "Invalid:", $0}' /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.fasta.fai /mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf3d7_drugresamplicons.bed
# If nothing prints to terminal, you're good.


# 7 Make sure there is no samclip as this gets rid of all amplicon reads
# Mapping of amplicon data
python /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/plasmodium_capeverde/sophie_nanopore_amplicon_script_minimap2_v2_pf.py \
    --index-file sample_file.csv \
    --ref /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.fasta \
    --gff /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.modified.new.gff3 \
    --bed /mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf3d7_drugresamplicons.bed \
    --min-base-qual 20 \
    --threads 10 > amplicon_log.txt 2>&1 \
    --snpeff-db Plasmodium_falciparum

# snpeff did not work so annotating separately using bcftools csq

bcftools view snps.vcf.gz | bcftools csq -p a \
-f /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.fasta \
-g /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.modified.new.gff3 \
-Oz -o snps.ann.vcf.gz 2> csq.parse.log

(echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tTBCSQ" ; bcftools query -f "[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%AD\t%DP\t%TBCSQ\\n]" snps.ann.vcf.gz) > combined_snps_bcftoolscsq_trans.txt

# now running for indels

bcftools view --threads 10 -v indels combined.genotyped_filtered_FMTDP10.vcf.gz -Oz -o indels.vcf.gz

bcftools view indels.vcf.gz | bcftools csq -p a \
-f /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.fasta \
-g /mnt/storage11/sophie/reference_genomes/pfal_3D7_ref_genome/Pfalciparum.genome.modified.new.gff3 \
-Oz -o indels.ann.vcf.gz 2> indels.csq.parse.log

(echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tTBCSQ" ; bcftools query -f "[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%AD\t%DP\t%TBCSQ\\n]" indels.ann.vcf.gz) > combined_indels_bcftoolscsq_trans.txt


#8 Create amplicon x sample coverage matrix by combining coverage files (the sample_coverage_mean.txt files) to see how good they are across the sample set
# Use script 2.create_ampliconxsample

# These amplicons had very high coverage, and the filter criteria were different to the malaria profiler results
# Therefore, I am downsampling the bam files to be maximum read depth of 1000 per amplicon
# Then, I will pass these through malaria profiler using specific criteria:
# Any variants that fail QC, I will check separately.

# Downsample the reads
python /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/plasmodium_capeverde/downsample_to_1000reads.py \
--highcov highcov_amplicons.tsv \
--bed /mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf3d7_drugresamplicons.bed \
--bam-dir . \
--out-dir downsampled_bams \
--max-reads 1000 \
--seed 47 \
> downsample_to_1000reads.log 2>&1

# Bams are now in downsampled_bams folder, and any amplicons which were not too high coverage have not been downsampled, but these amplicons remain in the bams.

# Check the coverage of the new downsampled bams
for bam in *.bam; do
    sample=${bam%.bam}
    echo "Processing $sample ..."
    bedtools coverage -a pf3d7_drugresamplicons.bed -b "$bam" -mean > "${sample}_coverage_mean.txt"
done

# BAMs have downsampled correctly.
# Now running BAMS through the malaria profiler with edited parameters, depth = 10, AF = 0.01,0.05, so now any variants with AF > 0.05 will be kept, and those between 0.01 and 0.05 will fail QC
# Installed malaria profiler using Jody's GitHub

conda activate malaria-profiler

for bam in *.bam; do
    sample=$(basename "$bam" .bam)

    malaria-profiler \
        profile \
        -a "$bam" \
        -p "$sample" \
        -t 10 \
        --txt \
        --depth 0,10 \
        --af 0.01,0.05 \
        --strand 0,3 \
        --resistance_db Plasmodium_falciparum
done

malaria-profiler collate 



