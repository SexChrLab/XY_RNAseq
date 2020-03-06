import os

configfile: "blood.config.json"

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
bwa_path = "bwa"
samtools_path = "samtools"
bbduksh_path = "bbduk.sh"
hisat_path = "hisat2"
bamtools_path = "bamtools"
picard_path = "picard"
featureCounts_path = "featureCounts"
trimmomatic_path = "trimmomatic"
star_path = "star"

#
REF_TYPE = ["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]
REF_TYPE_HISAT = ["Ref_GRCh38_Y_HardMasked_HISAT_index","Ref_GRCh38_Y_PARsMasked_HISAT_index"]

# Directory
fastq_directory = "fastq_files/"

rule all:
    input: 
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_XX.sam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_XY.sam", male_sample = config["males"]),
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_XX_DEF.sam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_XY_DEF.sam", male_sample = config["males"]),

        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_XX.bam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_XY.bam", male_sample = config["males"]),
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_XX_DEF.bam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_XY_DEF.bam", male_sample = config["males"]),
        
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_XX.bam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_XY.bam", male_sample = config["males"]),
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_XX_DEF.bam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_XY_DEF.bam", male_sample = config["males"]),
        
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam", male_sample = config["males"]),
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam", female_sample = config["females"]),
        expand("HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_mkdup_XY_DEF.bam", male_sample = config["males"]),
        expand("HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_mkdup_XX_DEF.bam", female_sample = config["females"]),

        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam", male_sample = config["males"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam", female_sample = config["females"]),
        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam", male_sample = config["males"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam", female_sample = config["females"]),

        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam.bai", male_sample = config["males"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam.bai", female_sample = config["females"]),
        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam.bai", male_sample = config["males"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam.bai", female_sample = config["females"]),
     
        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX.txt", female_sample = config["females"]),
        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX_DEF.txt", female_sample = config["females"]),

        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chr8.bam", male_sample = config["males"]),
        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chrX.bam", male_sample = config["males"]),
        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chrY.bam", male_sample = config["males"]),

        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chr8.bam", male_sample = config["males"]),
        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrX.bam", male_sample = config["males"]),
        expand("HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrY.bam", male_sample = config["males"]),

        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chr8.bam", female_sample = config["females"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chrX.bam", female_sample = config["females"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chrY.bam", female_sample = config["females"]),

        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chr8.bam", female_sample = config["females"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrX.bam", female_sample = config["females"]),
        expand("HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrY.bam", female_sample = config["females"]),

        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX.txt", female_sample = config["females"]),
        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX_DEF.txt", female_sample = config["females"]),

        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX.txt", female_sample = config["females"]),
        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX_DEF.txt", female_sample = config["females"]),

        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX.txt", female_sample = config["females"]),
        expand("stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX_DEF.txt", female_sample = config["females"]),
                 
        expand("FeatureCounts/HISAT/{male_sample}_HISAT_geneCounts_XY.txt", male_sample = config["males"]),
        expand("FeatureCounts/HISAT/{female_sample}_HISAT_geneCounts_XX.txt", female_sample = config["females"]),
        expand("FeatureCounts/HISAT/{male_sample}_HISAT_geneCounts_XY_DEF.txt", male_sample = config["males"]),
        expand("FeatureCounts/HISAT/{female_sample}_HISAT_geneCounts_XX_DEF.txt", female_sample = config["females"]),
        
        expand("STAR/tmp/{female_sample}_STAR_pair_trim_XX.bam", female_sample = config["females"]),
        expand("STAR/tmp/{male_sample}_STAR_pair_trim_XY.bam", male_sample = config["males"]),
        expand("STAR/tmp/{female_sample}_STAR_pair_trim_XX_DEF.bam", female_sample = config["females"]),
        expand("STAR/tmp/{male_sample}_STAR_pair_trim_XY_DEF.bam", male_sample = config["males"]),
        
        expand("STAR/tmp/{female_sample}_STAR_pair_trim_sort_XX.bam", female_sample = config["females"]),
        expand("STAR/tmp/{male_sample}_STAR_pair_trim_sort_XY.bam", male_sample = config["males"]),
        expand("STAR/tmp/{female_sample}_STAR_pair_trim_sort_XX_DEF.bam", female_sample = config["females"]),
        expand("STAR/tmp/{male_sample}_STAR_pair_trim_sort_XY_DEF.bam", male_sample = config["males"]),
     
        expand("STAR/tmp/{male_sample}_STAR_pair_trim_sort_mkdup_XY.bam", male_sample = config["males"]),
        expand("STAR/tmp/{female_sample}_STAR_pair_trim_sort_mkdup_XX.bam", female_sample = config["females"]),
        expand("STAR/tmp/{male_sample}_STAR_pair_trim_sort_mkdup_XY_DEF.bam", male_sample = config["males"]),
        expand("STAR/tmp/{female_sample}_STAR_pair_trim_sort_mkdup_XX_DEF.bam", female_sample = config["females"]),

        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam", male_sample = config["males"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam", female_sample = config["females"]),
        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam", male_sample = config["males"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam", female_sample = config["females"]),

        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam.bai", male_sample = config["males"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam.bai", female_sample = config["females"]),
        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam.bai", male_sample = config["males"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam.bai", female_sample = config["females"]),

        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chr8.bam", male_sample = config["males"]),
        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chrX.bam", male_sample = config["males"]),
        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chrY.bam", male_sample = config["males"]),

        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chr8.bam", male_sample = config["males"]),
        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrX.bam", male_sample = config["males"]),
        expand("STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrY.bam", male_sample = config["males"]),

        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chr8.bam", female_sample = config["females"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chrX.bam", female_sample = config["females"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chrY.bam", female_sample = config["females"]),

        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chr8.bam", female_sample = config["females"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrX.bam", female_sample = config["females"]),
        expand("STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrY.bam", female_sample = config["females"]),
        
        expand("stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY.txt", male_sample = config["males"]),
        expand("stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX.txt", female_sample = config["females"]),
        expand("stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX_DEF.txt", female_sample = config["females"]),

        expand("stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY.txt", male_sample = config["males"]),
        expand("stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX.txt", female_sample = config["females"]),
        expand("stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX_DEF.txt", female_sample = config["females"]),

        expand("stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY.txt", male_sample = config["males"]),
        expand("stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX.txt", female_sample = config["females"]),
        expand("stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY_DEF.txt", male_sample = config["males"]),
        expand("stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX_DEF.txt", female_sample = config["females"]),        

        expand("FeatureCounts/STAR/{male_sample}_STAR_geneCounts_XY.txt", male_sample = config["males"]),
        expand("FeatureCounts/STAR/{female_sample}_STAR_geneCounts_XX.txt", female_sample = config["females"]),
        expand("FeatureCounts/STAR/{male_sample}_STAR_geneCounts_XY_DEF.txt", male_sample = config["males"]),
        expand("FeatureCounts/STAR/{female_sample}_STAR_geneCounts_XX_DEF.txt", female_sample = config["females"]),

    input:
        expand("fastq_files/{sample_name}_1.fastq", sample_name = config["sample_names"]),
        expand("fastq_files/{sample_name}_2.fastq", sample_name = config["sample_names"]),
        expand("trimmed_fastqs/{sample_name}_trimmed_1.fastq", sample_name = config["sample_names"]),
        expand("trimmed_fastqs/{sample_name}_trimmed_2.fastq", sample_name = config["sample_names"])
    input:
        expand("fastq_files/{sample_name}_1.fastq", sample_name = config["sample_names"])
    input:
        expand("refs/{ref_type}.fa.fai", ref_type = REF_TYPE),
        expand("refs/{ref_type}.fa.amb", ref_type = REF_TYPE),
        expand("refs/{ref_type}.dict", ref_type = REF_TYPE)

rule prep_refs_mk_sy_ln:
    input:
        ref = lambda wildcards: config["genome_paths"][wildcards.ref_type]
    output:
        ref_sy_ln = "refs/{ref_type}.fa"
    shell:
        """
        ln -s {input.ref} {output.ref_sy_ln}
        """

rule prep_refs:
    input:
        ref = "refs/{ref_type}.fa"
    output:
        fai = "refs/{ref_type}.fa.fai",
        amb = "refs/{ref_type}.fa.amb",
        dict = "refs/{ref_type}.dict"
    params:
        samtools = samtools_path,
        bwa = bwa_path
    run:
        # faidx
        shell("{params.samtools} faidx {input.ref}")

        # .dict
        shell("{params.samtools} dict -o {output.dict} {input.ref}")

        # bwa
        shell("{params.bwa} index {input.ref}")


rule mk_sy_ln_fastqs:
    input:
        original_1 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq1"],
        original_2 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq2"]
    output:
        R1_out = "fastq_files/{sample_name}_1.fastq",
        R2_out = "fastq_files/{sample_name}_2.fastq"
    shell:
        """
        ln -s {input.original_1} {output.R1_out};
        ln -s {input.original_2} {output.R2_out}
        """

rule trimmomatic:
    input:
        fq1 = lambda wildcards: os.path.join(fastq_directory, config[wildcards.sample_name]["fq1_sy"]),
        fq2 = lambda wildcards: os.path.join(fastq_directory, config[wildcards.sample_name]["fq2_sy"]),
        ADAPTER_FASTA = "/mnt/storage/SAYRES/GTEx/DE_genomeFilters/00_tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10"
    output:
        out_fq1 = "trimmed_fastqs/{sample_name}_trimmed_1.fastq",
        out_fq2 = "trimmed_fastqs/{sample_name}_trimmed_2.fastq",
        out_fq1_unpair = "trimmed_fastqs/{sample_name}_trimmed_unpaired_1.fastq",
        out_fq2_unpair = "trimmed_fastqs/{sample_name}_trimmed_unpaired_2.fastq",
        logfile = "trimmed_fastqs/logfiles/{sample_name}_trimmomatic.log"
    params:
        threads = 4,
        leading = 10,
        trailing = 25,
        winsize = 4,
        winqual = 25,
        minlen = 40   
    shell:
        "trimmomatic PE -phred33 -threads {params.threads} -trimlog {output.logfile} "
        "{input.fq1} {input.fq2} {output.out_fq1} {output.out_fq1_unpair} "
        "{output.out_fq2} {output.out_fq2_unpair} "
        "ILLUMINACLIP:{input.ADAPTER_FASTA} "
        "LEADING:{params.leading} TRAILING:{params.trailing} "
        "SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}"

rule HISAT_paired_males:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{male_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{male_sample}_trimmed_2.fastq"
    output:
        out_1 = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY.sam"
    params:
        HISAT_Index_male = config["HG38_Transcriptome_Index_HISAT_Path_male"],
    shell:
        "hisat2 --dta -q --phred33 -p 8 -x {params.HISAT_Index_male} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"

rule HISAT_paired_females:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{female_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{female_sample}_trimmed_2.fastq"
    output:
        out_1 = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX.sam"
    params:
        HISAT_Index_female = config["HG38_Transcriptome_Index_HISAT_Path_female"],
    shell:
        "hisat2 --dta -q --phred33 -p 8 -x {params.HISAT_Index_female} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"

# males
rule HISAT_samtools_view_males:
    input:
        SAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY.sam"
    output:
        BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY.bam"
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule HISAT_bam_sort_males:    
    input:
        IN_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY.bam"
    output:
        sort_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_XY.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule HISAT_MarkDups_males:
    input:
        sort_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_XY.bam"
    output:
        BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam",
        metrics = "stats/HISAT/PICARD/{male_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_AddReadGrps_males:
    input:
        Read_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam"
    output:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam",
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_index_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        BAI = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule HISAT_stats_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_samSub_chr8_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        subBAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule HISAT_stats_chr8_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chr8.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_samSub_chrX_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        subBAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule HISAT_stats_chrX_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chrX.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_samSub_chrY_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        subBAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule HISAT_stats_chrY_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_chrY.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_feautreCounts_gene_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        counts = "FeatureCounts/HISAT/{male_sample}_HISAT_geneCounts_XY.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

        
# females
rule HISAT_samtools_view_females:
    input:
        SAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX.sam"
    output:
        BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX.bam"
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule HISAT_bam_sort_females:  
    input:
        IN_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX.bam"
    output:
        sort_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_XX.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule HISAT_MarkDups_females:
    input:
        sort_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_XX.bam"
    output:
        BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam",
        metrics = "stats/HISAT/PICARD/{female_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_AddReadGrps_females:
    input:
        Read_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam"
    output:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam",
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_index_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        BAI = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule HISAT_stats_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"
 
rule HISAT_samSub_chr8_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        subBAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule HISAT_stats_chr8_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chr8.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_samSub_chrX_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        subBAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule HISAT_stats_chrX_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chrX.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_samSub_chrY_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        subBAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule HISAT_stats_chrY_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_chrY.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}" 
 
rule HISAT_feautreCounts_gene_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        counts = "FeatureCounts/HISAT/{female_sample}_HISAT_geneCounts_XX.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"


# HISAT DEFAULT
rule HISAT_DEF_paired_males:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{male_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{male_sample}_trimmed_2.fastq"
    output:
        out_1 = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY_DEF.sam"
    params:
        HISAT_Index_DEF = config["HG38_Transcriptome_Index_HISAT_Path_DEF"],
    shell:
        "hisat2 --dta -q --phred33 -p 8 -x {params.HISAT_Index_DEF} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"

rule HISAT_DEF_paired_females:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{female_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{female_sample}_trimmed_2.fastq"
    output:
        out_1 = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX_DEF.sam"
    params:
        HISAT_Index_DEF = config["HG38_Transcriptome_Index_HISAT_Path_DEF"],
    shell:
        "hisat2 --dta -q --phred33 -p 8 -x {params.HISAT_Index_DEF} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"

# males
rule HISAT_DEF_samtools_view_males:
    input:
        SAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY_DEF.sam"
    output:
        BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY_DEF.bam"
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule HISAT_DEF_bam_sort_males:    
    input:
        IN_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_XY_DEF.bam"
    output:
        sort_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_XY_DEF.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule HISAT_DEF_MarkDups_males:
    input:
        sort_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_XY_DEF.bam"
    output:
        BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_mkdup_XY_DEF.bam",
        metrics = "stats/HISAT/PICARD/{male_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_DEF_AddReadGrps_males:
    input:
        Read_BAM = "HISAT/tmp/{male_sample}_HISAT_pair_trim_sort_mkdup_XY_DEF.bam"
    output:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam",
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_DEF_index_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        BAI = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule HISAT_DEF_stats_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule HISAT_DEF_samSub_chr8_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        subBAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule HISAT_DEF_stats_chr8_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chr8.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_DEF_samSub_chrX_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        subBAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule HISAT_DEF_stats_chrX_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrX.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_DEF_samSub_chrY_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        subBAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule HISAT_DEF_stats_chrY_bam_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrY.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule HISAT_DEF_DEF_feautreCounts_gene_males:
    input:
        BAM = "HISAT/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        counts = "FeatureCounts/HISAT/{male_sample}_HISAT_geneCounts_XY_DEF.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

        
# females
rule HISAT_DEF_samtools_view_females:
    input:
        SAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX_DEF.sam"
    output:
        BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX_DEF.bam"
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule HISAT_DEF_bam_sort_females:  
    input:
        IN_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_XX_DEF.bam"
    output:
        sort_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_XX_DEF.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule HISAT_DEF_MarkDups_females:
    input:
        sort_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_XX_DEF.bam"
    output:
        BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_mkdup_XX_DEF.bam",
        metrics = "stats/HISAT/PICARD/{female_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_DEF_AddReadGrps_females:
    input:
        Read_BAM = "HISAT/tmp/{female_sample}_HISAT_pair_trim_sort_mkdup_XX_DEF.bam"
    output:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam",
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule HISAT_DEF_index_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        BAI = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule HISAT_DEF_stats_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_DEF_samSub_chr8_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        subBAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule HISAT_DEF_stats_chr8_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chr8.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_DEF_samSub_chrX_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        subBAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule HISAT_DEF_stats_chrX_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrX.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule HISAT_DEF_samSub_chrY_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        subBAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule HISAT_DEF_stats_chrY_bam_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrY.bam"
    output:
        stats = "stats/HISAT/BAM_stats/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

 
rule HISAT_DEF_feautreCounts_gene_females:
    input:
        BAM = "HISAT/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        counts = "FeatureCounts/HISAT/{female_sample}_HISAT_geneCounts_XX_DEF.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

#--------------
#  STAR
#--------------
rule STAR_paired_males:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{male_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{male_sample}_trimmed_2.fastq",
    output:
        out_1 = "STAR/tmp/{male_sample}_STAR_pair_trim_XY.bam"
    params:
        STAR_Index_male = config["HG38_Transcriptome_Index_STAR_Path_male"],
        STAR_DEF_GTF = config["ensemble.gtf"],
    shell:
        "STAR --genomeDir {params.STAR_Index_male} --sjdbGTFfile {params.STAR_DEF_GTF} --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn {input.Trimmed_FASTQ1} {input.Trimmed_FASTQ2} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1} --runThreadN 8"

rule STAR_paired_females:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{female_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{female_sample}_trimmed_2.fastq"
    output:
        out_1 = "STAR/tmp/{female_sample}_STAR_pair_trim_XX.bam"
    params:
        STAR_Index_female = config["HG38_Transcriptome_Index_STAR_Path_female"],
        STAR_DEF_GTF = config["ensemble.gtf"],
    shell:
        "STAR --genomeDir {params.STAR_Index_female} --sjdbGTFfile {params.STAR_DEF_GTF} --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn {input.Trimmed_FASTQ1} {input.Trimmed_FASTQ2} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1} --runThreadN 8"

rule STAR_bam_sort_males:    
    input:
        IN_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_XY.bam"
    output:
        sort_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_XY.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule STAR_MarkDups_males:
    input:
        sort_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_XY.bam"
    output:
        BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_mkdup_XY.bam",
        metrics = "stats/STAR/PICARD/{male_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule STAR_AddReadGrps_males:
    input:
        Read_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_mkdup_XY.bam"
    output:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam",
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule STAR_index_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        BAI = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule STAR_stats_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule STAR_samSub_chr8_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        subBAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule STAR_stats_chr8_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chr8.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_samSub_chrX_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        subBAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule STAR_stats_chrX_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chrX.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_samSub_chrY_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        subBAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule STAR_stats_chrY_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_chrY.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule STAR_feautreCounts_gene_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        counts = "FeatureCounts/STAR/{male_sample}_STAR_geneCounts_XY.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

        
# females

rule STAR_bam_sort_females:  
    input:
        IN_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_XX.bam"
    output:
        sort_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_XX.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule STAR_MarkDups_females:
    input:
        sort_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_XX.bam"
    output:
        BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_mkdup_XX.bam",
        metrics = "stats/STAR/PICARD/{female_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule STAR_AddReadGrps_females:
    input:
        Read_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_mkdup_XX.bam"
    output:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam",
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule STAR_index_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        BAI = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule STAR_stats_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule STAR_samSub_chr8_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        subBAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule STAR_stats_chr8_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chr8.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_samSub_chrX_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        subBAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule STAR_stats_chrX_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chrX.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_samSub_chrY_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        subBAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule STAR_stats_chrY_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_chrY.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"
 
rule STAR_feautreCounts_gene_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        counts = "FeatureCounts/STAR/{female_sample}_STAR_geneCounts_XX.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"


#------- DEFAULT

rule STAR_DEF_STAR_paired_males:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{male_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{male_sample}_trimmed_2.fastq"
    output:
        out_1 = "STAR/tmp/{male_sample}_STAR_pair_trim_XY_DEF.bam"
    params:
        STAR_Index_DEF = config["HG38_Transcriptome_Index_STAR_Path_DEF"],
        STAR_DEF_GTF = config["ensemble.gtf"],
    shell:
        "STAR --genomeDir {params.STAR_Index_DEF} --sjdbGTFfile {params.STAR_DEF_GTF} --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn {input.Trimmed_FASTQ1} {input.Trimmed_FASTQ2} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1} --runThreadN 8"

rule STAR_DEF_STAR_paired_females:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{female_sample}_trimmed_1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{female_sample}_trimmed_2.fastq"
    output:
        out_1 = "STAR/tmp/{female_sample}_STAR_pair_trim_XX_DEF.bam"
    params:
        STAR_Index_DEF = config["HG38_Transcriptome_Index_STAR_Path_DEF"],
        STAR_DEF_GTF = config["ensemble.gtf"],
    shell:
        "STAR --genomeDir {params.STAR_Index_DEF} --sjdbGTFfile {params.STAR_DEF_GTF} --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn {input.Trimmed_FASTQ1} {input.Trimmed_FASTQ2} --outSAMtype BAM Unsorted --outFileNamePrefix {output.out_1} --runThreadN 8"

rule STAR_DEF_bam_sort_males:    
    input:
        IN_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_XY_DEF.bam"
    output:
        sort_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_XY_DEF.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"


rule STAR_DEF_MarkDups_males:
    input:
        sort_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_XY_DEF.bam"
    output:
        BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_mkdup_XY_DEF.bam",
        metrics = "stats/STAR/PICARD/{male_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule STAR_DEF_AddReadGrps_males:
    input:
        Read_BAM = "STAR/tmp/{male_sample}_STAR_pair_trim_sort_mkdup_XY_DEF.bam"
    output:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam",
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule STAR_DEF_index_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        BAI = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule STAR_DEF_stats_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule STAR_DEF_samSub_chr8_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        subBAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule STAR_DEF_stats_chr8_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chr8.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_DEF_samSub_chrX_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        subBAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule STAR_DEF_stats_chrX_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrX.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_DEF_samSub_chrY_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        subBAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule STAR_DEF_stats_chrY_bam_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF_chrY.bam"
    output:
        stats = "stats/STAR/BAM_stats/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XY_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_DEF_feautreCounts_gene_males:
    input:
        BAM = "STAR/{male_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XY_DEF.bam"
    output:
        counts = "FeatureCounts/STAR/{male_sample}_STAR_geneCounts_XY_DEF.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

        
# females
rule STAR_DEF_bam_sort_females:  
    input:
        IN_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_XX_DEF.bam"
    output:
        sort_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_XX_DEF.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule STAR_DEF_MarkDups_females:
    input:
        sort_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_XX_DEF.bam"
    output:
        BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_mkdup_XX_DEF.bam",
        metrics = "stats/STAR/PICARD/{female_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule STAR_DEF_AddReadGrps_females:
    input:
        Read_BAM = "STAR/tmp/{female_sample}_STAR_pair_trim_sort_mkdup_XX_DEF.bam"
    output:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam",
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule STAR_DEF_index_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        BAI = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule STAR_DEF_stats_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"
 
rule STAR_DEF_samSub_chr8_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        subBAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chr8.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} 8 > {output.subBAM}"

rule STAR_DEF_stats_chr8_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chr8.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chr8_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_DEF_samSub_chrX_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        subBAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} X > {output.subBAM}"

rule STAR_DEF_stats_chrX_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrX.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrX_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"


rule STAR_DEF_samSub_chrY_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        subBAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.BAM} Y > {output.subBAM}"

rule STAR_DEF_stats_chrY_bam_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF_chrY.bam"
    output:
        stats = "stats/STAR/BAM_stats/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_stats_chrY_XX_DEF.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule STAR_DEF_feautreCounts_gene_females:
    input:
        BAM = "STAR/{female_sample}_STAR_pair_trim_sort_mkdup_rdgrp_XX_DEF.bam"
    output:
        counts = "FeatureCounts/STAR/{female_sample}_STAR_geneCounts_XX_DEF.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 1 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

