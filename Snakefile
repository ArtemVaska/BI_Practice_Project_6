URL_y00_1="http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941816/SRR941816.fastq.gz"
URL_y00_2="http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941817/SRR941817.fastq.gz"
URL_y30_1="http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941818/SRR941818.fastq.gz"
URL_y30_2="http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941819/SRR941819.fastq.gz"
URL_ref_fna="http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
URL_ref_gff="http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz"


rule sample_download:  # downloads sequencing results
	output:
		"raw_data/{sample}.fastq.gz"
	params:
	    URL=lambda wildcards: globals()["URL_{sample}".format(sample=wildcards.sample)]
	shell:
		"wget -O {output} {params.URL}"  # -O allows you to save results to a file


rule reference_download:  # downloads reference
	output:
		"reference/Saccharomyces_cerevisiae.fna.gz"
	shell:
		"wget -O {output} {URL_ref_fna}"


rule gff_download:  # downloads annotation
	output:
		"reference/Saccharomyces_cerevisiae.gff.gz"
	shell:
		"wget -O {output} {URL_ref_gff}"


rule reference_unzip:  # unpacks reference
	input:
		fna="reference/Saccharomyces_cerevisiae.fna.gz",
		gff="reference/Saccharomyces_cerevisiae.gff.gz"
	output:
		fna="reference/Saccharomyces_cerevisiae.fna",
		gff="reference/Saccharomyces_cerevisiae.gff"
	shell:
		"""
		gunzip -c {input.fna} > {output.fna} && gunzip -c {input.gff} > {output.gff}
		"""


rule fastqc:  # runs FastQC for all samples (8 files)
	input:
		"raw_data/{sample}.fastq.gz"
	output:
		"results/fastqc/{sample}_fastqc.html",
		"results/fastqc/{sample}_fastqc.zip"
	shell:
	    "fastqc -o results/fastqc/ --noextract {input}"


rule multiqc:  # runs MultiQC for all samples
	input:
		"results/fastqc"
	output:
		"results/multiqc/multiqc_report.html"
	shell:
		"""
		multiqc {input} -o results/multiqc/
		"""


rule build_index:  # builds genome indices for hisat2
    input:
        "reference/Saccharomyces_cerevisiae.fna"
    output:
        expand("results/hisat/index/genome_index.{index}.ht2", index=range(1, 9))
    params:
        index_prefix="genome_index"
    shell:
        "hisat2-build {input} results/hisat/index/{params.index_prefix}"


rule run_hisat_se:  # executes aligning in single-end mode with hisat2 
    input:
        fastq=rules.sample_download.output,
        index=rules.build_index.output
    output:
        "results/BAM/{sample}.bam"
    params:
        index_prefix="genome_index",
        threads=4
    shell:
        "hisat2 -p {params.threads} -x results/hisat/index/{params.index_prefix} -U {input.fastq} | samtools sort > {output}"


#rule gff2gtf:  # converts file with .gff format to .gtf
#    input:
#        "reference/Saccharomyces_cerevisiae.gff"
#    output:
#        "reference/Saccharomyces_cerevisiae.gtf"
#    shell:
#        "gffread {input} -T -o {output}"


rule featurecounts:  # calculates the count of features using featureCounts tool
    input:
        bams=expand("results/BAM/{sample}.bam", sample=["y00_1", "y00_2", "y30_1", "y30_2"]),
        ref="reference/Saccharomyces_cerevisiae.gff"
    output:
        tsv="results/count/all_samples.tsv",
        summary="results/count/all_samples.tsv.summary"
    log:
        "logs/feature_all_samples.log"
    shell:
        "featureCounts -t gene -g ID -a {input.ref} -o {output.tsv} {input.bams} 2> {log}"
        
        
rule simplify_the_counts:  # does subset of the data
    input:
        "results/count/all_samples.tsv"
    output:
        "results/count/simple_counts.txt"
    shell:
        "cat {input} | cut -f 1,7-10 > {output}"


rule run_deseq2:  # calculates metrics using r.script
	input:
		"results/count/simple_counts.txt"
	output:
		"results/deseq2/result.txt",
		"results/deseq2/norm-matrix-deseq2.txt"
	shell:
		"cat {input} | R -f r_scripts/deseq2.r"


rule draw_heatmap:  # draws heatmap using r.script
	input:
		"results/deseq2/norm-matrix-deseq2.txt"
	output:
		"results/heatmaps/output.pdf"
	shell:
		"cat {input} | R -f r_scripts/draw-heatmap.r"


rule cut_50_genes:  # cuts first 50 genes from result file
	input:
		"results/deseq2/result.txt"
	output:
		"results/deseq2/genes.txt"
	shell:
		"""
		head -n 50 {input} | cut -f 1 | cut -d "-" -f 2 > {output}
		"""
