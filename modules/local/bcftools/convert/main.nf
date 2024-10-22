process BCFTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(chromosomes)
    path(fasta_ref_37)

    output:
    tuple val(meta), path("${meta.id}.gq20_hs37d5.vcf.gz"), emit: converted_vcf
    tuple val(meta), path("${meta.id}.gq20_hs37d5.vcf.gz.tbi"), emit: converted_tbi

    script:
    def memory = task.memory.toMega()  // Use memory in megabytes
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view \\
        --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \\
        $vcf \\
	    --no-version \\
        --threads ${task.cpus} \\
    | bcftools annotate \\
        --rename-chrs $chromosomes \\
	    --no-version \\
    | bcftools convert \\
	    --gvcf2vcf \\
	    --fasta-ref $fasta_ref_37  \\
	    --no-version \\
    | bcftools norm \\
        --check-ref s \\
        --fasta-ref $fasta_ref_37 \\
        --multiallelics +any \\
	    --no-version \\
	    --threads ${task.cpus} \\
    | bcftools sort \\
        --max-mem ${memory}M \\
        -Oz \\
        -o ${prefix}.gq20_hs37d5.vcf.gz \\
    && bcftools index --tbi ${prefix}.gq20_hs37d5.vcf.gz --threads ${task.cpus}
    """
}