process BCFTOOLS_FILTER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), emit: filtered_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools norm \\
        --multiallelic -any \\
        $vcf \\
        --threads ${task.cpus} \\
	    --no-version \\
    | bcftools filter \\
        -i 'FORMAT/GQ>20 && FORMAT/DP>1.65' \\
        --set-GTs . \\
        --threads ${task.cpus} \\
	--no-version \\
    | bcftools view \\
        --trim-alt-alleles \\
        --exclude-uncalled \\
    | bcftools annotate \\
        -x ^FORMAT/GT,^FORMAT/GQ,^FORMAT/DP,^INFO/END \\
        --threads ${task.cpus} \\
	    --no-version \\
    | bcftools norm \\
        --multiallelics +any \\
        --threads ${task.cpus} \\
	    --no-version \\
        -Oz \\
        -o ${prefix}.filtered.vcf.gz
    """
}