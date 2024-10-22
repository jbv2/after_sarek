process LIFTOVER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf)
    tuple path(fasta_ref_19), path(fai)
    path(chain)

    output:
    tuple val(meta), path("${meta.id}.lifted.vcf.gz"), emit: lifted_vcf
    tuple val(meta), path("${meta.id}.tmp.vcf.gz"),    emit: tmp_vcf
    tuple val(meta), path("${meta.id}.tmp.vcf.gz.tbi"),    emit: tmp_tbi

    script:
    def memory = task.memory.toMega()  // Use memory in megabytes
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    CrossMap.py gvcf \\
	$chain \\
	$vcf \\
	$fasta_ref_19 \\
	${prefix}.lifted.vcf.gz \\
	--no-comp-alleles \\
	--compress \\
    && bcftools reheader \\
        --fai $fai \\
        --threads ${task.cpus} \\
        ${prefix}.lifted.vcf.gz \\
    | bcftools sort \\
        --max-mem ${memory}M \\
        -Oz \\
        -o ${prefix}.tmp.vcf.gz \\
    && bcftools index --tbi ${prefix}.tmp.vcf.gz --threads ${task.cpus}
    """
}