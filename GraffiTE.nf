params.vcf = "variants.vcf"
params.reads = "reads.csv"
params.assemblies = "assemblies.csv"
params.reference = "reference.fa"
params.threads = 1
params.memory = "100G"
params.out = "out"

Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path)]}.view().set{reads_ch}
Channel.fromPath(params.assemblies).splitCsv(header:true).map{row -> [row.sample, file(row.path)]}.view().set{asm_ch}

Channel.fromPath(params.vcf).view().set{vcf_ch}
Channel.fromPath(params.reference).view().into{ref_geno_ch; ref_asm_ch}

reads_ch.combine(vcf_ch).combine(ref_geno_ch).view().set{input_ch}
asm_ch.combine(ref_asm_ch).view().set{svim_in_ch}

process svim_asm {
    cpus params.threads
    memory params.memory
    publishDir "${params.out}", mode: 'copy'

    input:
    set val(asm_name), file(asm), file(ref) from svim_in_ch

    output:
    set val(asm_name), file("${asm_name}.vcf.gz") into svim_out_ch

    script:
    """
    mkdir asm
    minimap2 -a -x asm5 --cs -r2k -t ${params.threads} ${ref} ${asm} | samtools sort -m4G -@4 -o asm/asm.sorted.bam -
    samtools index asm/asm.sorted.bam
    svim-asm haploid --min_sv_size 100 --types INS --sample ${asm_name} asm/ asm/asm.sorted.bam ${ref}
    mv asm/variants.vcf ${asm_name}.vcf
    bgzip ${asm_name}.vcf
    """
}

process pangenie {
	  cpus params.threads
	  memory params.memory
	  publishDir "${params.out}", mode: 'copy'

    input:
	  set val(sample_name), file(sample_reads), file(vcf), file(ref) from input_ch

	  output:
	  file("${sample_name}_genotyping.vcf")

    script:
    """
    PanGenie -t ${params.threads} -j ${params.threads} -s ${sample_name} -i ${sample_reads} -r ${ref} -v ${vcf}  -o ${sample_name}
    """
}
