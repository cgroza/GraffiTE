params.vcf = "variants.vcf"
params.reads = "reads.csv"
params.reference = "reference.fa"
params.threads = 1
params.memory = "100G"
params.out = "out"

Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path)]}.view().set{reads_ch}

Channel.fromPath(params.vcf).view().set{vcf_ch}
Channel.fromPath(params.reference).view().set{ref_ch}

reads_ch.combine(vcf_ch).combine(ref_ch).view().set{input_ch}

process pangenie {
	  cpus params.threads
	  memory params.memory
	  publishDir "${params.out}", mode: 'copy'

    input:
	  set val(sample_name), file(sample_reads), file(vcf), file(ref) from input_ch

    output:
	  output:
	  file("${sample_name}_genotyping.vcf")

    script:
    """
    PanGenie -t ${params.threads} -j ${params.threads} -s ${sample_name} -i ${sample_reads} -r ${ref} -v ${vcf}  -o ${sample_name}
    """
}
