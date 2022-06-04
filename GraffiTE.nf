params.vcf = false
params.reads = "reads.csv"
params.assemblies = "assemblies.csv"
params.reference = "reference.fa"
params.TE_library = "TE_library.fa"
params.threads = 1
params.memory = "100G"
params.out = "out"

Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path)]}.set{reads_ch}
Channel.fromPath(params.reference).into{ref_geno_ch; ref_asm_ch; ref_make_vcf_ch}

if(!params.vcf) {
  Channel.fromPath(params.assemblies).splitCsv(header:true).map{row -> [row.sample, file(row.path)]}.set{asm_ch}
  asm_ch.combine(ref_asm_ch).set{svim_in_ch}
  Channel.fromPath(params.TE_library).set{TE_library_ch}

  process svim_asm {
    cpus params.threads
    memory params.memory
    publishDir "${params.out}", mode: 'copy'

    input:
    set val(asm_name), file(asm), file(ref) from svim_in_ch

    output:
    set val(asm_name), file("${asm_name}.vcf") into svim_out_ch

    script:
    """
    mkdir asm
    minimap2 -a -x asm5 --cs -r2k -t ${params.threads} ${ref} ${asm} | samtools sort -m4G -@4 -o asm/asm.sorted.bam -
    samtools index asm/asm.sorted.bam
    svim-asm haploid --min_sv_size 100 --types INS,DEL --sample ${asm_name} asm/ asm/asm.sorted.bam ${ref}
    sed 's/svim_asm\\./${asm_name}\\.svim_asm\\./g' asm/variants.vcf > ${asm_name}.vcf
    """
  }

  process make_vcf {
    cpus params.threads
    memory params.memory
    publishDir "${params.out}", mode: 'copy'

    input:
    file(vcfs) from svim_out_ch.map{sample -> sample[1]}.collect()
    file(TE_library) from TE_library_ch
    file(ref_fasta) from ref_make_vcf_ch

    output:
    file("pangenie.vcf") into vcf_ch

    script:
    """
    ls *.vcf > vcfs.txt
    SURVIVOR merge vcfs.txt 0.1 0 0 0 0 100 genotypes.vcf
    repmask_vcf.sh genotypes.vcf genotypes_repmasked.vcf.gz ${TE_library}
    bcftools view -G genotypes_repmasked.vcf.gz | \
    awk -v FS='\t' -v OFS='\t' \
    '{if(\$0 ~ /#CHROM/) {\$9 = "FORMAT"; \$10 = "ref"; print \$0} else if(substr(\$0, 1, 1) == "#") {print \$0} else {\$9 = "GT"; \$10 = "1|0"; print \$0}}' | \
    bcftools view -i 'INFO/match_span > 0.80'  -o pangenie_temp.vcf
    fix_deletions.py --ref ${ref_fasta} --vcf_in pangenie_temp.vcf --vcf_out pangenie.vcf
    """

  }
} else {
  Channel.fromPath(params.vcf).set{vcf_ch}
}


reads_ch.combine(vcf_ch).combine(ref_geno_ch).set{input_ch}

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
