params.vcf = false
params.genotype = true
params.reads = "reads.csv"
params.assemblies = "assemblies.csv"
params.reference = "reference.fa"
params.TE_library = "TE_library.fa"
params.out = "out"
params.tsd_win = 30 // add default value for TSD window search

Channel.fromPath(params.reference).into{ref_geno_ch; ref_asm_ch; ref_repeatmasker_ch; ref_tsd_ch; ref_tsd_search_ch}

if(!params.vcf) {
    Channel.fromPath(params.assemblies).splitCsv(header:true).map{row ->
        [row.sample, file(row.path, checkIfExists:true)]}.set{asm_ch}
  asm_ch.combine(ref_asm_ch).set{svim_in_ch}
  Channel.fromPath(params.TE_library).set{TE_library_ch}

  process svim_asm {
    cpus params.svim_asm_threads
    memory params.svim_asm_memory
    publishDir "${params.out}", mode: 'copy'

    input:
    set val(asm_name), file(asm), file(ref) from svim_in_ch

    output:
    set val(asm_name), file("${asm_name}.vcf") into svim_out_ch

    script:
    """
    mkdir asm
    minimap2 -a -x asm5 --cs -r2k -t ${params.svim_asm_threads} ${ref} ${asm} | samtools sort -m4G -@4 -o asm/asm.sorted.bam -
    samtools index asm/asm.sorted.bam
    svim-asm haploid --min_sv_size 100 --types INS,DEL --sample ${asm_name} asm/ asm/asm.sorted.bam ${ref}
    sed 's/svim_asm\\./${asm_name}\\.svim_asm\\./g' asm/variants.vcf > ${asm_name}.vcf
    """
  }

  process repeatmasker {
    cpus params.repeatmasker_threads
    memory params.repeatmasker_memory
    publishDir "${params.out}", mode: 'copy'

    input:
    file(vcfs) from svim_out_ch.map{sample -> sample[1]}.collect()
    file(TE_library) from TE_library_ch
    file(ref_fasta) from ref_repeatmasker_ch

    output:
    file("genotypes_repmasked_filtered.vcf") into tsd_ch, tsd_search_ch
    path("repeatmasker_dir/") into tsd_RM_ch, tsd_search_RM_ch // my hope is to export the output within their folder

    script:
    """
    ls *.vcf > vcfs.txt
    SURVIVOR merge vcfs.txt 0.1 0 0 0 0 100 genotypes.vcf
    repmask_vcf.sh genotypes.vcf genotypes_repmasked.vcf.gz ${TE_library}
    bcftools view -G genotypes_repmasked.vcf.gz | \
    awk -v FS='\t' -v OFS='\t' \
    '{if(\$0 ~ /#CHROM/) {\$9 = "FORMAT"; \$10 = "ref"; print \$0} else if(substr(\$0, 1, 1) == "#") {print \$0} else {\$9 = "GT"; \$10 = "1|0"; print \$0}}' | \
    awk 'NR==1{print; print "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"} NR!=1' | \
    bcftools view -i 'INFO/total_match_span > 0.80' -o genotypes_repmasked_temp.vcf
    fix_vcf.py --ref ${ref_fasta} --vcf_in genotypes_repmasked_temp.vcf --vcf_out genotypes_repmasked_filtered.vcf
    """

  }

  process tsd_prep {
    // cpus params.tsd_search_threads
    // memory params.tsd_search_memory
    publishDir "${params.out}", mode: 'copy'

    input:
    file("genotypes_repmasked_filtered.vcf") from tsd_ch
    path("repeatmasker_dir/*") from tsd_RM_ch // my hope is that this will pull the files from the repeatmasker_dir to the working directory
    file(ref_fasta) from ref_tsd_ch

    output:
    file("indels.txt") into tsd_search_input
    file("SV_sequences_L_R_trimmed_WIN.fa") into tsd_search_SV
    file("flanking_sequences.fasta") into tsd_search_flanking
    //file("pangenie.vcf") into vcf_ch
    //path("TSD_*.txt") into tsd_out_ch

    script:
    """
    cp repeatmasker_dir/repeatmasker_dir/* .
    prepTSD.sh ${ref_fasta} ${params.tsd_win}
    """
  }
  
  process tsd_search {
    // cpus params.tsd_search_threads
    // memory params.tsd_search_memory
    //publishDir "${params.out}", mode: 'copy'

    input:
    //file("genotypes_repmasked_filtered.vcf") from tsd_ch
    file("indels.txt") from tsd_search_input.splitText() // check how to split by params.cpu
    file("genotypes_repmasked_filtered.vcf") from tsd_search_ch.toList()
    file("SV_sequences_L_R_trimmed_WIN.fa") from tsd_search_SV.toList()
    file("flanking_sequences.fasta") from tsd_search_flanking.toList()
    path("repeatmasker_dir/*") from tsd_search_RM_ch.toList() // my hope is that this will pull the files from the repeatmasker_dir to the working directory
    file(ref_fasta) from ref_tsd_search_ch.toList()

    output:
    path("TSD_summary.txt") into tsd_out_ch
    path("TSD_full_log.txt") into tsd_full_out_ch
    //file("pangenie.vcf") into vcf_ch

    script:
    """
    cp repeatmasker_dir/repeatmasker_dir/* .
    findTSD.sh
    """
  }

  process tsd_report {
    // cpus params.tsd_search_threads
    // memory params.tsd_search_memory
    //publishDir "${params.out}", mode: 'copy'

    input:
    file("TSD_summary.txt") from tsd_out_ch.collect()
    file("TSD_full_log.txt") from tsd_full_out_ch.collect()

    output:
    path("TSD_summary.txt") into tsd_sum_group_ch
    path("TSD_full_log.txt") into tsd_full_group_ch
    //file("pangenie.vcf") into vcf_ch

    script:
    """
    cat TSD_summary.txt* > TSD_summary.txt
    cat TSD_full_log.txt* > TSD_full_log.txt
    """
  }

} else {
  Channel.fromPath(params.vcf).set{vcf_ch}
}

if(params.genotype) {

  Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:true)]}.set{reads_ch}
  reads_ch.combine(vcf_ch).combine(ref_geno_ch).set{input_ch}

  process pangenie {
    cpus params.pangenie_threads
    memory params.pangenie_memory
    publishDir "${params.out}", mode: 'copy'

    input:
    set val(sample_name), file(sample_reads), file(vcf), file(ref) from input_ch

    output:
    file("${sample_name}_genotyping.vcf")

    script:
    """
    PanGenie -t ${params.pangenie_threads} -j ${params.pangenie_threads} -s ${sample_name} -i ${sample_reads} -r ${ref} -v ${vcf}  -o ${sample_name}
    """
  }
}
