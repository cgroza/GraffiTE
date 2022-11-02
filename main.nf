params.vcf        = false
params.genotype   = true
params.reads      = "reads.csv"
params.assemblies = "assemblies.csv"
params.reference  = "reference.fa"
params.TE_library = "TE_library.fa"
params.out        = "out"
params.tsd_win    = 30 // add default value for TSD window search
params.cores      = false // set to an integer
params.mammal     = false
params.mini_K     = "4G"
params.stSort_m   = "4G"
params.stSort_t   = 4
params.version    = "0.0 beta (2022)"


// SAY HELLO

log.info """

  ▄████  ██▀███   ▄▄▄        █████▒ █████▒██▓▄▄▄█████▓▓█████
 ██▒ ▀█▒▓██ ▒ ██▒▒████▄    ▓██   ▒▓██           ██▒ ▓▒▓█   ▀
▒██░▄▄▄░▓██ ░▄█ ▒▒██  ▀█▄  ▒████ ░▒████ ░▒██▒▒ ▓██░ ▒░▒███
░▓█  ██▓▒██▀▀█▄  ░██▄▄▄▄██ ░▓█▒  ░░▓█▒  ░░██░░ ▓██▓ ░ ▒▓█  ▄
░▒▓███▀▒░██▓ ▒██▒ ▓█   ▓██▒░▒█░   ░▒█░   ░██░  ▒██▒ ░ ░▒████▒
 ░▒   ▒ ░ ▒▓ ░▒▓░ ▒▒   ▓▒█░ ▒ ░    ▒ ░   ░▓    ▒ ░░   ░░ ▒░ ░
  ░   ░   ░▒ ░ ▒░  ▒   ▒▒ ░ ░      ░      ▒ ░    ░     ░ ░  ░
░ ░   ░   ░░   ░   ░   ▒    ░ ░    ░ ░    ▒ ░  ░         ░
      ░    ░           ░  ░               ░              ░  ░

                  V . ${params.version}

Find and Genotype Transposable Elements Insertion Polymorphisms
      in Genome Assemblies using a Pangenomic Approach           

Authors: Cristian Groza and Clément Goubert
Bug/issues: https://github.com/cgroza/GraffiTE/issues

"""

// if user uses global preset for number of cores
if(params.cores) {
    repeatmasker_threads = params.cores
    svim_asm_threads     = params.cores
    pangenie_threads     = params.cores
} else {
    repeatmasker_threads = params.repeatmasker_threads
    svim_asm_threads     = params.svim_asm_threads
    pangenie_threads     = params.pangenie_threads
}

Channel.fromPath(params.reference).into{ref_geno_ch; ref_asm_ch; ref_repeatmasker_ch; ref_tsd_ch; ref_tsd_search_ch}

if(!params.vcf) {
    Channel.fromPath(params.assemblies).splitCsv(header:true).map{row ->
        [row.sample, file(row.path, checkIfExists:true)]}.set{asm_ch}
  asm_ch.combine(ref_asm_ch).set{svim_in_ch}
  Channel.fromPath(params.TE_library).set{TE_library_ch}

  process svim_asm {
    cpus svim_asm_threads
    memory params.svim_asm_memory
    publishDir "${params.out}/1_SV_search", mode: 'copy'

    input:
    set val(asm_name), file(asm), file(ref) from svim_in_ch

    output:
    set val(asm_name), file("${asm_name}.vcf") into svim_out_ch

    script:
    """
    mkdir asm
    minimap2 -a -x asm5 --cs -r2k -t ${svim_asm_threads} -K ${params.mini_K} ${ref} ${asm} | samtools sort -m${params.stSort_m} -@${params.stSort_t} -o asm/asm.sorted.bam -
    samtools index asm/asm.sorted.bam
    svim-asm haploid --min_sv_size 100 --types INS,DEL --sample ${asm_name} asm/ asm/asm.sorted.bam ${ref}
    sed 's/svim_asm\\./${asm_name}\\.svim_asm\\./g' asm/variants.vcf > ${asm_name}.vcf
    """
  }

  process repeatmasker {
    cpus repeatmasker_threads
    memory params.repeatmasker_memory
    publishDir "${params.out}/2_Repeat_Filtering", mode: 'copy'

    input:
    file(vcfs) from svim_out_ch.map{sample -> sample[1]}.collect()
    file(TE_library) from TE_library_ch
    file(ref_fasta) from ref_repeatmasker_ch

    output:
    file("genotypes_repmasked_filtered.vcf") into tsd_ch, tsd_search_ch, tsd_gather_ch
    path("repeatmasker_dir/") into tsd_RM_ch, tsd_search_RM_ch // my hope is to export the output within their folder

    script:
    if(params.mammal)
    """
    ls *.vcf > vcfs.txt
    SURVIVOR merge vcfs.txt 0.1 0 0 0 0 100 genotypes.vcf
    repmask_vcf.sh genotypes.vcf genotypes_repmasked.vcf.gz ${TE_library} MAM
    bcftools view -G genotypes_repmasked.vcf.gz | \
    awk -v FS='\t' -v OFS='\t' \
    '{if(\$0 ~ /#CHROM/) {\$9 = "FORMAT"; \$10 = "ref"; print \$0} else if(substr(\$0, 1, 1) == "#") {print \$0} else {\$9 = "GT"; \$10 = "1|0"; print \$0}}' | \
    awk 'NR==1{print; print "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"} NR!=1' | \
    bcftools view -i 'INFO/total_match_span > 0.80' -o genotypes_repmasked_temp.vcf
    fix_vcf.py --ref ${ref_fasta} --vcf_in genotypes_repmasked_temp.vcf --vcf_out genotypes_repmasked_filtered.vcf
    """
    else
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
    // publishDir "${params.out}", mode: 'copy'

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

    input:
    val indels from tsd_search_input.splitText() 
    file("genotypes_repmasked_filtered.vcf") from tsd_search_ch.toList()
    file("SV_sequences_L_R_trimmed_WIN.fa") from tsd_search_SV.toList()
    file("flanking_sequences.fasta") from tsd_search_flanking.toList()
    path("repeatmasker_dir/*") from tsd_search_RM_ch.toList() // my hope is that this will pull the files from the repeatmasker_dir to the working directory
    file(ref_fasta) from ref_tsd_search_ch.toList()

    output:
    path('*TSD_summary.txt') into tsd_out_ch
    path('*TSD_full_log.txt') into tsd_full_out_ch

    script:
    """
    cp repeatmasker_dir/repeatmasker_dir/* .
    findTSD.sh ${indels}
    name="\$(head -n 1 <(echo ${indels}))"
    mv TSD_summary.txt \${name}.TSD_summary.txt
    mv TSD_full_log.txt \${name}.TSD_full_log.txt
    """
  }

  process tsd_report {
    publishDir "${params.out}/3_TSD_search", mode: 'copy'

    input:
    path(x) from tsd_out_ch.collect()
    path(y) from tsd_full_out_ch.collect()
    path("genotypes_repmasked_filtered.vcf") from tsd_gather_ch

    output:
    path("TSD_summary.txt") into tsd_sum_group_ch
    path("TSD_full_log.txt") into tsd_full_group_ch
    path("pangenie.vcf") into vcf_ch

    script:
    """
    cat ${x} > TSD_summary.txt
    cat ${y} > TSD_full_log.txt
    join -13 -21 <(grep -v "#" genotypes_repmasked_filtered.vcf | cut -f 1-3 | sort -k3,3) <(grep 'PASS' TSD_summary.txt | awk '{print \$1"\t"\$(NF-2)","\$(NF-1)}' | sort -k1,1) | \
    awk '{print \$2"\t"\$3"\t"\$1"\t"\$4}' | \
    sort -k1,1 -k2,2n > TSD_annotation
    HDR_FILE=\$(mktemp)
    echo -e '##INFO=<ID=TSD,Number=1,Type=String,Description="Target site duplication sequence passing filters">' >> \${HDR_FILE}
    TSD_FILE=TSD_annotation 
    bgzip \${TSD_FILE}
    tabix -s1 -b2 -e2 \${TSD_FILE}.gz 
    bcftools annotate -a \${TSD_FILE}.gz -h \${HDR_FILE} -c CHROM,POS,ID,INFO/TSD genotypes_repmasked_filtered.vcf | bcftools view > pangenie.vcf
    """
  }

} else {
  // if a vcf is provided as parameter, skip discovery and go directly to genotyping
  Channel.fromPath(params.vcf).set{vcf_ch}
}

if(params.genotype) {

  Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:true)]}.set{reads_ch}
  reads_ch.combine(vcf_ch).combine(ref_geno_ch).set{input_ch}

  process pangenie {
    cpus pangenie_threads
    memory params.pangenie_memory
    publishDir "${params.out}/4_Genotyping", mode: 'copy'

    input:
    set val(sample_name), file(sample_reads), file(vcf), file(ref) from input_ch

    output:
    file("${sample_name}_genotyping.vcf.gz*") into indexed_vcfs

    script:
    """
    PanGenie -t ${pangenie_threads} -j ${pangenie_threads} -s ${sample_name} -i ${sample_reads} -r ${ref} -v ${vcf} -o ${sample_name}
    bgzip ${sample_name}_genotyping.vcf
    tabix -p vcf ${sample_name}_genotyping.vcf.gz
    """
  }

  process mergeVcfs {
  publishDir "${params.out}/4_Genotyping", mode: 'copy'

  input:
  file vcfFiles from indexed_vcfs.collect()
   
  output:
  file "GraffiTE.merged.genotypes.vcf" into typeref_outputs
  
  script:
  """
  ls *vcf.gz > vcf.list
  bcftools merge -l vcf.list > GraffiTE.merged.genotypes.vcf
  """
  }
}
