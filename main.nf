params.graffite_vcf  = false
params.vcf           = false
params.RM_vcf        = false // mainly for debug. Requires --RM_dir
params.RM_dir        = false // mainly for debug. Requires --RM_vcf
params.genotype      = true
params.graph_method  = "pangenie" // or giraffe or graphaligner
params.reads         = "reads.csv"
params.assemblies    = "assemblies.csv"
params.reference     = "reference.fa"
params.TE_library    = "TE_library.fa"
params.out           = "out"
params.tsd_win       = 30 // add default value for TSD window search
params.cores         = false // set to an integer
params.mammal        = false
params.mini_K        = "500M"
params.stSort_m      = "4G"
params.stSort_t      = 4
params.version       = "0.2 beta (11-11-2022)"

// ideally, we should have defaults relative to genome size
params.svim_asm_memory     = null
params.repeatmasker_memory = null
params.pangenie_memory     = null
params.make_graph_memory   = null
params.graph_align_memory  = null
params.vg_call_memory      = null
params.min_mapq            = 0
params.min_support         = "2,4"

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
    graph_threads        = params.cores
} else {
    repeatmasker_threads = params.repeatmasker_threads
    svim_asm_threads     = params.svim_asm_threads
    pangenie_threads     = params.pangenie_threads
    graph_threads        = params.graph_threads
}

// initiate channels that will provide the reference genome to processes
Channel.fromPath(params.reference).into{ref_geno_ch; ref_asm_ch; ref_repeatmasker_ch; ref_tsd_ch; ref_tsd_search_ch}

// if the user provide alternative assemblies, align and call SV
// initiate channels 
if(!params.graffite_vcf && !params.vcf && !params.RM_vcf) {
    Channel.fromPath(params.assemblies).splitCsv(header:true).map{row ->
        [row.sample, file(row.path, checkIfExists:true)]}.set{asm_ch}
  asm_ch.combine(ref_asm_ch).set{svim_in_ch}
  Channel.fromPath(params.TE_library).set{TE_library_ch}
// This process perform alternative assemblies alignments and SV calling for each file in --assemblies
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

}
// if the user doesn't provide a VCF already made by GraffiTE with --graffite_vcf, use RepeatMasker to annotate repeats
if(!params.graffite_vcf) {
  // in case the RepeatMasker input is an already available sequence resolved VCF  
  if(params.vcf){
    Channel.fromPath(params.vcf).set{raw_vcf_ch}
    Channel.fromPath(params.TE_library).set{TE_library_ch}

  process repeatmasker_fromVCF {
    cpus repeatmasker_threads
    memory params.repeatmasker_memory
    publishDir "${params.out}/2_Repeat_Filtering", mode: 'copy'

    input:
    file("genotypes.vcf") from raw_vcf_ch
    file(TE_library) from TE_library_ch
    file(ref_fasta) from ref_repeatmasker_ch

    output:
    file("genotypes_repmasked_filtered.vcf") into tsd_ch, tsd_search_ch, tsd_gather_ch
    path("repeatmasker_dir/") into tsd_RM_ch, tsd_search_RM_ch

    script:
    if(params.mammal)
    """
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
    // otherwise, merge individuals' VCFs from the svim_asm process and annotate with RepeatMasker
  } else if(!params.RM_vcf){

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
    path("repeatmasker_dir/") into tsd_RM_ch, tsd_search_RM_ch

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

  }

  // if --RM_vcf is given, starts here and set the input channel
  if(params.RM_vcf){
    Channel.fromPath(params.RM_vcf).into{tsd_ch; tsd_search_ch; tsd_gather_ch}
    Channel.fromPath(params.RM_dir).into{tsd_RM_ch; tsd_search_RM_ch}
    }
  // this first TSD process stage the list of SV to scan
  process tsd_prep {

    input:
    file("genotypes_repmasked_filtered.vcf") from tsd_ch
    path("repeatmasker_dir/*") from tsd_RM_ch
    file(ref_fasta) from ref_tsd_ch

    output:
    file("indels.txt") into tsd_search_input
    file("SV_sequences_L_R_trimmed_WIN.fa") into tsd_search_SV
    file("flanking_sequences.fasta") into tsd_search_flanking

    script:
    """
    cp repeatmasker_dir/repeatmasker_dir/* .
    prepTSD.sh ${ref_fasta} ${params.tsd_win}
    """
  }
  // this second process actually search for TSDs
  process tsd_search {

    input:
    file indels from tsd_search_input.splitText( by: 2 )
    file("genotypes_repmasked_filtered.vcf") from tsd_search_ch.toList()
    file("SV_sequences_L_R_trimmed_WIN.fa") from tsd_search_SV.toList()
    file("flanking_sequences.fasta") from tsd_search_flanking.toList()
    path("repeatmasker_dir/*") from tsd_search_RM_ch.toList()
    file(ref_fasta) from ref_tsd_search_ch.toList()

    output:
    path('*TSD_summary.txt') into tsd_out_ch
    path('*TSD_full_log.txt') into tsd_full_out_ch

    script:
    """
    cp repeatmasker_dir/repeatmasker_dir/* .
    TSD_Match_v2.sh SV_sequences_L_R_trimmed_WIN.fa flanking_sequences.fasta ${indels}
    """
  }
  // this process combine the TSD search into single output files
  process tsd_report {
    publishDir "${params.out}/3_TSD_search", mode: 'copy'

    input:
    path(x) from tsd_out_ch.collect()
    path(y) from tsd_full_out_ch.collect()
    path("genotypes_repmasked_filtered.vcf") from tsd_gather_ch

    output:
    path("TSD_summary.txt") into tsd_sum_group_ch
    path("TSD_full_log.txt") into tsd_full_group_ch
    path("pangenome.vcf") into vcf_ch, vcf_merge_ch

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
    bcftools annotate -a \${TSD_FILE}.gz -h \${HDR_FILE} -c CHROM,POS,~ID,INFO/TSD genotypes_repmasked_filtered.vcf | bcftools view > pangenome.vcf
    """
  }

} else {
  // if a vcf is provided as parameter, skip discovery and go directly to genotyping
    Channel.fromPath(params.graffite_vcf).into{vcf_ch; vcf_merge_ch}
}

if(params.genotype) {
    Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:true)]}.set{reads_ch}
    if(params.graph_method == "pangenie") {
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
            PanGenie -t ${pangenie_threads} -j ${pangenie_threads} -s ${sample_name} -i <(zcat -f ${sample_reads}) -r ${ref} -v ${vcf} -o ${sample_name}
            bgzip ${sample_name}_genotyping.vcf
            tabix -p vcf ${sample_name}_genotyping.vcf.gz
            """
        }
    }

    else if(params.graph_method != "pangenie") {
        // will be useful to know the graph filenames later
        String graph =  ""
        switch(params.graph_method) {
            case "giraffe":
                graph = "index.giraffe.gbz"
                break
            case "graphaligner":
                graph = "index.vg"
                break
        }

        process makeGraph {
            cpus graph_threads
            memory params.make_graph_memory
            input:
            file vcf from vcf_ch
            file fasta from ref_geno_ch

            output:
            file "index" into graph_index_ch, vg_index_call_ch

            script:
            prep = """
            bcftools sort -Oz -o sorted.vcf.gz ${vcf}
            tabix sorted.vcf.gz
            mkdir index
            """
            finish = """
            vg snarls index/${graph} > index/index.pb
            """
            switch(params.graph_method) {
                case "giraffe":
                    prep + """
                    vg autoindex --tmp-dir \$PWD  -p index/index -w giraffe -v sorted.vcf.gz -r ${fasta}
                    """ + finish
                    break
                case "graphaligner":
                    prep + """
                    export TMPDIR=$PWD
                    vg construct -a  -r ${fasta} -v ${vcf} -m 1024 > index/index.vg
                    """ + finish
                    break
            }
        }

        reads_ch.combine(graph_index_ch).set{reads_align_ch}
        process graphAlignReads {
            cpus graph_threads
            memory params.graph_align_memory
            input:
            set val(sample_name), file(sample_reads), file("index") from reads_align_ch

            output:
            set val(sample_name), file("${sample_name}.gam"), file("${sample_name}.pack") into aligned_ch

            script:
            pack =  """
            vg pack -x index/${graph} -g ${sample_name}.gam -o ${sample_name}.pack -Q ${params.min_mapq}
            """

            switch(params.graph_method) {
                case "giraffe":
                    """
                    vg giraffe -t ${graph_threads} -Z index/index.giraffe.gbz -m index/index.min -d index/index.dist -i -f ${sample_reads} > ${sample_name}.gam
                    """ + pack
                    break
                case "graphaligner":
                    """
                    GraphAligner -t ${graph_threads} -x vg -g index/index.vg -f ${sample_reads} -a ${sample_name}.gam
                    """ + pack
                    break
            }
        }

        aligned_ch.combine(vg_index_call_ch).set{graph_pack_ch}
        process vgCall {
            cpus graph_threads
            memory params.vg_call_memory

            input:
            set val(sample_name), file(gam), file(pack), file("index") from graph_pack_ch

            output:
            file("${sample_name}.vcf.gz*") into indexed_vcfs

            script:
            """
            vg call -a -m ${params.min_support} -r index/index.pb -s ${sample_name} -k ${pack} index/${graph} > ${sample_name}.vcf
            bgzip ${sample_name}.vcf
            tabix ${sample_name}.vcf.gz
            """
        }
    }

process mergeVcfs {
  publishDir "${params.out}/4_Genotyping", mode: 'copy', glob: 'GraffiTE.merged.genotypes.vcf'

  input:
  file vcfFiles from indexed_vcfs.collect()
  path pangenome_vcf from vcf_merge_ch

  output:
  file "GraffiTE.merged.genotypes.vcf" into typeref_outputs

  script:
  """
  ls *vcf.gz > vcf.list
  bcftools merge -l vcf.list > GraffiTE.merged.genotypes.vcf
  bgzip GraffiTE.merged.genotypes.vcf
  tabix -p vcf GraffiTE.merged.genotypes.vcf.gz
  grep '#' ${pangenome_vcf} > P_header
  grep -v '#' ${pangenome_vcf} | sort -k1,1 -k2,2n > P_sorted_body
  cat P_header P_sorted_body > pangenome.sorted.vcf
  bgzip pangenome.sorted.vcf
  tabix -p vcf pangenome.sorted.vcf.gz
  bcftools annotate -a pangenome.sorted.vcf.gz -c CHROM,POS,ID,INFO GraffiTE.merged.genotypes.vcf.gz > GraffiTE.merged.genotypes.vcf
  """
  }
}
