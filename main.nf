// SAY HELLO

log.info """

▄████  ██▀███   ▄▄▄        █████▒ █████▒██▓▄▄▄█████▓▓█████
██▒ ▀█▒▓██ ▒ ██▒▒████▄    ▓██   ▒▓██           ██▒ ▓▒▓█   ▀
▒██░▄▄▄░▓██ ░▄█ ▒▒██  ▀█▄  ▒████ ░▒████ ░▒██▒▒ ▓██░ ▒░▒███
░▓█  ██▓▒██▀▀█▄  ░██▄▄▄▄██ ░▓█▒  ░░▓█▒  ░░██░░ ▓██▓ ░ ▒▓█  ▄
░▒▓███▀▒░██▓ ▒██▒  █   ▓██▒░▒█░   ░▒█░   ░██░  ▒██▒ ░ ░▒████▒
░▒   ▒ ░ ▒▓ ░▒▓░ ▒▒   ▓▒█░ ▒ ░    ▒ ░   ░▓    ▒ ░░   ░░ ▒░ ░
░   ░   ░▒ ░ ▒░  ▒   ▒▒ ░ ░      ░      ▒ ░    ░     ░ ░  ░
░ ░   ░   ░░   ░   ░   ▒    ░ ░    ░ ░    ▒ ░  ░         ░
░    ░           ░  ░               ░              ░  ░

V . ${workflow.commitId}

Find and Genotype Transposable Elements Insertion Polymorphisms
in Genome Assemblies using a Pangenomic Approach

Authors: Cristian Groza and Clément Goubert
Bug/issues: https://github.com/cgroza/GraffiTE/issues

"""

// if user uses global preset for number of cores

if (params.cores) {
  map_longreads_threads = params.cores
  map_asm_threads      = params.cores
  repeatmasker_threads = params.cores
  svim_asm_threads     = params.cores
  pangenie_threads     = params.cores
  graph_align_threads  = params.cores
  vg_call_threads      = params.cores
  sniffles_threads     = params.cores
}
else {
  map_longreads_threads = params.map_longreads_threads
  map_asm_threads       = params.map_asm_threads
  repeatmasker_threads  = params.repeatmasker_threads
  svim_asm_threads      = params.svim_asm_threads
  pangenie_threads      = params.pangenie_threads
  graph_align_threads   = params.graph_align_threads
  vg_call_threads       = params.vg_call_threads
  sniffles_threads      = params.sniffles_threads
}

String graph =  ""
switch(params.graph_method) {
  case "giraffe":
    graph = "index.giraffe.gbz"
    break
  case "graphaligner":
    graph = "index.vg"
    break
}

process break_scaffold {
  cpus 1

  input:
  tuple val(asm_name), path(asm)

  output:
  tuple val(asm_name), path("broken/${asm_base_name}.fa.gz")

  script:
  asm_base_name = asm.getName()
  """
  mkdir broken
  breakgaps.py ${asm} | gzip > broken/${asm_base_name}.fa.gz
  """
}

process map_asm {
  cpus map_asm_threads
  memory params.map_asm_memory
  time params.map_asm_time

  input:
  tuple val(asm_name), path(asm), path(ref)

  output:
  tuple val(asm_name), path("asm.sorted.bam"), path(ref), emit: map_asm_ch

  script:
  if(params.aligner == "minimap2") {
    """
    minimap2 -a -x ${params.asm_divergence} --cs -r2k -t ${map_asm_threads} -K ${params.mini_K} ${ref} ${asm} | samtools sort -m${params.stSort_m} -@${params.stSort_t} -o asm.sorted.bam -
      """
  }
  else if(params.aligner == "winnowmap") {
    """
    meryl count k=19 output merylDB ${ref}
    meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt
    winnowmap -a  -x ${params.asm_divergence} --cs -r2k -t ${map_asm_threads} -K ${params.mini_K} -W repetitive_k19.txt ${ref} ${asm} | samtools sort -m${params.stSort_m} -@${params.stSort_t} -o asm.sorted.bam -
    """
  }
}

process map_longreads {
  cpus map_longreads_threads
  memory params.map_longreads_memory
  time params.map_longreads_time

  input:
  tuple val(sample_name), path(longreads), val(type), path(ref)

  output:
  tuple val(sample_name), path("${sample_name}.bam"), path(ref), emit: map_longreads_ch

  script:
  read_preset = "map-${type}"
  if(type == "lr:hq") {
    read_preset = "${type}"
  }

  if(params.aligner == "minimap2") {
    """
    minimap2 -t ${map_longreads_threads} -ax ${read_preset} ${ref} ${longreads} | samtools sort -m${params.stSort_m} -@${params.stSort_t} -o ${sample_name}.bam  -
      """
  }
  else if(params.aligner == "winnowmap") {
    """
    meryl count k=15 output merylDB ${ref}
    meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

    winnowmap -W repetitive_k15.txt -t ${map_longreads_threads} -ax ${read_preset} ${ref} ${longreads} | samtools sort -m${params.stSort_m} -@${params.stSort_t} -o ${sample_name}.bam  -
    """
  }

}

process sniffles_sample_call {
  cpus sniffles_threads
  memory params.sniffles_memory
  time params.sniffles_time
  publishDir "${params.out}/1_SV_search/sniffles2_individual_VCFs", mode: 'copy'

  input:
  tuple val(sample_name), path(longreads_bam), path(ref)

  output:
  path("${sample_name}.snf"), emit: sniffles_sample_call_out_ch
  path("${sample_name}.vcf")

  script:
  """
  samtools index ${longreads_bam}
  sniffles --minsvlen 100 --threads ${sniffles_threads} --reference ${ref} --input ${longreads_bam} --snf ${sample_name}.snf --vcf ${sample_name}.vcf
  """
}

process sniffles_population_call {
  cpus sniffles_threads
  memory params.sniffles_memory
  publishDir "${params.out}/1_SV_search", mode: 'copy'

  input:
  path(snfs)
  path(ref)

  output:
  path("*variants.vcf"), emit: sn_variants_ch
  path("snfs.tsv")

  """
  ls *.snf > snfs.tsv
  sniffles --minsvlen 100  --threads ${sniffles_threads} --reference ${ref} --input snfs.tsv --vcf genotypes_unfiltered.vcf
  bcftools filter -i 'INFO/SVTYPE == "INS" | INFO/SVTYPE == "DEL"' genotypes_unfiltered.vcf | awk '\$5 != "<INS>" && \$5 != "<DEL>"' > sniffles2_variants.vcf
  """
}

process svim_asm {
  cpus svim_asm_threads
  memory params.svim_asm_memory
  time params.svim_asm_time
  publishDir "${params.out}/1_SV_search/svim-asm_individual_VCFs/", mode: 'copy'

  input:
  tuple val(asm_name), path(asm_bam), path(ref)

  output:
  tuple val(asm_name), path("${asm_name}.vcf"), emit: svim_out_ch

  script:
  """
  mkdir asm
  samtools index ${asm_bam}
  svim-asm haploid --min_sv_size 100 --types INS,DEL --sample ${asm_name} asm/ ${asm_bam} ${ref}
  sed 's/svim_asm\\./${asm_name}\\.svim_asm\\./g' asm/variants.vcf > ${asm_name}.vcf
  """
}

process survivor_merge {
  cpus svim_asm_threads
  memory params.svim_asm_memory
  publishDir "${params.out}/1_SV_search", mode: 'copy'

  input:
  path(vcfs)

  output:
  path("svim-asm_variants.vcf"), emit: sv_variants_ch
  path("vcfs.txt")

  script:
  """
  ls *.vcf > vcfs.txt
  SURVIVOR merge vcfs.txt 0.1 0 1 0 0 100 svim-asm_variants.vcf
  """
}

process merge_svim_sniffles2 {
  cpus params.merge_svim_sniffles2_threads
  memory params.merge_svim_sniffles2_memory
  time params.merge_svim_sniffles2_time

  publishDir "${params.out}/1_SV_search", mode: 'copy'

  input:
  path(svim_vcf)
  path(sniffles_vcf)

  output:
  path("svim-sniffles_merged_variants.vcf"), emit: sv_sn_variants_ch

  script:
  """
  ls sniffles2_variants.vcf svim-asm_variants.vcf > svim-sniffles2.vcfs.txt
  SURVIVOR merge svim-sniffles2.vcfs.txt 0.1 0 1 0 0 100 svim-sniffles2_merge_genotypes.vcf

  # header part to keep
  HEADERTOP=\$(grep '#' svim-sniffles2_merge_genotypes.vcf | grep -v 'CHROM')
  # modify last header line to fit content
  HEADERLINE=\$(grep '#CHROM' svim-sniffles2_merge_genotypes.vcf | awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tFORMAT\tGT"}')
  # add new info fields
  HEADERMORE=header_more
  echo -e '##INFO=<ID=sniffles2_SUPP,Number=1,Type=String,Description="Support vector from sniffle2-population calls">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=sniffles2_SVLEN,Number=1,Type=Integer,Description="SV length as called by sniffles2-population">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=sniffles2_SVTYPE,Number=1,Type=String,Description="Type of SV from sniffle2-population calls">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=sniffles2_ID,Number=1,Type=String,Description="ID from sniffle2-population calls">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=svim-asm_SUPP,Number=1,Type=String,Description="Support vector from svim-asm calls">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=svim-asm_SVLEN,Number=1,Type=Integer,Description="SV length as called by svim-asm">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=svim-asm_SVTYPE,Number=1,Type=String,Description="Type of SV from svim-asm calls">' >> \${HEADERMORE}
  echo -e '##INFO=<ID=svim-asm_ID,Number=1,Type=String,Description="ID from svim-asm calls">' >> \${HEADERMORE}
  # arrange the body part
  BODY=body_file
  paste -d ";" <(grep -v '#' svim-sniffles2_merge_genotypes.vcf | \
    cut -f 1-8) <(grep -v '#' svim-sniffles2_merge_genotypes.vcf | \
    cut -f 10 | sed 's/:/\t/g' | \
    awk '{print "sniffles2_SUPP="\$2";sniffles2_SVLEN="\$3";sniffles2_SVTYPE="\$7";sniffles2_ID="\$8}') <(grep -v '#' svim-sniffles2_merge_genotypes.vcf | \
    cut -f 11 | sed 's/:/\t/g' | awk '{print "svim-asm_SUPP="\$2";svim-asm_SVLEN="\$3";svim-asm_SVTYPE="\$7";svim-asm_ID="\$8"\t\\.\t\\."}') >> \${BODY}
  # concatenate and save "variants" file
  cat <(echo "\${HEADERTOP}") \${HEADERMORE} <(echo "\${HEADERLINE}") \${BODY} > svim-sniffles_merged_variants.vcf
  """

}

process repeatmask_VCF {
  cpus repeatmasker_threads
  memory params.repeatmasker_memory
  time params.repeatmasker_time
  publishDir "${params.out}/2_Repeat_Filtering", mode: 'copy'

  input:
  path("genotypes.vcf")
  path(TE_library)
  path(ref_fasta)

  output:
  path("genotypes_repmasked_filtered.vcf"), emit: RM_vcf_ch
  path("repeatmasker_dir/"), emit: RM_dir_ch
  //path("vcf_annotation.gz")
  //path("vcf_annotation_1")

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
  memory params.tsd_memory
  time params.tsd_time

  input:
  path("genotypes_repmasked_filtered.vcf")
  path("repeatmasker_dir/*")
  path(ref_fasta)

  output:
  path("indels.txt"), emit: tsd_search_input

  path("SV_sequences_L_R_trimmed_WIN.fa"), emit: tsd_search_SV
  path("flanking_sequences.fasta"), emit: tsd_search_flanking

  script:
  """
  cp repeatmasker_dir/repeatmasker_dir/* .
  prepTSD.sh ${ref_fasta} ${params.tsd_win}
  #wc -l indels.txt > indel_len
  """
}

process tsd_search {
  memory params.tsd_memory
  time params.tsd_time

  input:
  file indels
  path("genotypes_repmasked_filtered.vcf")
  path("SV_sequences_L_R_trimmed_WIN.fa")
  path("flanking_sequences.fasta")
  path("repeatmasker_dir/*")
  path(ref_fasta)

  output:
  path('*TSD_summary.txt'), emit: tsd_out_ch
  path('*TSD_full_log.txt'), emit: tsd_full_out_ch

  script:
  """
  cp repeatmasker_dir/repeatmasker_dir/* .
  TSD_Match_v2.sh SV_sequences_L_R_trimmed_WIN.fa flanking_sequences.fasta ${indels}
  """
}

process tsd_report {
  memory params.tsd_memory
  time params.tsd_time
  publishDir "${params.out}/3_TSD_search", mode: 'copy'

  input:
  path(x)
  path(y)
  path("genotypes_repmasked_filtered.vcf")

  output:
  path("TSD_summary.txt"), emit: tsd_sum_group_ch
  path("TSD_full_log.txt"), emit: tsd_full_group_ch
  path("pangenome.vcf"), emit: vcf_ch

  script:
  """
  cat ${x} > TSD_summary.txt
  cat ${y} > TSD_full_log.txt
  join -13 -21 <(grep -v "#" genotypes_repmasked_filtered.vcf | cut -f 1-5 | sort -k3,3) <(grep 'PASS' TSD_summary.txt | awk '{print \$1"\t"\$(NF-2)","\$(NF-1)}' | sort -k1,1) | \
    awk '{print \$2"\t"\$3"\t"\$1"\t"\$4"\t"\$5"\t"\$6}' | \
    sort -k1,1 -k2,2n > TSD_annotation
  HDR_FILE=header_file
  echo -e '##INFO=<ID=TSD,Number=1,Type=String,Description="Target site duplication sequence passing filters">' >> \${HDR_FILE}
  TSD_FILE=TSD_annotation
  bgzip \${TSD_FILE}
  tabix -s1 -b2 -e2 \${TSD_FILE}.gz
  bcftools annotate -a \${TSD_FILE}.gz -h \${HDR_FILE} -c CHROM,POS,~ID,REF,ALT,INFO/TSD genotypes_repmasked_filtered.vcf | bcftools view > pangenome.vcf
  """
}

process pangenie {
  cpus pangenie_threads
  memory params.pangenie_memory
  time params.pangenie_time
  publishDir "${params.out}/4_Genotyping", mode: 'copy'

  input:
  tuple val(sample_name), path(sample_reads), val(preset), path(vcf), path(ref)

  output:
  path("${sample_name}_genotyping.vcf.gz*"), emit: indexed_vcfs

  script:
  """
  PanGenie -t ${pangenie_threads} -j ${pangenie_threads} -s ${sample_name} -i <(zcat -f ${sample_reads}) -r ${ref} -v ${vcf} -o ${sample_name}
  bgzip ${sample_name}_genotyping.vcf
  tabix -p vcf ${sample_name}_genotyping.vcf.gz
  """
}

process make_graph {
  cpus params.make_graph_threads
  memory params.make_graph_memory
  time params.make_graph_time
  publishDir "${params.out}/GraffiTE_graph/", mode: 'copy'
  input:
  path(vcf)
  path(fasta)

  output:
  path("index"), emit: graph_index_ch

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

process graph_align_reads {
  cpus graph_align_threads
  memory params.graph_align_memory
  time params.graph_align_time
  errorStrategy 'finish'

  input:
  tuple val(sample_name), path(sample_reads), val(preset), path("index")

  output:
  tuple val(sample_name), path("${sample_name}.gam"), path("${sample_name}.pack"), emit: aligned_ch

  script:
  pack =  """
  vg pack -x index/${graph} -g ${sample_name}.gam -o ${sample_name}.pack -Q ${params.min_mapq}
  """

  interleaved = "-i"
  if(preset != "default") {
    interleaved = ""
  }

  switch(params.graph_method) {
    case "giraffe":
      """
      vg giraffe --parameter-preset ${preset} -t ${graph_align_threads} --index-basename index/index ${interleaved} -f ${sample_reads} > ${sample_name}.gam
      """ + pack
      break
    case "graphaligner":
      """
      GraphAligner -t ${graph_align_threads} -x vg -g index/index.vg -f ${sample_reads} -a ${sample_name}.gam
      """ + pack
      break
  }
}

process vg_call {
  cpus vg_call_threads
  memory params.vg_call_memory
  time params.vg_call_time

  input:
  tuple val(sample_name), path(gam), path(pack), path("index")

  output:
  path("${sample_name}.vcf.gz*"), emit: indexed_vcfs

  script:
  """
  vg call -a -m ${params.min_support} -r index/index.pb -s ${sample_name} -k ${pack} index/${graph} > ${sample_name}.vcf
  bgzip ${sample_name}.vcf
  tabix ${sample_name}.vcf.gz
  """
}

process merge_VCFs {
  memory params.merge_vcf_memory
  time params.merge_vcf_time
  publishDir "${params.out}/4_Genotyping", mode: 'copy', glob: 'GraffiTE.merged.genotypes.vcf'

  input:
  path(vcfFiles)
  path(pangenome_vcf)

  output:
  path("GraffiTE.merged.genotypes.vcf.gz"), emit: typeref_outputs

  script:
  """
  find . -name "*vcf.gz" | sort > vcf.list
  bcftools merge -l vcf.list > GraffiTE.merged.genotypes.vcf
  bgzip GraffiTE.merged.genotypes.vcf
  tabix -p vcf GraffiTE.merged.genotypes.vcf.gz
  grep '#' ${pangenome_vcf} > P_header
  grep -v '#' ${pangenome_vcf} | sort -k1,1 -k2,2n > P_sorted_body
  cat P_header P_sorted_body > pangenome.sorted.vcf
  bgzip pangenome.sorted.vcf
  tabix -p vcf pangenome.sorted.vcf.gz
  bcftools annotate -a pangenome.sorted.vcf.gz -c CHROM,POS,ID,REF,ALT,INFO GraffiTE.merged.genotypes.vcf.gz > GraffiTE.merged.genotypes.vcf
  rm -f GraffiTE.merged.genotypes.vcf.gz
  bgzip GraffiTE.merged.genotypes.vcf
  """
}

workflow {
  // initiate channels that will provide the reference genome to processes
  Channel.fromPath(params.reference, checkIfExists:true).set{ref_asm_ch}

  if(!params.graffite_vcf && !params.vcf && !params.RM_vcf) {
    if(params.longreads || params.bams) {
      sniffles_reads_in_ch = channel.of()
      sniffles_bams_in_ch = channel.of()

      if(params.longreads) {
        Channel.fromPath(params.longreads).splitCsv(header:true).map{row ->
          [row.sample, file(row.path, checkIfExists:true), row.type]}.combine(ref_asm_ch).set{map_longreads_in_ch}
        map_longreads(map_longreads_in_ch)
        sniffles_reads_in_ch = map_longreads.out.map_longreads_ch
      }

      if(params.bams) {
        sniffles_bams_in_ch = Channel.fromPath(params.bams).splitCsv(header:true).map{row ->
          [row.sample, file(row.path, checkIfExists:true)]}.combine(ref_asm_ch)
      }

      sniffles_sample_call(sniffles_reads_in_ch.concat(sniffles_bams_in_ch))
      sniffles_population_call(sniffles_sample_call.out.sniffles_sample_call_out_ch.collect(), ref_asm_ch)
    }

    if(params.assemblies) {
      Channel.fromPath(params.assemblies).splitCsv(header:true).map{row ->
        [row.sample, file(row.path, checkIfExists:true)]}.set{map_asm_in_ch}
      if(params.break_scaffolds) {
        map_asm_in_ch = break_scaffold(map_asm_in_ch)
      }
      map_asm(map_asm_in_ch.combine(ref_asm_ch))
      svim_asm(map_asm.out.map_asm_ch)
      survivor_merge(svim_asm.out.svim_out_ch.map{sample -> sample[1]}.collect())
    }

    if(params.assemblies && (params.longreads || params.bams)) {
      merge_svim_sniffles2(survivor_merge.out.sv_variants_ch, sniffles_population_call.out.sn_variants_ch)
    }
  }

  // if the user doesn't provide a VCF already made by GraffiTE with --graffite_vcf, use RepeatMasker to annotate repeats
  if(!params.graffite_vcf) {
    // except if --RM_vcf is given, in which case skip RepeatMasker here and set the input channel
    if(params.RM_vcf){
      Channel.fromPath(params.RM_vcf, checkIfExists:true).set{RM_vcf_ch}
      Channel.fromPath(params.RM_dir, checkIfExists:true).set{RM_dir_ch}
    } else {
      Channel.fromPath(params.TE_library, checkIfExists:true).set{TE_library_ch}
      // we need to set the vcf input depending what was given
      if((params.longreads || params.bams) && !params.assemblies){
        sniffles_population_call.out.sn_variants_ch.set { raw_vcf_ch }
      } else if (params.assemblies && !(params.longreads || params.bams)){
        survivor_merge.out.sv_variants_ch.set { raw_vcf_ch }
      } else if (params.assemblies && (params.longreads || params.bams)){
        merge_svim_sniffles2.out.sv_sn_variants_ch.set { raw_vcf_ch }
      } else if(params.vcf){
        Channel.fromPath(params.vcf, checkIfExists : true).set{raw_vcf_ch}
      } else {
        error "No --longreads, --assemblies, --vcf or --RM_vcf parameters passed to GraffiTE."
      }
      repeatmask_VCF(raw_vcf_ch, TE_library_ch, ref_asm_ch)
      repeatmask_VCF.out.RM_vcf_ch.set{RM_vcf_ch}
      repeatmask_VCF.out.RM_dir_ch.set{RM_dir_ch}
    }
    tsd_prep(RM_vcf_ch, RM_dir_ch, ref_asm_ch)
    tsd_search(tsd_prep.out.tsd_search_input.splitText( by: params.tsd_batch_size),
               RM_vcf_ch.toList(),
               tsd_prep.out.tsd_search_SV.toList(),
               tsd_prep.out.tsd_search_flanking.toList(),
               RM_dir_ch.toList(),
               ref_asm_ch.toList())
    tsd_report(tsd_search.out.tsd_out_ch.collect(),
               tsd_search.out.tsd_full_out_ch.collect(),
               RM_vcf_ch)
    tsd_report.out.vcf_ch.set{vcf_ch}
  } else {
    // if a vcf is provided as parameter, skip discovery and go directly to genotyping
    Channel.fromPath(params.graffite_vcf).set{vcf_ch}
  }

  if(params.genotype) {
    Channel.fromPath(params.reads).splitCsv(header:true).map{ row ->
      def parameter_preset = null
      switch(row.type) {
        case "pb":
          parameter_preset = "hifi"
          break
        case "hifi":
          parameter_preset = "hifi"
          break
        case "ont":
          parameter_preset = "r10"
          break
        default:
          parameter_preset = "default"
          break
      }
      [row.sample, file(row.path, checkIfExists:true), parameter_preset]
    }.set{reads_ch}

    if(params.graph_method == "pangenie") {
      reads_ch.combine(vcf_ch).combine(ref_asm_ch).set{input_ch}
      pangenie(input_ch)
      pangenie.out.indexed_vcfs.set{indexed_vcfs}
    }

    else if(params.graph_method == "giraffe" || params.graph_method == "graphaligner") {
      make_graph(vcf_ch, ref_asm_ch)
      reads_ch.combine(make_graph.out.graph_index_ch).set{reads_align_ch}
      graph_align_reads(reads_align_ch)
      graph_align_reads.out.aligned_ch.combine(make_graph.out.graph_index_ch).set{graph_pack_ch}
      vg_call(graph_pack_ch)
      vg_call.out.indexed_vcfs.set{indexed_vcfs}
    } else {
      error "Unsupported --graph_method. --graph_method must be pangenie, giraffe or graphaligner."
    }

    merge_VCFs(indexed_vcfs.collect(), vcf_ch)
  }
}
