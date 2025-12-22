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
  input:
  tuple val(asm_name), path(asm), path(ref)

  output:
  tuple val(asm_name), path("asm.sorted.bam"), path(ref)

  script:
  if(params.aligner == "minimap2") {
    """
    minimap2 -a -x ${params.asm_divergence} --cs -r2k -t ${task.cpus} -K ${params.mini_K} ${ref} ${asm} | \
      samtools sort -m${params.stSort_m} -@${params.stSort_t} -o asm.sorted.bam -
    """
  }
  else if(params.aligner == "winnowmap") {
    """
    meryl count k=19 output merylDB ${ref}
    meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt
    winnowmap -a  -x ${params.asm_divergence} --cs -r2k -t ${task.cpus} -K ${params.mini_K} -W repetitive_k19.txt ${ref} ${asm} | \
      samtools sort -m${params.stSort_m} -@${params.stSort_t} -o asm.sorted.bam -
    """
  }
}

process map_longreads {
  input:
  tuple val(sample_name), path(longreads), val(type), path(ref)

  output:
  tuple val(sample_name), path("${sample_name}.bam"), path(ref)

  script:
  read_preset = "map-${type}"
  if(type == "lr:hq") {
    read_preset = "${type}"
  }

  if(params.aligner == "minimap2") {
    """
    minimap2 -t ${task.cpus} -ax ${read_preset} ${ref} ${longreads} | \
      samtools sort -m${params.stSort_m} -@${params.stSort_t} -o ${sample_name}.bam  -
      """
  }
  else if(params.aligner == "winnowmap") {
    """
    meryl count k=15 output merylDB ${ref}
    meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

    winnowmap -W repetitive_k15.txt -t ${task.cpus} -ax ${read_preset} ${ref} ${longreads} | \
      samtools sort -m${params.stSort_m} -@${params.stSort_t} -o ${sample_name}.bam  -
    """
  }

}

process sniffles_sample_call {
  input:
  tuple val(sample_name), path(longreads_bam), path(ref)

  output:
  tuple path("${sample_name}.snf"), path("${sample_name}.vcf")

  script:
  """
  samtools index ${longreads_bam}
  sniffles --minsvlen 100 --threads ${task.cpus} --reference ${ref} --input ${longreads_bam} --snf ${sample_name}.snf --vcf ${sample_name}.vcf
  """
}

process sniffles_population_call {
  publishDir "${params.out}/1_SV_search", mode: 'copy'

  input:
  path(snfs)
  path(ref)

  output:
  path("sniffles2_individual_VCFs/*.vcf.gz")

  """
  ls *.snf > snfs.tsv
  sniffles --minsvlen 100  --threads ${task.cpus} --reference ${ref} --input snfs.tsv --vcf genotypes_unfiltered.vcf
  bcftools filter -i 'INFO/SVTYPE == "INS" | INFO/SVTYPE == "DEL"' genotypes_unfiltered.vcf | awk '\$5 !~ "<INS>" && \$5 !~ "<DEL>"' | \
    bcftools sort -Oz -o sniffles2_variants.vcf.gz
  mkdir sniffles2_individual_VCFs
  bcftools +split sniffles2_variants.vcf.gz -Oz -o  sniffles2_individual_VCFs
  """
}

process svim_asm {
  publishDir "${params.out}/1_SV_search/svim-asm_individual_VCFs/", mode: 'copy'

  input:
  tuple val(asm_name), path(asm_bam), path(ref)

  output:
  tuple val(asm_name), path("${asm_name}.vcf.gz")

  script:
  """
  mkdir asm
  samtools index ${asm_bam}
  svim-asm haploid --min_sv_size 100 --types INS,DEL --sample ${asm_name} asm/ ${asm_bam} ${ref}
  sed 's/svim_asm\\./${asm_name}\\.svim_asm\\./g' asm/variants.vcf | bcftools sort -Oz -o ${asm_name}.vcf.gz
  """
}

process truvari_merge {
  publishDir "${params.out}/1_SV_search", mode: 'copy'

  input:
  path(vcfs)
  path(ref)

  output:
  path("SVs.vcf")

  script:
  """
  for f in ${vcfs}
  do
  tabix \${f}
  done

  bcftools merge -Oz -m none -o merged.vcf.gz *.vcf.gz
  tabix merged.vcf.gz
  truvari collapse --chain -P 0.5 -p 0.5 -S -1 -k common -i merged.vcf.gz -o truvari_merged.vcf
  bcftools +setGT truvari_merged.vcf -- -t . -n 0 | bcftools norm -f ${ref} > truvari_merged_filled.vcf
  shorten_ids.py --vcf_in  truvari_merged_filled.vcf --vcf_out SVs.vcf
  """
}


process split_repeatmask {
  input:
  path(vcf)

  output:
  path("*.vcf")

  script:
  """
  bcftools sort -Oz -o ${vcf}.gz ${vcf}
  tabix ${vcf}.gz
  bcftools index -s ${vcf}.gz | cut -f 1 | while read C; do bcftools view -O v -o \${C}.vcf ${vcf}.gz "\${C}" ; done
  """
}

process concat_repeatmask {
  publishDir "${params.out}/3_TSD_search", mode: 'copy'
  input:
  path("tsd_pangenome_*.vcf")
  path("TSD_full_log_*.txt")
  path("TSD_summary_*.txt")
  path(ref_fasta)

  output:
  path("pangenome.vcf"), emit: vcf_ch
  path("TSD_summary.txt")
  path("TSD_full_log.txt")

  script:
  """
  cat TSD_summary_*.txt > TSD_summary.txt
  cat TSD_full_log_*.txt > TSD_full_log.txt
  bcftools concat tsd_pangenome_*.vcf | \
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "LC_ALL=C sort -k1,1 -k2,2n"}' | \
    bcftools view -Ov -o pangenome_temp.vcf -i 'INFO/total_match_span > 0.80'
  fix_vcf.py --ref ${ref_fasta} --vcf_in pangenome_temp.vcf --vcf_out pangenome.vcf
  """
}

process repeatmask_VCF {
  publishDir "${params.out}/2_Repeat_Filtering/${task.index}", mode: 'copy'

  input:
  tuple path("genotypes.vcf"), path(TE_library), path(ref_fasta)

  output:
  tuple path("genotypes_repmasked_filtered.vcf"), path("repeatmasker_dir/")

  script:
  def mammal = ""
  if(params.mammal) {
    mammal = "MAM"
  }
  """
  repmask_vcf.sh genotypes.vcf genotypes_repmasked.vcf.gz ${TE_library} ${mammal}
  bcftools view -Ov -o genotypes_repmasked_filtered.vcf -i 'INFO/total_match_span > 0.80' genotypes_repmasked.vcf.gz
  """
}

process tsd_prep {
  input:
  tuple path("genotypes_repmasked_filtered.vcf"), path("repeatmasker_dir/*"), path(ref_fasta)

  output:
  tuple path("genotypes_repmasked_filtered.vcf"), path("repeatmasker_dir/repeatmasker_dir"), path(ref_fasta),
    path("indels.txt"), path("SV_sequences_L_R_trimmed_WIN.fa"), path("flanking_sequences.fasta")

  script:
  """
  cp repeatmasker_dir/repeatmasker_dir/* .
  prepTSD.sh ${ref_fasta} ${params.tsd_win}
  """
}

process tsd_search {
  input:
  tuple path("genotypes_repmasked_filtered.vcf"), path("repeatmasker_dir/*"), path(ref_fasta), file(indels),
    path("SV_sequences_L_R_trimmed_WIN.fa"), path("flanking_sequences.fasta")

  output:
  tuple path('*TSD_summary.txt'), path('*TSD_full_log.txt'), path("genotypes_repmasked_filtered.vcf"), env("CHROM")

  script:
  """
  CHROM=\$(bcftools view -H genotypes_repmasked_filtered.vcf | cut -f1 | uniq)
  cp repeatmasker_dir/repeatmasker_dir/* .
  TSD_Match_v2.sh SV_sequences_L_R_trimmed_WIN.fa flanking_sequences.fasta ${indels}
  """
}

process tsd_report {
  input:
  tuple path(x), path(y), path("genotypes_repmasked_filtered.vcf"), val(chrom)

  output:
  path("TSD_summary.txt"), emit: tsd_sum_group_ch
  path("TSD_full_log.txt"), emit: tsd_full_group_ch
  path("pangenome.vcf"), emit: vcf_ch

  script:
  """
  cat ${x} > TSD_summary.txt
  cat ${y} > TSD_full_log.txt
  join -13 -21 <(grep -v "#" genotypes_repmasked_filtered.vcf | cut -f 1-5 | \
    sort -k3,3) <(grep 'PASS' TSD_summary.txt | \
    awk '{print \$1"\t"\$(NF-2)","\$(NF-1)}' | sort -k1,1) | \
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

process pangenie_index {
  input:
  tuple path(vcf), path(ref)

  output:
  path("pangenie_index")

  script:
  """
  bcftools view -G ${vcf} | \
    awk -v FS='\t' -v OFS='\t' \
    '{if(\$0 ~ /#CHROM/) {\$9 = "FORMAT"; \$10 = "ref"; print \$0} else if(substr(\$0, 1, 1) == "#") {print \$0} else {\$9 = "GT"; \$10 = "1|0"; print \$0}}' | \
    awk 'NR==1{print; print "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"} NR!=1' | \
    bcftools norm -m+ -Ov -o graph.vcf
  merge_vcfs.py merge -r ${ref} -v graph.vcf -ploidy 2 > graph_merged.vcf
  mkdir pangenie_index
  PanGenie-index -v graph_merged.vcf -r ${ref} -t ${task.cpus} -o pangenie_index/pangenie_index
  """
}

process pangenie {
  publishDir "${params.out}/4_Genotyping", mode: 'copy'

  input:
  tuple val(sample_name), path(sample_reads), val(preset), path(index)
  path(ref)

  output:
  tuple val(sample_name), path("${sample_name}_genotyping.vcf.gz*")

  script:
  """
  PanGenie -t ${task.cpus} -j ${task.cpus} -s ${sample_name} -i <(zcat -f ${sample_reads}) -f ${index}/pangenie_index -o ${sample_name}
  bgzip ${sample_name}_genotyping.vcf
  tabix ${sample_name}_genotyping.vcf.gz
  bcftools norm -f ${ref} -m- -Oz -o ${sample_name}.vcf.gz ${sample_name}_genotyping.vcf.gz
  tabix -p vcf ${sample_name}.vcf.gz
  """
}

process make_graph {
  publishDir "${params.out}/GraffiTE_graph/", mode: 'copy'
  input:
  path(vcf)
  path(fasta)

  output:
  path("index")

  script:
  prep = """
  mkdir index
  """
  switch(params.graph_method) {
    case "giraffe":
      prep + """
      vg autoindex --tmp-dir \$PWD  -p index/index -w sr-giraffe -w lr-giraffe -v ${vcf} -r ${fasta}
      vg convert --vg-algorithm -f index/${graph} > index/index.gfa
      vg snarls index/${graph} > index/index.pb
      """
      break
    case "graphaligner":
      prep + """
      export TMPDIR=$PWD
      vg construct -a  -r ${fasta} -v ${vcf} -m 1024 > index/index.vg
      vg convert --vg-algorithm -f index/index.vg > index/index.gfa
      vg snarls index/index.gfa > index/index.pb
      """
      break
  }
}

process bam_to_fastq {
  input:
  tuple val(sample_name), path(sample_reads), val(preset)
  output:
  tuple val(sample_name), path("${sample_reads.baseName}.fq.gz"), val(preset)

  script:
  """
  samtools sort -n -@ ${task.cpus} ${sample_reads} | samtools fastq -@ ${task.cpus} - | pigz > ${sample_reads.baseName}.fq.gz
  """
}

process graph_align_reads {
  input:
  tuple val(sample_name), path(sample_reads), val(preset), path("index")

  output:
  tuple val(sample_name), path("${sample_name}.gaf.gz"), path("${sample_name}.pack")

  script:

  interleaved = "-i"
  if(preset != "default") {
    interleaved = ""
  }

  switch(params.graph_method) {
    case "giraffe":
      """
      vg giraffe --parameter-preset ${preset} -o gam -t ${graph_align_threads} --index-basename index/index ${interleaved} -f ${sample_reads} > ${sample_name}.gam
      vg pack -x index/${graph} -g ${sample_name}.gam -o ${sample_name}.pack -Q ${params.min_mapq}
      vg convert -G ${sample_name}.gam index/${graph} | subset_gaf.py | gzip > ${sample_name}.gaf.gz
      rm ${sample_name}.gam
      """
      break
    case "graphaligner":
      """
      GraphAligner -t ${task.cpus} -x vg -g index/index.gfa -f ${sample_reads} -a ${sample_name}.gam
      vg pack -x index/index.gfa -g ${sample_name}.gam -o ${sample_name}.pack -Q ${params.min_mapq}
      vg convert -G ${sample_name}.gam index/index.gfa | subset_gaf.py | gzip > ${sample_name}.gaf.gz
      rm ${sample_name}.gam
      """
      break
  }
}

process vg_call {
  input:
  tuple val(sample_name), path(gaf), path(pack), path("index"), path(ref)

  output:
  tuple val(sample_name), path("${sample_name}.vcf.gz*")

  script:
  """
  vg call -a -m ${params.min_support} -r index/index.pb -s ${sample_name} -k ${pack} index/${graph} | \
    bcftools norm -f ${ref} -m-  | \
    bcftools sort -Oz -o ${sample_name}.vcf.gz
  tabix ${sample_name}.vcf.gz
  """
}

process merge_VCFs {
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
