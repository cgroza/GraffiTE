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


process pav_asm {
  publishDir "${params.out}/1_SV_search/pav_individual_VCFs/", mode: 'copy'

  input:
  tuple val(sample_name), path(haps), path(ref)

  output:
  path("sv_${sample_name}.vcf.gz")

  script:
  """
  export XDG_CACHE_HOME=\$(pwd)
  echo "{\\"reference\\": \\"${ref}\\"}" > config.json

  printf 'NAME' > assemblies.tsv
  i=1
  for hap in ${haps}; do
    printf '\tHAP%s' "\${i}" >> assemblies.tsv
    ((i++))
  done
  printf '\n' >> assemblies.tsv

  printf '%s' "${sample_name}" >> assemblies.tsv
  for hap in ${haps}; do
    printf '\t%s' "\${hap}" >> assemblies.tsv
  done
  printf '\n' >> assemblies.tsv


  /opt/pav/files/docker/run -c ${task.cpus}
  bcftools filter -i 'ABS(INFO/SVLEN) > 50' -Oz -o sv_${sample_name}.vcf.gz ${sample_name}.vcf.gz
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
  val(from_vcf)

  output:
  path("SVs.vcf")

  script:
  """
  if [[ "${from_vcf}" == "true" ]]; then
    # User-supplied --vcf: a single VCF, no truvari collapse needed.
    # Pass IDs through unchanged (shorten_ids.py is only for collapse output).
    if [[ "${vcfs}" == *.gz ]]; then
      gunzip --force --stdout ${vcfs} > SVs.vcf
    else
      cp ${vcfs} SVs.vcf
    fi
  else

  for f in ${vcfs}
  do
  tabix \${f}
  done

  num_files=\$(ls -1q ${vcfs} | wc -l)

  if [[ "\$num_files" -eq "1" ]]; then
    gunzip --force --stdout ${vcfs} > SVs.vcf
  else

  for f in *.vcf.gz
  do
    bcftools annotate -x INFO \${f} -Oz -o stripped_\${f}
    tabix stripped_\${f}
    done

    bcftools merge -Oz -m none -o merged.vcf.gz stripped_*.vcf.gz
    tabix merged.vcf.gz

    mkdir -p collapsed
    truvari divide -T ${task.cpus} merged.vcf.gz shards

    for shard in shards/*.vcf.gz; do
        [[ -f "\${shard}.tbi" ]] || tabix "\${shard}"
    done

    printf '%s\n' shards/*.vcf.gz | \
    xargs -P "${task.cpus}" -I{} bash -c '
    shard="{}"
    base=\$(basename "\${shard}" .vcf.gz)
    truvari collapse \
      --chain -P 0.5 -p 0.5 -S -1 -k common \
      -i "\${shard}" \
      -o "collapsed/\${base}.vcf" && \
      bcftools sort \
      -o "collapsed/\${base}.sorted.vcf.gz" \
      -O z \
      "collapsed/\${base}.vcf" && \
      bcftools index --tbi "collapsed/\${base}.sorted.vcf.gz" && \
      rm "collapsed/\${base}.vcf"
    '

    bcftools concat -Oz -o unsorted.vcf.gz \
      \$(find collapsed/ -name '*.vcf.gz' | sort -V)
    bcftools sort  -Oz -o truvari_merged.vcf.gz unsorted.vcf.gz
    tabix truvari_merged.vcf.gz
    rm unsorted.vcf.gz

    bcftools +setGT truvari_merged.vcf.gz -- -t . -n 0 | bcftools norm -f ${ref} | \
    bcftools +fill-tags - -Ov -o truvari_merged_filled.vcf -- -t 'SVLEN=strlen(ALT)-strlen(REF)'
    shorten_ids.py --vcf_in  truvari_merged_filled.vcf --vcf_out SVs.vcf
  fi
  fi
  """
}



// HERV-K (HML-2) allele-state annotation and locus layer. --human only.
//
// Reads the raw RepeatMasker tables rather than INFO/repeat_ids: annotate_vcf.R
// collapses each RepeatMasker link group to one name plus "(x)", which erases
// the LTR-INT-LTR architecture this needs. It then masks a window of the
// reference at each candidate to establish what the REF allele actually holds.
//
// pangenome.vcf is READ ONLY here. It induces the graph, so it must stay
// byte-identical; only the human subset is annotated. Candidate calling and
// locus grouping still run over every candidate in pangenome.vcf, because the
// --human filter requires FILTER="PASS" and would otherwise be able to hide one
// member of a locus behind an unrelated caller flag.
process hervk_annotate {
  publishDir "${params.out}/3_TSD_search", mode: 'copy', overwrite: true

  input:
  path(pangenome_vcf, stageAs: 'in.pangenome.vcf')
  path(human_vcf, stageAs: 'in.pangenome.human.vcf')
  path(human_tsv, stageAs: 'in.pangenome.presence-absence_human.tsv')
  path("rmdir_*")
  path(ref_fasta)
  path(te_library)

  output:
  path("pangenome.human.vcf"), emit: human_vcf_ch
  path("pangenome.presence-absence_human.tsv")
  path("hervk_loci.tsv"), emit: loci_ch
  path("hervk_calls.tsv"), emit: calls_ch
  path("hervk_arch.tsv")
  path("hervk_refstate.tsv")
  path("hervk_candidates.vcf"), emit: hervk_candidates_ch
  path("hervk_polymorphism_summary.md")
  path("pangenome.human.consolidated.vcf"), emit: human_consolidated_ch
  path("hervk_discovery_consolidation_report.md")

  script:
  def cfg_arg = params.hervk_config ? "--config ${params.hervk_config}" : ""
  def strict_arg = params.hervk_strict ? "--strict" : ""
  """
  REF="${ref_fasta}"
  if [[ "\$REF" == *.gz ]]; then
      if ! (file -L "\$REF" | grep -q "BGZF"); then
          zcat "\$REF" | bgzip -c > ref.fa.gz
          REF=ref.fa.gz
      fi
  fi
  samtools faidx "\$REF"

  # Every LTR/ERVK record in the discovery VCF, not just the human subset.
  #
  # The |SVLEN| cap has to be applied HERE, at candidacy, not only in
  # hervk_classify.py's DEFAULTS (max_svlen=25000). hervk_candidate.ids is what
  # hervk_ref_state.py masks against, so a cap that only bites at classify time
  # drops the record from the calls *after* its reference window has already
  # been masked. On the CaG set that is one 25.3 Mb window
  # (chr1-120594342-DEL-25264467, a zero-query-footprint ALNTRUNC artefact)
  # carrying 95.8% of the 26.4 Mb sent to RepeatMasker -- 1h39m of masking plus
  # 2h14m of ProcessRepeats, which timed out two 4h jobs. With the cap here the
  # first pass is ~1.1 Mb, as RERUN_2.md predicts.
  #
  bcftools view -H -i 'matching_classes="LTR/ERVK" & abs(SVLEN)<=${params.hervk_max_svlen}' in.pangenome.vcf \\
    | cut -f3 > hervk_candidate.ids

  hervk_arch.py --rm-out rmdir_* --ids hervk_candidate.ids --out hervk_arch.tsv

  if [[ -n "${params.hervk_ref_annotation ?: ''}" ]]; then
    hervk_ref_state.py --vcf in.pangenome.vcf --ids hervk_candidate.ids \\
        --reference "\$REF" --rm-annotation ${params.hervk_ref_annotation} \\
        --flank ${params.hervk_ref_flank} --out hervk_refstate.tsv
  else
    hervk_ref_state.py --vcf in.pangenome.vcf --ids hervk_candidate.ids \\
        --reference "\$REF" --te-library ${te_library} \\
        --flank ${params.hervk_ref_flank} --threads ${task.cpus} \\
        --out hervk_refstate.tsv
  fi

  # Calls over the full candidate set, and a VCF of it. The --human pME
  # filter is narrower than the HERV-K candidate list on purpose -- it defines
  # the paper's TE set and is not ours to widen -- so a locus can lose
  # members to it. chr7:4.70 Mb loses two of three, including the one carrying
  # the common allele. This file is where those records keep their annotation
  # and their discovery genotypes; hervk_loci.tsv flags the split with
  # LOCUS_SPLIT_BY_HUMAN_FILTER.
  hervk_classify.py ${cfg_arg} --max-svlen ${params.hervk_max_svlen} \\
      --vcf-in in.pangenome.vcf \\
      --vcf-out hervk_candidates.vcf --vcf-out-candidates-only \\
      --calls-out hervk_calls.tsv \\
      --arch hervk_arch.tsv --ref-state hervk_refstate.tsv \\
      --summary hervk_polymorphism_summary.md

  hervk_classify.py ${cfg_arg} ${strict_arg} --max-svlen ${params.hervk_max_svlen} \\
      --vcf-in in.pangenome.human.vcf --vcf-out human.hervk.vcf \\
      --arch hervk_arch.tsv --ref-state hervk_refstate.tsv \\
      --tsv-in in.pangenome.presence-absence_human.tsv \\
      --tsv-out pangenome.presence-absence_human.tsv

  hervk_reconcile.py flag \\
      --calls hervk_calls.tsv \\
      --vcf-in human.hervk.vcf --vcf-out pangenome.human.vcf \\
      --loci-out hervk_loci.tsv --ref-state hervk_refstate.tsv \\
      --window ${params.hervk_locus_window}

  # The same loci merged into one record each, as a separate file.
  #
  # pangenome.human.vcf is what induces the graph, so its record structure is
  # fixed and a HERV-K locus reaches a reader there as a scatter of records.
  # This is that VCF with each locus collapsed onto one multi-allelic record
  # carrying its allele set, its counts from the assemblies, and whether it is
  # an insertion polymorphism. Genotyping still reads the unmerged file.
  #
  # Nothing is masked here. These genotypes are haplotype-resolved assembly
  # alignments, which resolve a tandem array directly; the graph consolidation
  # downstream is the one that has to withhold copy-number alleles.
  hervk_reconcile.py consolidate --source discovery \\
      --vcf-in    pangenome.human.vcf \\
      --loci      hervk_loci.tsv \\
      --calls     hervk_calls.tsv \\
      --reference "\$REF" \\
      --out-vcf   pangenome.human.consolidated.vcf \\
      --report    hervk_discovery_consolidation_report.md

  for v in pangenome.human.vcf pangenome.human.consolidated.vcf; do
    awk -v v="${params.graffite_version}" 'NR==1 && /^##fileformat/ {print; print "##GraffiTE_version="v; next} {print}' \\
        "\$v" > "\$v".tmp && mv "\$v".tmp "\$v"
  done
  """
}


// Stage E: the human merged genotypes VCF, with HERV-K loci consolidated.
// --human only, and only after graph genotyping.
//
// GraffiTE.merged.genotypes.vcf.gz (the full call set) is never rewritten; the
// human subset is a separate output and is where consolidation lands.
process hervk_reconcile {
  publishDir "${params.out}/4_Genotyping", mode: 'copy', overwrite: true

  input:
  path(merged_vcf)
  path(human_vcf)
  path(loci_tsv)
  path(calls_tsv)
  path(ref_fasta)
  val(genotyper)

  output:
  path("GraffiTE.merged.genotypes.human.vcf.gz"), emit: human_gt_ch
  path("GraffiTE.merged.genotypes.human.vcf.gz.tbi")
  path("hervk_unconsolidated_records.vcf")
  path("hervk_reconciliation_report.md")

  script:
  // hervk_mask_tandem is the old name for this switch; honour it while
  // anything is still passing it.
  def legacy   = params.hervk_mask_tandem
  def mask_cnv = (legacy == null) ? params.hervk_mask_graph_gt_at_cnv : legacy
  def mask_arg = mask_cnv ? "" : "--no-mask-cnv-gt"
  """
  # Subset the merged genotypes to the human candidate set, by ID. Records whose
  # ID did not survive merge_VCFs' `bcftools annotate` keep a raw snarl ID and
  # simply will not match -- the reconciler reports any locus member it cannot
  # locate rather than emitting a partial locus.
  bcftools query -f '%ID\\n' ${human_vcf} > human.ids
  bcftools view -i 'ID=@human.ids' -Ov -o merged.human.vcf ${merged_vcf}

  # Carry the HERV-K annotation across from the discovery VCF.
  #
  # Genotyping does not preserve it. merge_VCFs transfers INFO from
  # pangenome.vcf, which is deliberately left un-annotated because it induces
  # the graph, so the genotyped records arrive without HERVK_LOCUS or the locus
  # flags. In the CaG run that left 4 of 29 HERV-K records carrying a locus id
  # and none that could be filtered on HERVK_MEI, which is the field a user
  # needs to tell an insertion polymorphism from structural variation in an
  # element every haplotype carries.
  #
  # The tag list comes from the source header rather than being written out
  # here, so a new HERVK_* field travels without another edit.
  bcftools sort -Oz -o disc.annot.vcf.gz ${human_vcf}
  tabix -f -p vcf disc.annot.vcf.gz
  TAGS=\$(bcftools view -h disc.annot.vcf.gz \\
      | sed -n 's|^##INFO=<ID=\\(HERVK_[^,>]*\\).*|INFO/\\1|p' | paste -sd, -)
  if [[ -n "\$TAGS" ]]; then
    # bcftools annotate reads the target through htslib's indexed reader, so
    # the target has to be bgzipped and indexed even though it is only
    # streamed.
    bcftools sort -Oz -o merged.human.sorted.vcf.gz merged.human.vcf
    tabix -f -p vcf merged.human.sorted.vcf.gz
    bcftools annotate -a disc.annot.vcf.gz -c "\$TAGS" \\
        -Ov -o merged.human.vcf merged.human.sorted.vcf.gz
  fi

  hervk_reconcile.py consolidate ${mask_arg} \\
      --genotyped-vcf merged.human.vcf \\
      --loci          ${loci_tsv} \\
      --calls         ${calls_tsv} \\
      --discovery-vcf ${human_vcf} \\
      --reference     ${ref_fasta} \\
      --genotyper     ${genotyper} \\
      --out-vcf       GraffiTE.merged.genotypes.human.vcf \\
      --out-archive   hervk_unconsolidated_records.vcf \\
      --report        hervk_reconciliation_report.md

  bgzip -f GraffiTE.merged.genotypes.human.vcf
  tabix -p vcf GraffiTE.merged.genotypes.human.vcf.gz
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
  path("pangenome.trusted.vcf"), optional: true
  path("pangenome.human.vcf"), emit: human_vcf_ch, optional: true
  path("pangenome.presence-absence.tsv")
  path("pangenome.presence-absence_trusted.tsv"), optional: true
  path("pangenome.presence-absence_human.tsv"), emit: human_tsv_ch, optional: true
  path("human_filter_summary.txt"), optional: true
  path("TSD_summary.txt")
  path("TSD_full_log.txt")

  script:
  def trusted_filter = "n_hits==1 & abs(SVLEN)>=${params.trusted_min_svlen} & (ULTRA_TR_span<${params.trusted_max_ultra_span} | matching_classes=\"Simple_repeat\") & ((matching_classes!~\"LINE\" & matching_classes!~\"SINE\" & matching_classes!~\"Retroposon\") | polyA=\"TRUE\")"
  def trusted_filter_full = params.trusted_ignore_filter ? trusted_filter : "(${trusted_filter}) & FILTER=\"PASS\""

  // --human pME filter, applied directly to pangenome.vcf (not to the trusted
  // subset, which is a species-agnostic heuristic for non-model organisms).
  // Each whitelist param is a comma-separated list of regexes because bcftools
  // regexes have no alternation; ~ is matched element-wise on these Number=.
  // fields, and anchoring with ^ applies per element. Prefixes only: repeat_ids
  // carry "(x)" and "(VNTR_only)" suffixes from bin/annotate_vcf.R.
  def orIds = { csv -> '(' + csv.toString().split(',').collect{ "repeat_ids~\"${it.trim()}\"" }.join(' | ') + ')' }
  def grp   = { cls, csv -> csv?.toString()?.trim() ? "(matching_classes=\"${cls}\" & ${orIds(csv)})" : "matching_classes=\"${cls}\"" }
  def human_ids = [grp('SINE/Alu',       params.human_alu_ids),
                   grp('LINE/L1',        params.human_l1_ids),
                   grp('Retroposon/SVA', params.human_sva_ids),
                   grp('Simple_repeat',  params.human_sva_ids),
                   grp('LTR/ERVK',       params.human_hervk_ids)].join(' | ')
  def human_size   = "abs(SVLEN)>=${params.human_min_svlen} & (ULTRA_TR_span<${params.human_max_ultra_span} | matching_classes=\"Simple_repeat\")"
  // polyA (TPRT signature) is required for Alu/L1/SVA but not for HML-2 or for
  // SVA-VNTR expansions. Stated positively: bcftools "!~" does not negate
  // reliably on these Number=. fields (matching_classes!~"SINE" is true for
  // every record, including SINE/Alu ones), so the trusted filter's
  // "!~LINE & !~SINE & !~Retroposon" idiom must not be reused here.
  def human_single = "n_hits==1 & (matching_classes=\"LTR/ERVK\" | matching_classes=\"Simple_repeat\" | polyA=\"TRUE\")"
  // HML-2 proviral SVs where RepeatMasker splits a small SVA hit off the LTR
  // (SVA/LTR5_Hs homology). OR-ed at the n_hits level so it also bypasses the
  // polyA requirement that "Retroposon" in matching_classes would trigger.
  //
  // The hit count is a cap, not an equality, because RepeatMasker also splits
  // the internal region of a degraded or rearranged provirus. On the CaG set
  // two of the three records at chr7:4,699,714 come back with three hits
  // (LTR5_Hs,SVA_A,HERVK-int and HERVK-int,SVA_A,HERVK-int) where the third is
  // that split, not a second element. At n_hits==2 the filter kept one record
  // of that locus and dropped the two carrying its common allele, so the locus
  // reported 39/1 for a two-unit against three-unit difference where the truth
  // is 26/12/2 across three states.
  //
  // The cap does the work here, and the rest of the clause is what keeps it
  // honest: LTR/ERVK beside Retroposon/SVA, an ^HERVK-int id, and a length no
  // greater than one provirus. Relaxing n_hits on its own instead, anywhere in
  // human_single, admits 18 more records on this cohort that are an LTR5
  // fragment beside something else -- COMP-subunit_FAM90A, alpha satellite,
  // L1PA10, and HERVK9, a lineage human_hervk_ids deliberately excludes.
  def hervk_pair   = "n_hits<=${params.hervk_pair_max_hits} & matching_classes=\"LTR/ERVK\" & matching_classes=\"Retroposon/SVA\" & repeat_ids~\"^HERVK-int\" & abs(SVLEN)<=${params.hervk_pair_max_svlen}"
  def human_hits   = params.hervk_sva_pair ? "((${human_single}) | (${hervk_pair}))" : "(${human_single})"
  def human_filter_base = "(${human_ids}) & ${human_size} & ${human_hits}"
  def human_filter = params.human_ignore_filter ? human_filter_base : "(${human_filter_base}) & FILTER=\"PASS\""
  """
  cat TSD_summary_*.txt > TSD_summary.txt
  cat TSD_full_log_*.txt > TSD_full_log.txt
  bcftools concat tsd_pangenome_*.vcf | \
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "LC_ALL=C sort -k1,1 -k2,2n"}' | \
    bcftools view -Ov -o pangenome_temp.vcf -i 'INFO/total_repeat_span > ${params.repeat_span_cutoff}'
  # htslib can't index gzip-compressed fasta; re-compress with bgzip if needed
  REF="${ref_fasta}"
  if [[ "\$REF" == *.gz ]]; then
      if ! (file -L "\$REF" | grep -q "BGZF"); then
          zcat "\$REF" | bgzip -c > ref.fa.gz
          REF=ref.fa.gz
      fi
  fi
  fix_vcf.py --ref "\$REF" --vcf_in pangenome_temp.vcf --vcf_out pangenome_nopa.vcf
  add_polyA.py pangenome_nopa.vcf -o pangenome_raw.vcf

  # pangenome.vcf retains the original FILTER values from upstream.
  cp pangenome_raw.vcf pangenome.vcf

  # presence-absence TSV for the full callset
  vcf_to_pa_tsv.py pangenome.vcf -o pangenome.presence-absence.tsv

  if [[ "${params.human}" != "true" ]]; then
    # trusted subset: variants matching the trusted criteria. By default
    # also requires existing FILTER=="PASS"; bypass with --trusted_ignore_filter.
    bcftools view -Ov -o pangenome.trusted.vcf -i '${trusted_filter_full}' pangenome.vcf
    vcf_to_pa_tsv.py pangenome.trusted.vcf -o pangenome.presence-absence_trusted.tsv
  else
    # --human replaces the trusted subset with a polymorphic-MEI subset built
    # directly from pangenome.vcf: young subfamilies only (AluY*, L1HS,
    # SVA_D/E/F, HML-2), single RepeatMasker hit, plus the HERVK-int+SVA
    # proviral exception.
    bcftools view -Ov -o pangenome.human.vcf -i '${human_filter}' pangenome.vcf
    vcf_to_pa_tsv.py pangenome.human.vcf -o pangenome.presence-absence_human.tsv

    {
      echo "# GraffiTE --human pME filter"
      echo
      echo "filter expression:"
      echo '  ${human_filter}'
      echo
      printf 'records in pangenome.vcf       : %s\\n' "\$(bcftools view -H pangenome.vcf | wc -l | tr -d ' ')"
      printf 'records in pangenome.human.vcf : %s (HERV-K annotation is added downstream by hervk_annotate)\\n' "\$(bcftools view -H pangenome.human.vcf | wc -l | tr -d ' ')"
      echo
      echo "kept (count, matching_classes, repeat_ids):"
      bcftools query -f '%INFO/matching_classes\\t%INFO/repeat_ids\\n' pangenome.human.vcf | sort | uniq -c | sort -rn
      echo
      echo "dropped pME-class records (count, matching_classes, repeat_ids):"
      bcftools query -e '${human_filter}' -f '%INFO/matching_classes\\t%INFO/repeat_ids\\n' pangenome.vcf | \\
        awk -F'\\t' '\$1 ~ /Alu|L1|SVA|Simple_repeat|ERVK/' | sort | uniq -c | sort -rn
    } > human_filter_summary.txt

  fi

  # Stamp GraffiTE version into the header of each published VCF
  for VCF in pangenome.vcf pangenome.trusted.vcf pangenome.human.vcf; do
    [ -f "\$VCF" ] || continue
    awk -v v="${params.graffite_version}" 'NR==1 && /^##fileformat/ {print; print "##GraffiTE_version="v; next} {print}' "\$VCF" > "\$VCF.tmp" && mv "\$VCF.tmp" "\$VCF"
  done
  """
}

process repeatmask_VCF {
  publishDir "${params.out}/2_Repeat_Filtering/${task.index}", mode: 'copy'

  input:
  tuple path("genotypes.vcf"), path(TE_library), path(ref_fasta)

  output:
  tuple path("genotypes_repmasked_filtered.vcf"), path("repeatmasker_dir/"), emit: vcf
  path("ultra_out.bed"), emit: ultra_bed
  path("ultra_out.span"), emit: ultra_span
  path("genotypes_repmasked.vcf.gz"), emit: repmasked_vcf_debug
  path("vcf_annotation.bak.txt"), emit: vcf_annotation_debug
  path("union.bp"), emit: union_bp_debug
  path("total_repeat_span.tsv"), emit: total_repeat_span_debug
  path("combined.stats"), emit: combined_stats_debug
  path("ultra_out.stats"), emit: ultra_stats_debug

  script:
  def mammal = ""
  if(params.mammal) {
    mammal = "MAM"
  }
  """
  repmask_vcf.sh genotypes.vcf genotypes_repmasked.vcf.gz ${TE_library} ${mammal}
  bcftools view -Ov -o genotypes_repmasked_filtered.vcf -i 'INFO/total_repeat_span > ${params.repeat_span_cutoff}' genotypes_repmasked.vcf.gz
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
  prepTSD.sh ${ref_fasta} ${params.tsd_win} ${task.cpus}
  """
}

process tsd_search {
  input:
  tuple path("genotypes_repmasked_filtered.vcf"), path("repeatmasker_dir/*"), path(ref_fasta), path(indels),
    path("SV_sequences_L_R_trimmed_WIN.fa"), path("flanking_sequences.fasta")

  output:
  tuple path('*TSD_summary.txt'), path('*TSD_full_log.txt'), path("genotypes_repmasked_filtered.vcf"), path("chrom.txt")

  script:
  """
  bcftools view -H genotypes_repmasked_filtered.vcf | cut -f1 | uniq > chrom.txt
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
  tsd_annotate_vcf.sh genotypes_repmasked_filtered.vcf TSD_summary.txt pangenome.vcf
  """
}

// pangenie_graph_variants.tsv lists each ALT allele of pangenome.vcf with its
// ID in the graph and whether merge_vcfs.py kept it. PanGenie writes the graph
// ID to INFO/ID of the genotyped VCFs.
process pangenie_index {
  publishDir "${params.out}/4_Genotyping", mode: 'copy', pattern: 'pangenie_graph_variants.tsv'

  input:
  tuple path(vcf), path(ref)

  output:
  path("pangenie_index"), emit: index
  path("pangenie_graph_variants.tsv"), emit: table

  script:
  """
  bcftools view -G -Ov -o sites.vcf ${vcf}
  pangenie_graph_vcf.py prepare sites.vcf graph_input.vcf pangenie_graph_variants.tsv
  bcftools sort graph_input.vcf | bcftools norm -m+ -Ov -o graph.vcf
  merge_vcfs.py merge -r ${ref} -v graph.vcf -ploidy 2 > graph_merged.vcf
  pangenie_graph_vcf.py report pangenie_graph_variants.tsv graph_merged.vcf
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
  val(graph_method)

  output:
  path("index")

  script:
  prep = """
  mkdir index
  bcftools +setGT ${vcf} -- -t a -n u > unphased.vcf
  """
  switch(graph_method) {
    case "giraffe":
      prep + """
      vg autoindex --tmp-dir \$PWD  -p index/index -w sr-giraffe -w lr-giraffe -v unphased.vcf -r ${fasta}
      vg convert --vg-algorithm -f index/index.giraffe.gbz > index/index.gfa
      vg snarls index/index.giraffe.gbz > index/index.pb
      """
      break
    case "graphaligner":
      prep + """
      export TMPDIR=$PWD
      vg construct -a  -r ${fasta} -v unphased.vcf -m 1024 > index/index.vg
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

  samtools view -h ${sample_reads} \
    | awk 'BEGIN{OFS="\t"} /^@/{print; next} {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11}' \
    | samtools view -bS - \
    | samtools sort -n -@ ${task.cpus} - \
    | samtools fastq -@ ${task.cpus} - \
    | pigz > ${sample_reads.baseName}.fq.gz
  """
}

process graph_align_reads {
  publishDir "${params.out}/GraffiTE_alignments/", mode: 'copy'
  input:
  tuple val(sample_name), path(sample_reads), val(preset), path("index")
  val(graph_method)

  output:
  tuple val(sample_name), path("${sample_name}.gaf.gz"), path("${sample_name}.pack")

  script:

  interleaved = "-i"
  if(preset != "default") {
    interleaved = ""
  }

  switch(graph_method) {
    case "giraffe":
      """
      vg giraffe --parameter-preset ${preset} -o gam -t ${task.cpus} --index-basename index/index ${interleaved} -f ${sample_reads} > ${sample_name}.gam
      vg pack -x index/index.giraffe.gbz -g ${sample_name}.gam -o ${sample_name}.pack -Q ${params.min_mapq}
      vg convert -G ${sample_name}.gam index/index.giraffe.gbz | subset_gaf.py | sort -k1b,1 | gzip > ${sample_name}.gaf.gz
      rm ${sample_name}.gam
      """
      break
    case "graphaligner":
      """
      GraphAligner -t ${task.cpus} -x vg -g index/index.gfa -f ${sample_reads} -a ${sample_name}.gam
      vg pack -x index/index.gfa -g ${sample_name}.gam -o ${sample_name}.pack -Q ${params.min_mapq}
      vg convert -G ${sample_name}.gam index/index.gfa | subset_gaf.py | sort -k1b,1 | gzip > ${sample_name}.gaf.gz
      rm ${sample_name}.gam
      """
      break
  }
}

process vg_call {
  input:
  tuple val(sample_name), path(gaf), path(pack), path("index")
  val(graph_method)

  output:
  tuple val(sample_name), path("${sample_name}.vcf.gz*")

  script:
  def graph = graph_method == "giraffe" ? 'index.giraffe.gbz' : 'index.gfa'
  """
  vg call -a -A --threads ${task.cpus} -R chrX:1,chrY:1 -m ${params.min_support} -r index/index.pb -s ${sample_name} -k ${pack} index/${graph} | \
    bcftools norm -m-  | \
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
  bcftools merge -m none -l vcf.list > GraffiTE.merged.genotypes.vcf
  bgzip GraffiTE.merged.genotypes.vcf
  tabix -p vcf GraffiTE.merged.genotypes.vcf.gz
  grep '#' ${pangenome_vcf} > P_header
  grep -v '#' ${pangenome_vcf} | sort -k1,1 -k2,2n > P_sorted_body
  cat P_header P_sorted_body > pangenome.sorted.vcf
  bgzip pangenome.sorted.vcf
  tabix -p vcf pangenome.sorted.vcf.gz
  bcftools annotate -a pangenome.sorted.vcf.gz -c CHROM,POS,ID,REF,ALT,INFO GraffiTE.merged.genotypes.vcf.gz > GraffiTE.merged.genotypes.vcf
  awk -v v="${params.graffite_version}" 'NR==1 && /^##fileformat/ {print; print "##GraffiTE_version="v; next} {print}' GraffiTE.merged.genotypes.vcf > GraffiTE.merged.genotypes.vcf.tmp && mv GraffiTE.merged.genotypes.vcf.tmp GraffiTE.merged.genotypes.vcf
  rm -f GraffiTE.merged.genotypes.vcf.gz
  bgzip GraffiTE.merged.genotypes.vcf
  """
}
