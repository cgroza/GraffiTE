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

include { index_graph; bamtags_to_bed; epigenome_to_CSV; annotate_VCF } from './panmethyl/module/'

include { break_scaffold; map_asm; map_longreads; sniffles_sample_call; sniffles_population_call;
         svim_asm; truvari_merge; split_repeatmask; concat_repeatmask; repeatmask_VCF; tsd_prep;
         tsd_search; tsd_report; pangenie_index; pangenie; make_graph; bam_to_fastq;
         graph_align_reads; vg_call; merge_VCFs } from './module'

workflow {
  // initiate channels that will provide the reference genome to processes
  Channel.fromPath(params.reference, checkIfExists:true).set{ref_asm_ch}

  if(!params.graffite_vcf && !params.vcf && !params.RM_dir) {
    svim_variants_ch = channel.empty()
    sn_variants_ch = channel.empty()

    if(params.longreads || params.bams) {
      sniffles_reads_in_ch = channel.empty()
      sniffles_bams_in_ch = channel.empty()

      if(params.longreads) {
        Channel.fromPath(params.longreads).splitCsv(header:true).map{row ->
          [row.sample, file(row.path, checkIfExists:true), row.type]}.combine(ref_asm_ch).set{map_longreads_in_ch}
        map_longreads(map_longreads_in_ch).set{sniffles_reads_in_ch}
      }

      if(params.bams) {
        sniffles_bams_in_ch = Channel.fromPath(params.bams).splitCsv(header:true).map{row ->
          [row.sample, file(row.path, checkIfExists:true)]}.combine(ref_asm_ch)
      }


      sniffles_population_call(
        sniffles_sample_call(
          sniffles_reads_in_ch.concat(sniffles_bams_in_ch)).map{it -> it[0]}.collect(),
        ref_asm_ch).flatten().set{sn_variants_ch}
    }

    if(params.assemblies) {
      Channel.fromPath(params.assemblies).splitCsv(header:true).map{row ->
        [row.sample, file(row.path, checkIfExists:true)]}.set{map_asm_in_ch}
      if(params.break_scaffolds) {
        map_asm_in_ch = break_scaffold(map_asm_in_ch)
      }

      svim_asm(map_asm(map_asm_in_ch.combine(ref_asm_ch))).map{sample -> sample[1]}.set{svim_variants_ch}
    }

    truvari_merge(svim_variants_ch.mix(sn_variants_ch).collect(), ref_asm_ch).set{sv_variants_ch}
  }

  // if the user doesn't provide a VCF already made by GraffiTE with --graffite_vcf, use RepeatMasker to annotate repeats
  if(!params.graffite_vcf) {
    // except if --RM_dir is given, in which case skip RepeatMasker here and set the input channel
    RM_ch = channel.empty()
    if(params.RM_dir){
      channel.fromPath("${params.RM_dir}/*", type: "dir").
      map{p -> [file("${p}/genotypes_repmasked_filtered.vcf", checkIfExists: true), file("${p}/repeatmasker_dir", checkIfExists: true)]}.
      map{v -> [v[0], v[1]]}.set{RM_ch}
    } else {
      Channel.fromPath(params.TE_library, checkIfExists:true).set{TE_library_ch}
      // we need to set the vcf input depending what was given
      if(params.longreads || params.bams || params.assemblies){
        sv_variants_ch.set{raw_vcf_ch}
      } else if(params.vcf){
        Channel.fromPath(params.vcf, checkIfExists : true).set{raw_vcf_ch}
      } else {
        error "No --longreads, --assemblies, --vcf or --RM_dir parameters passed to GraffiTE."
      }
      repeatmask_VCF(split_repeatmask(raw_vcf_ch).flatten().combine(TE_library_ch).combine(ref_asm_ch)).set{RM_ch}
    }
    tsd_report(tsd_search(tsd_prep(RM_ch.combine(ref_asm_ch)).
                          splitText(elem: 3, by: params.tsd_batch_size, file: true)).
               groupTuple(by: 3).
               map{v -> tuple(v[0], v[1], v[2][0], v[3])}
    )
    concat_repeatmask(tsd_report.out.vcf_ch.collect(),
                      tsd_report.out.tsd_full_group_ch.collect(),
                      tsd_report.out.tsd_sum_group_ch.collect(),
                      ref_asm_ch)
    concat_repeatmask.out.vcf_ch.set{vcf_ch}
  } else {
    // if a vcf is provided as parameter, skip discovery and go directly to genotyping
    Channel.fromPath(params.graffite_vcf).set{vcf_ch}
  }

  if(params.genotype) {
    Channel.fromPath(params.genotype_with).splitCsv(header:true).map{ row ->
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
    }.branch{ it ->
        bam: it[1].extension == "bam"
        fastq: it[1].extension != "bam"
    }.set{reads_input_ch}

    reads_input_ch.fastq.mix(bam_to_fastq(reads_input_ch.bam)).set{reads_ch};

    indexed_vcfs = channel.empty()
    if(params.graph_method == "pangenie") {
      reads_ch.combine(pangenie_index(vcf_ch.combine(ref_asm_ch))).set{input_ch}
      pangenie(input_ch, ref_asm_ch).set{indexed_vcfs}
    } else if(params.graph_method == "giraffe" || params.graph_method == "graphaligner") {
      graph_method = channel.value(params.graph_method)

      make_graph(vcf_ch, ref_asm_ch, graph_method).set{graph_index_ch}
      reads_ch.combine(graph_index_ch).set{reads_align_ch}
      graph_align_reads(reads_align_ch, graph_method).set{aligned_ch}
      aligned_ch.combine(graph_index_ch).set{graph_pack_ch}
      vg_call(graph_pack_ch, graph_method).set{indexed_vg_call_vcfs}

      if(params.epigenomes) {
        reads_input_ch.bam.map{row -> [row[0], row[1]]}.set{epigenome_ch}

        index_graph(graph_index_ch.map(p -> p / 'index.gfa'),
                    channel.value(params.motif)).set{indexed_graph_ch}

        bamtags_to_bed(
          epigenome_ch.combine(aligned_ch.map{it -> [it[0], it[1]]}, by: 0)
            .combine(indexed_graph_ch),
          channel.value(params.tag),
          channel.value(params.missing_modifications)).set{mods_ch}

        epigenome_to_CSV(mods_ch.combine(indexed_graph_ch)).set{mods_csv_ch}
        annotate_VCF(indexed_vg_call_vcfs.map{v -> [v[0], v[1][0]]}.combine(mods_csv_ch, by: 0)).map{it -> [it[0], it[1]]}.set{indexed_vcfs}
      } else {
        indexed_vg_call_vcfs.set{indexed_vcfs}
      }
    } else {
      error "Unsupported --graph_method. --graph_method must be pangenie, giraffe or graphaligner."
    }

    merge_VCFs(indexed_vcfs.map{v -> v[1]}.collect(), vcf_ch)
  }
}
