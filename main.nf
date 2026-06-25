// SAY HELLO

// 1. Read the version from the local file
def versionFile = file("${baseDir}/version.txt")
def pipelineVersion = versionFile.exists() && versionFile.text.trim() ? versionFile.text.trim() : '1.1.0'

// Expose version to process scripts (used for stamping VCF headers)
params.graffite_version = pipelineVersion

// 2. Define the revision (branch name)
def pipelineRevision = workflow.revision ?: 'main'


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

V. ${pipelineVersion} - ${pipelineRevision}

Pangenomic Toolbox for the Analysis of Transposable Element Insertion Polymorphisms

Authors: Cristian Groza and Clément Goubert
Bug/issues: https://github.com/cgroza/GraffiTE/issues

"""

include { index_graph; bamtags_to_BED; lift_epigenome; annotate_VCF; annotate_BED; merge_BED; BED_to_graph; merge_CSV } from './panmethyl/module/'

include { break_scaffold; map_asm; map_longreads; sniffles_sample_call; sniffles_population_call;
         svim_asm; pav_asm; truvari_merge; split_repeatmask; concat_repeatmask; repeatmask_VCF; tsd_prep;
         tsd_search; tsd_report; pangenie_index; pangenie; make_graph; bam_to_fastq;
         graph_align_reads; vg_call; merge_VCFs } from './module'

workflow {
  // initiate channels that will provide the reference genome to processes
  Channel.fromPath(params.reference, checkIfExists:true).set{ref_asm_ch}

  if(!params.graffite_vcf && !params.vcf && !params.RM_dir) {
    svim_variants_ch = channel.empty()
    pav_variants_ch = channel.empty()
    sn_variants_ch = channel.empty()
    vcfs_variants_ch = channel.empty()

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

    if(params.pav) {
      Channel.fromPath(params.pav).splitCsv(header:false, skip:1).map{row ->
        [row[0], row[1..-1].collect({ file(it, checkIfExists:true) })]}.set{pav_in_ch}
      pav_asm(pav_in_ch.combine(ref_asm_ch)).set{pav_variants_ch}
    }

    if(params.svs) {
      Channel.fromPath(params.svs).splitCsv(header:true).map{row ->
        [row.sample, file(row.path, checkIfExists:true)]}.map{sample -> sample[1]}.set{vcfs_variants_ch}
    }

    truvari_merge(svim_variants_ch.mix(sn_variants_ch).mix(vcfs_variants_ch).mix(pav_variants_ch).collect(), ref_asm_ch, false).set{sv_variants_ch}
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
      if(params.longreads || params.bams || params.assemblies || params.pav || params.svs){
        sv_variants_ch.set{raw_vcf_ch}
      } else if(params.vcf){
        truvari_merge(Channel.fromPath(params.vcf, checkIfExists : true), ref_asm_ch, true).set{raw_vcf_ch}
      } else {
        error "No --longreads, --assemblies, --pav, --vcf or --RM_dir parameters passed to GraffiTE."
      }
      repeatmask_VCF(split_repeatmask(raw_vcf_ch).flatten().combine(TE_library_ch).combine(ref_asm_ch))
      repeatmask_VCF.out.vcf.set{RM_ch}
    }
    tsd_report(tsd_search(tsd_prep(RM_ch.combine(ref_asm_ch)).
                          splitText(elem: 3, by: params.tsd_batch_size, file: true)).
               map{it -> [it[0], it[1], it[2], it[3].getText()]}.
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
    } else if(params.graph_method == "giraffe" || params.graph_method == "graphaligner" || params.graph_method == "precomputed") {
      graph_method = channel.value(params.graph_method)

      graph_index_ch = channel.empty()
      if (params.graph) {
        Channel.fromPath(params.graph).set{graph_index_ch}
      } else {
        make_graph(vcf_ch, ref_asm_ch, graph_method).set{graph_index_ch}
      }

      indexed_vg_call_vcfs = channel.empty()

      if (params.vcfs) {
        Channel.fromPath(params.vcfs).splitCsv(header : true).map{
          row -> [row.sample, file(row.path, checkIfExists: true).toSorted()]}.set{indexed_vg_call_vcfs}
      } else {
        reads_ch.combine(graph_index_ch).set{reads_align_ch}
        graph_align_reads(reads_align_ch, graph_method).set{aligned_ch}
        aligned_ch.combine(graph_index_ch).set{graph_pack_ch}
        vg_call(graph_pack_ch, graph_method).set{indexed_vg_call_vcfs}
      }

      if(params.epigenomes) {
        index_graph(graph_index_ch.map(p -> p / 'index.gfa'),
                    channel.value(params.motif)).set{indexed_graph_ch}

        lifted_mods_ch = channel.empty()
        if (params.lifted) {
          Channel.fromPath(params.lifted).splitCsv(header : true)
            .map{row -> [row.sample, file(row.path, checkIfExists : true)]}.set{lifted_mods_ch}
        }
        else {
          reads_input_ch.bam.map{row -> [row[0], row[1]]}.set{epigenome_ch}

          bamtags_to_BED(epigenome_ch, channel.value(params.code)).set{mods_ch}

          lift_epigenome(mods_ch.combine(aligned_ch.map{it -> [it[0], it[1]]}, by: 0).combine(indexed_graph_ch)).set{lifted_mods_ch}

        }

        merge_CSV(lifted_mods_ch.groupTuple(by: 0).combine(indexed_graph_ch)).set{mods_csv_ch}
        annotate_VCF(indexed_vg_call_vcfs.map{v -> [v[0], v[1][0]]}.combine(mods_csv_ch, by: 0)).map{it -> [it[0], it[1]]}.set{indexed_vcfs}

        if(params.bed) {
          BED_to_graph(graph_index_ch.map{it -> it / "index.gfa"}.combine(Channel.fromPath(params.bed))).set{bed_ch}
          merge_BED(annotate_BED(mods_csv_ch.combine(bed_ch).combine(indexed_graph_ch.map{it[0]})).map{it[1]}.collect())

        }

      } else {
        indexed_vg_call_vcfs.set{indexed_vcfs}
      }
    } else {
      error "Unsupported --graph_method. --graph_method must be pangenie, giraffe or graphaligner."
    }

    merge_VCFs(indexed_vcfs.map{v -> v[1]}.collect(), vcf_ch)
  }
}
