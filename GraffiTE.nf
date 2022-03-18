params.vcf = "variants.vcf"
params.reads = "reads.csv"

Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [file(row.path), val(row.sample)]}.view().set{reads_ch}

