vg chunk -S Zea_mays.pb -p CM039158.1:13076986-13149277 -x Zea_mays.vg | vg view - > Zea_bz.gfa
awk '$1 == "S" {print $2, length($3)}' Zea_bz.gfa > node_sizes.csv
vg paths -F -x Zea_bz.gfa > Zea_bz.fasta
vg paths -A -x Zea_bz.gfa > Zea_bz.gaf
singularity exec --bind $(pwd):~ graffite_latest.sif RepeatMasker -s -nolow -lib NAM.EDTA2.0.0.MTEC02052020.TElib.fa Zea_bz.fasta
Rscript join_annotation.R  Zea_bz.fasta.out Zea_bz.gaf node_sizes.csv
