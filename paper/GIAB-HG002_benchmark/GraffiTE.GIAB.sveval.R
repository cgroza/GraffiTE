#BiocManager::install('jmonlong/sveval')
library(sveval)
library(ggplot2)
# load the vcfs paths

## Truth GIAB
truth.vcf='sveval_VCFS/GIAB.Tier1.Alu.L1.SVA.250bp.vcf'

## GT
# pangenomes
gt.svsn05x='sveval_VCFS/GT.svim-sniffles05x.Alu.L1.SVA.250.vcf'
gt.svsn10x='sveval_VCFS/GT.svim-sniffles10x.Alu.L1.SVA.250.vcf'
gt.svsn20x='sveval_VCFS/GT.svim-sniffles20x.Alu.L1.SVA.250.vcf'
gt.svsn30x='sveval_VCFS/GT.svim-sniffles30x.Alu.L1.SVA.250bp.vcf'
gt.sn05x="sveval_VCFS/GT.sniffles05x.Alu.L1.SVA.250bp.vcf"
gt.sn10x="sveval_VCFS/GT.sniffles10x.Alu.L1.SVA.250bp.vcf"
gt.sn20x="sveval_VCFS/GT.sniffles20x.Alu.L1.SVA.250bp.vcf"
gt.sn30x='sveval_VCFS/GT.sniffles30x.Alu.L1.SVA.250bp.vcf'
gt.sv='sveval_VCFS/GT.svim.Alu.L1.SVA.250bp.vcf'

# pangenie @ 30X
gt.svsn30x.pan30X='sveval_VCFS/HG002_mat-pat.svim_30Xhifi.sniffles_30Xhifi.pangenie.genotypes.Alu.L1.SVA.250bp.vcf'
gt.sn30x.pan30X='sveval_VCFS/HG002_30Xhifi.sniffles_30Xhifi.pangenie.genotypes.Alu.L1.SVA.250bp.vcf'
gt.sv.pan30X='sveval_VCFS/HG002_mat-pat.svim.pangenie.genotypes.Alu.L1.SVA.250bp.vcf'
# graphaligner @30X
gt.svsn30x.gra30X='sveval_VCFS/HG002_mat-pat.svim_30Xhifi.sniffles_30Xhifi.graphaligner.genotypes.Alu.L1.SVA.250bp.vcf'
gt.sn30x.gra30X='sveval_VCFS/HG002_30Xhifi.sniffles_30Xhifi.graphaligner.genotypes.Alu.L1.SVA.250bp.vcf'
gt.sv.gra30X='sveval_VCFS/HG002_mat-pat.svim.graphaligner.genotypes.Alu.L1.SVA.250bp.vcf'

## TLDR
tl.05X="sveval_VCFS/TLDR.05X_mean_.vcf"
tl.10X="sveval_VCFS/TLDR.10X_mean_.vcf"
tl.20X="sveval_VCFS/TLDR.20X_mean_.vcf"
tl.30X="sveval_VCFS/TLDR.30X_mean_.vcf"
## MELT2
m2.test05X="sveval_VCFS/M2.Alu.L1.SVA.05X.vcf"
m2.test10X="sveval_VCFS/M2.Alu.L1.SVA.10X.vcf"
m2.test20X="sveval_VCFS/M2.Alu.L1.SVA.20X.vcf"
m2.test30X="sveval_VCFS/M2.Alu.L1.SVA.30X.vcf"

## MEGAnE
mg.05X="sveval_VCFS/MG.Alu.L1.SVA.05X.vcf"
mg.10X="sveval_VCFS/MG.Alu.L1.SVA.10X.vcf"
mg.20X="sveval_VCFS/MG.Alu.L1.SVA.20X.vcf"
mg.30X="sveval_VCFS/MG.Alu.L1.SVA.30X.vcf"

########## pME detection evaluation (presence/absence)

# GT svim_mat.pat + snif_05X 
eval.svsn05x = svevalOl(gt.svsn05x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
# GT svim_mat.pat + snif_30X 
eval.svsn10x = svevalOl(gt.svsn10x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
# GT svim_mat.pat + snif_30X 
eval.svsn20x = svevalOl(gt.svsn20x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
# GT svim_mat.pat + snif_30X 
eval.svsn30x = svevalOl(gt.svsn30x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.svsn30x$eval # data.frame with results using all variants
# GT svim_mat.pat ONLY (from 2 alt asms)
eval.sv = svevalOl(gt.sv, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.sv$eval
# GT snif_05X ONLY (only from long reads)
eval.sn05 = svevalOl(gt.sn05x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.sn05$eval
# GT snif_10X ONLY (only from long reads)
eval.sn10 = svevalOl(gt.sn10x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.sn10$eval
# GT snif_20X ONLY (only from long reads)
eval.sn20 = svevalOl(gt.sn20x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.sn20$eval
# GT snif_30X ONLY (only from long reads)
eval.sn30 = svevalOl(gt.sn30x, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.sn30$eval

## TLDR
eval.tl.05=svevalOl(tl.05X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.tl.05$eval
eval.tl.10=svevalOl(tl.10X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.tl.10$eval
eval.tl.20=svevalOl(tl.20X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.tl.20$eval
eval.tl.30=svevalOl(tl.30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.tl.30$eval

## MELT
eval.m2.05=svevalOl(m2.test05X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.m2.05$eval
eval.m2.10=svevalOl(m2.test10X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.m2.10$eval
eval.m2.20=svevalOl(m2.test20X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.m2.20$eval
eval.m2.30=svevalOl(m2.test30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.m2.30$eval

## MEGAnE
eval.mg.05=svevalOl(mg.05X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.mg.05$eval
eval.mg.10=svevalOl(mg.10X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.mg.10$eval
eval.mg.20=svevalOl(mg.20X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.mg.20$eval
eval.mg.30=svevalOl(mg.30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed')
#eval.mg.30$eval

### Graph
# make a summary table
pMEevals=rbind(cbind(eval.sv$eval, data.frame(cov = c(0,0,0),
                                                method = rep("GT-sv", 3))
),
cbind(eval.svsn05x$eval, data.frame(cov = c(5,5,5),
                                    method = rep("GT-sv-sn", 3))
),
cbind(eval.svsn10x$eval, data.frame(cov = c(10,10,10),
                                 method = rep("GT-sv-sn", 3))
),
cbind(eval.svsn20x$eval, data.frame(cov = c(20,20,20),
                                 method = rep("GT-sv-sn", 3))
),
cbind(eval.svsn30x$eval, data.frame(cov = c(30,30,30),
                                 method = rep("GT-sv-sn", 3))
),
cbind(eval.sn05$eval, data.frame(cov = c(5,5,5),
                                 method = rep("GT-sn", 3))
),
cbind(eval.sn10$eval, data.frame(cov = c(10,10,10),
                                 method = rep("GT-sn", 3))
),
cbind(eval.sn20$eval, data.frame(cov = c(20,20,20),
                                 method = rep("GT-sn", 3))
),
cbind(eval.sn30$eval, data.frame(cov = c(30,30,30),
                                 method = rep("GT-sn", 3))
),
cbind(eval.m2.05$eval, data.frame(cov = c(5,5,5),
                                 method = rep("MELT2", 3))
),
cbind(eval.m2.10$eval, data.frame(cov = c(10,10,10),
                                  method = rep("MELT2", 3))
),
cbind(eval.m2.20$eval, data.frame(cov = c(20,20,20),
                                  method = rep("MELT2", 3))
),
cbind(eval.m2.30$eval, data.frame(cov = c(30,30,30),
                                  method = rep("MELT2", 3))
),
cbind(eval.mg.05$eval, data.frame(cov = c(5,5,5),
                                  method = rep("MEGAnE", 3))
),
cbind(eval.mg.10$eval, data.frame(cov = c(10,10,10),
                                  method = rep("MEGAnE", 3))
),
cbind(eval.mg.20$eval, data.frame(cov = c(20,20,20),
                                  method = rep("MEGAnE", 3))
),
cbind(eval.mg.30$eval, data.frame(cov = c(30,30,30),
                                  method = rep("MEGAnE", 3))
),
cbind(eval.tl.05$eval[2,], data.frame(cov = 5),
      method = "TLDR"
),
cbind(eval.tl.10$eval[2,], data.frame(cov = 10),
      method = "TLDR"
),
cbind(eval.tl.20$eval[2,], data.frame(cov = 20),
      method = "TLDR"
),
cbind(eval.tl.30$eval[2,], data.frame(cov = 30),
      method = "TLDR"
)
)

ggplot(pMEevals[pMEevals$cov > 0,], aes(x = precision, y = recall))+
  geom_point(aes(size = cov, fill = method), pch = 21, col = "white")+
  scale_fill_manual(values = c("#FF214B", "#AD000D", "#7C55E9", "#59206F", "#47E075"))+
  geom_path(aes(col = method))+
  scale_color_manual(values = c("#FF214B", "#AD000D", "#7C55E9", "#59206F", "#47E075"))+
  geom_point(data = pMEevals[pMEevals$cov == 0,], aes(x = precision, y = recall, fill = method), fill = "#FF90B9", pch = 24, size = 5, col = "white")+
  facet_wrap(~type)+
  theme_classic()

SVplot<-ggplot(pMEevals[pMEevals$cov > 0,], aes(x = precision, y = recall))+
  geom_point(aes(size = cov, fill = method), pch = 21, col = "white", show.legend = F)+
  scale_fill_manual(values = c("#FF214B", "#AD000D", "#7C55E9", "#59206F", "#47E075"))+
  geom_path(aes(col = method), show.legend = F)+
  scale_color_manual(values = c("#FF214B", "#AD000D", "#7C55E9", "#59206F", "#47E075"))+
  geom_point(data = pMEevals[pMEevals$cov == 0,], aes(x = precision, y = recall, fill = method), fill = "#FF90B9", pch = 24, size = 5, col = "white")+
  facet_wrap(~type)+
  theme_classic()

########### GENOTYPES

#GT pangenie 30X
eval.svsn.p30.geno = svevalOl(gt.svsn30x.pan30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                           geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.svsn.p30.geno$eval # data.frame with results using all variants
eval.sv.p30.geno = svevalOl(gt.sv.pan30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                              geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.sv.p30.geno$eval # data.frame with results using all variants
eval.sn.p30.geno = svevalOl(gt.sn30x.pan30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                            geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.sn.p30.geno$eval # data.frame with results using all variants
#GT graphaligner 30X
eval.svsn.g30.geno = svevalOl(gt.svsn30x.gra30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                              geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.svsn.g30.geno$eval # data.frame with results using all variants
eval.sv.g30.geno = svevalOl(gt.sv.gra30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                            geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.sv.g30.geno$eval # data.frame with results using all variants
eval.sn.g30.geno = svevalOl(gt.sn30x.gra30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                            geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.sn.g30.geno$eval # data.frame with results using all variants
# MELT2 30X
eval.m2.30.geno=svevalOl(m2.test30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                         geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.m2.30.geno$eval
# MEGAnE 30X
eval.mg.30.geno=svevalOl(mg.30X, truth.vcf, bed.regions = 'HG002_SVs_Tier1_v0.6.bed',
                         geno.eval=TRUE, method="bipartite", stitch.hets=TRUE, merge.hets=FALSE)
eval.mg.30.geno$eval

### Graph
# make a summary table
pMEevals.geno=rbind(cbind(eval.sn.p30.geno$eval, data.frame(geno = rep("Pangenie",3),
                                 SV_method = rep("GT-sniffles", 3))
),
cbind(eval.svsn.p30.geno$eval, data.frame(geno = rep("Pangenie",3),
                                        SV_method = rep("GT-svim-sniffles", 3))
),
cbind(eval.sv.p30.geno$eval, data.frame(geno = rep("Pangenie",3),
                                          SV_method = rep("GT-svim", 3))
),
cbind(eval.sn.g30.geno$eval, data.frame(geno = rep("GraphAligner",3),
                                        SV_method = rep("GT-sniffles", 3))
),
cbind(eval.svsn.g30.geno$eval, data.frame(geno = rep("GraphAligner",3),
                                          SV_method = rep("GT-svim-sniffles", 3))
),
cbind(eval.sv.g30.geno$eval, data.frame(geno = rep("GraphAligner",3),
                                        SV_method = rep("GT-svim", 3))
),
cbind(eval.m2.30.geno$eval, data.frame(geno = rep("MELT2",3),
                                        SV_method = rep("MELT2", 3))
),
cbind(eval.mg.30.geno$eval, data.frame(geno = rep("MEGAnE",3),
                                       SV_method = rep("MEGAnE", 3))
)
)

ggplot(pMEevals.geno, aes(x = precision, y = recall))+
  #geom_point(aes(pch = geno, col = SV_method), size =5)+
  geom_point(aes(col = paste(SV_method, geno, sep ="+"), pch = paste(SV_method, geno, sep ="+")), size =3)+
  scale_shape_manual(values = c(21,19,21,19,21,19,19,19))+
  scale_color_manual(values = c( "#FF214B", "#FF214B", "#AD000D", "#AD000D", "#FF90B9", "#FF90B9", "#7C55E9", "#59206F"))+
  #geom_path(aes(col = method))+
  facet_wrap(~type)+
  theme_classic()

GENOplot<-ggplot(pMEevals.geno, aes(x = precision, y = recall))+
  #geom_point(aes(pch = geno, col = SV_method), size =5)+
  geom_point(aes(col = paste(SV_method, geno, sep ="+"), pch = paste(SV_method, geno, sep ="+")), size =3, show.legend = F)+
  scale_shape_manual(values = c(21,19,21,19,21,19,19,19))+
  scale_color_manual(values = c( "#FF214B", "#FF214B", "#AD000D", "#AD000D", "#FF90B9", "#FF90B9", "#7C55E9", "#59206F"))+
  #geom_path(aes(col = method))+
  facet_wrap(~type)+
  theme_classic()

cowplot::plot_grid(SVplot, GENOplot, ncol = 1)
