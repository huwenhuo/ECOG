options(width = 222)
library(data.table)
library(ggpubr)
library(survminer)
library(survival)

setwd('~/gdrive/huw/chung/ecog/')

library(org.Hs.eg.db)


gpl = fread("GPL6602-433.txt", header=T)
setnames(gpl, 'ID', 'NimblGen')
#colnames(gpl) = c('NimblGen', 'ACCNUM')

res = select(org.Hs.eg.db, keys=gpl$GB_ACC, columns=c("SYMBOL", "GENENAME", "ENSEMBL", "ACCNUM"), keytype="ACCNUM")
head(res)
dim(res)
gpl.id = merge(gpl, res, by.x = 'GB_ACC', by.y = "ACCNUM") 
#gpl.id[SYMBOL == 'BCAT2',]
# SLC7A5 (LAT1), SLC3A2, SLC16A19, SLC16A1 (MCT1), SLC16A3 (MCT4), SLC25A44
genes = c('SLC7A5', 'SLC3A2', 'SLC16A19', 'SLC16A1', 'SLC16A3', 'SLC25A44')
gpl.id[SYMBOL %in% genes, ]
id = unique(gpl.id[SYMBOL %in% genes, NimblGen])
id

ecog.dt = fread('ecog1900_expression_outcomes.txt')
# genotype mutation
ecog.dt[1:10, 1:33]
# clinical information
ecog.dt[1:10, 34:59]
# gene expression, NimbelGen
ecog.dt[1:10, 60:89]
table(ecog.dt$tumor_normal)

mtx.1 = ecog.dt[6:324, 5:8]
mtx.2 = ecog.dt[6:324, 9:39]
mtx.3 = ecog.dt[6:324, id, with=F]
mtx = cbind(mtx.1, mtx.2, mtx.3)
mtx[1:10, ]
mtx.3[1:10, ]

lapply(genes, function(gene){genex.fun(gene, gpl.id = gpl.id, ecog.dt = ecog.dt)})

genex.fun = function(genex, gpl.id = gpl.id, ecog.dt = ecog.dt, fname.pre = NULL){

	if(is.null(fname.pre)){ fname.pre = paste0(getwd(), '/ecog') }

	probe.id = unique(gpl.id[SYMBOL == genex, NimblGen])
	if(length(probe.id) == 0){
		message('No genes found for gene ', genex)
		return(0)}
	probe.id
	tmp = ecog.dt[6:324, 9:39]
	tmp[1:10, ]
	for(ii in colnames(tmp)){
		tmp[get(ii) != 'WT', (ii) := 'MT']
		tmp[, (ii) := factor(get(ii), levels = c('WT', 'MT'))]
	}
	tmp2 = cbind(tmp, ecog.dt[6:324, probe.id, with=F])
	tmp2[1:10, ]
	tmp2

	genotype.genes = colnames(tmp2)[1:31]
	genotype.genes
	gg.list = list()
	for(genotype.gene in genotype.genes){
		gg = lapply(probe.id, function(xx){ ggboxplot(tmp2, legend.title = genotype.gene, 
							      x = genotype.gene, y = xx, color = genotype.gene, add = 'jitter') + 
			    stat_compare_means() + xlab(genotype.gene) })
		gg.list = c(gg.list, gg)
	}

	gg = gridExtra::marrangeGrob(gg.list, ncol=6, nrow=length(probe.id), as.table = FALSE)
	fname = paste0(fname.pre, '_', genex, '.pdf'); fname
	cat('plotting ', fname, '\n')
	ggplot2::ggsave(gg, file=fname, width=13, height=3*length(probe.id))

}


## BCAT1

{

	genex = 'BCAT1'
	probe.id = unique(gpl.id[SYMBOL == genex, NimblGen])
	tmp = ecog.dt[6:324, 9:39]
	tmp[1:10, ]
	for(ii in colnames(tmp)){
		tmp[get(ii) != 'WT', (ii) := 'MT']
		tmp[, (ii) := factor(get(ii), levels = c('WT', 'MT'))]
	}
	tmp2 = cbind(tmp, ecog.dt[6:324, id, with=F])
	tmp2[1:10, ]
	tmp2

	genes = colnames(tmp2)[1:31]
	genes
	gg.list = list()
	for(gene in genes){
		tmp = lapply(probe.id, function(xx){ ggboxplot(tmp2, legend.title = gene, x = gene, y = xx, color = gene, add = 'jitter') + stat_compare_means() + xlab(gene) })
		gg.list = c(gg.list, tmp)
	}

	gg = gridExtra::marrangeGrob(gg.list, ncol=6, nrow=length(probe.id), as.table = FALSE)
	fname = paste0('res/all_', genex, '.pdf'); fname
	cat('plotting ', fname, '\n')
	ggplot2::ggsave(gg, file=fname, width=13, height=3*length(probe.id))

}


# BCAT2
{

	genex = 'BCAT2'
	probe.id = unique(gpl.id[SYMBOL == genex, NimblGen])
	tmp = ecog.dt[6:324, 9:39]
	tmp[1:10, ]
	for(ii in colnames(tmp)){
		tmp[get(ii) != 'WT', (ii) := 'MT']
		tmp[, (ii) := factor(get(ii), levels = c('WT', 'MT'))]
	}
	tmp2 = cbind(tmp, ecog.dt[6:324, id, with=F])
	tmp2[1:10, ]
	tmp2

	genes = colnames(tmp2)[1:31]
	genes
	gg.list = list()
	for(gene in genes){
		tmp = lapply(probe.id, function(xx){ ggboxplot(tmp2, legend.title = gene, x = gene, y = xx, color = gene, add = 'jitter') + stat_compare_means() + xlab(gene) })
		gg.list = c(gg.list, tmp)
	}

	gg = gridExtra::marrangeGrob(gg.list, ncol=6, nrow=length(probe.id), as.table = FALSE)
	fname = paste0('res/all_', genex, '.pdf'); fname
	cat('plotting ', fname, '\n')
	ggplot2::ggsave(gg, file=fname, width=13, height=3*length(probe.id))

}

save.image()

fit <- survfit(Surv(OS, OSStatus) ~ sex, data = ecog.dt)


ggsurvplot(
  fit, 
  data = lung, 
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("Male", "Female"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

