##jmd
##10.18.11
##fix_bm_coeffs.r

smm <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.smm', na='NIL')
rownames(smm) <- smm$FRAME

##% mass to stoich
#(bm2 <- bm[bm$rxn=='SecondaryMetaboliteComposition' & bm$coeff<0,])
#(mass <- smm[bm2$comp, 'MASS']/1000)
#(cf2  <-  -bm2$coeff / mass); names(cf2) <- bm2$comp
#round(cf2,2) %*% mass

##make coeffs specifically for FullBiomassComposition
# (bm2=bm[intersect(which(bm$rxn=='biomass'), grep('_bm$', bm$comp)),])
# cf2=c(-bm2$coeff, 0.0084); names(cf2)=c(bm2$comp, 'Secondary-Metabolites_bm')
# cf2['Lipids_bm']=0.1116
# write.table(data.frame('( mix2t-27', names(cf2), paste(cf2, ')')), file='', quote=FALSE, row.names=FALSE)

##get mass
(y <- get.bm.mass(bm, smm))

##fix a particular set
set <- 'LipidComposition'
bm2 <- bm[bm$rxn %in% set & bm$coeff<0,]
cf <- -bm2$coeff
mets <- gsub('^Charged-|-tRNAs$', '', bm2$comp) 
mass <- smm[mets, 'MASS']/1000
#if mass is NA, assume it's 1, since we don't have any instances in biomass composition that aren't in smm
mass[is.na(mass)] <- 1
cf %*% mass

##get coeffs to make 1 gram
cf2 <- cf / (cf %*% mass); names(cf2) <- mets
round(cf2,2) %*% mass

##write to console for lisp
write.table(data.frame('( mix2t-27', names(cf2), paste(round(cf2,2), ')')), file='', quote=FALSE, row.names=FALSE)

#########################################################################################################################
####LISP!
#########################################################################################################################

(add-slot-value 'mix2t-27 'composition 'SPHINGOLIPIDS)

(setq composition 
    '( 
        ( mix2t-27 |Lipids| 0.1082 )
        ( mix2t-27 SPHINGOLIPIDS 0.00335 )
    )
)

(loop for (bm cpd coeff) in composition do
    (put-value-annot bm 'composition cpd 'coefficient coeff)
)
