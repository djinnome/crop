##jmd
##2.13.11
##make_biomass.r

options(stringsAsFactors=FALSE)

setwd('/home/radon00/jdreyf/farm')
#setwd('Z:/meta_recon')

#biomass unbalanced: needs h2o and h+, like ecoli_iaf

##make biomass
bm <- read.delim('biomass.txt', header=FALSE)
bm.grps <- gsub('-->| \\+', '', readLines('biomass_notes.txt', n=1))
bm.grps <- data.frame(matrix(unlist(strsplit(bm.grps, split=' ')), ncol=2, byrow=TRUE))
bm.grps[,1] <- as.numeric(gsub('x', '50', bm.grps[,1]))
#coeffs
match.inds <- match(tolower(substr(bm[,3],1,3)), tolower(substr(bm.grps[,2],1,3)))
coeff.grp <- bm.grps[match.inds, 1]
coeff.grp[is.na(coeff.grp)] <- 1 #these don't have grps, so coeff=1*bm[,2]
coeff <- coeff.grp*bm[,2]
#ecoli_iaf similarly had most obj coeffs < 1 in magnitude but atp, adp, & p (also h2o & h) > 50
coeff[bm[,1]=='ATP'] <- 50
#sparse S
bm.s <- data.frame(Rxn.name='biomass', Compound=c(bm[,1], 'ADP', 'Pi'), coefficient=c(-coeff, 50, 50))
#write.table(bm.s, 'biomass_s.txt', sep='\t', row.names=FALSE, quote=FALSE)
