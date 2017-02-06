##jmd
##11.28.11
##summarize_model.r
##for model overview table in paper
#want this to be consistent w/ cobra, maybe better to do some of this there

#run header of farm_script.r
source('/msc/neurospora/FBA/farm/farm_header.r')

comps <- unlist(unique(apply(as.matrix(c(rownames(s.al), colnames(s.al))), 1, FUN=function(x){ 
    if (length(grep('[', x, fixed=TRUE))>0){
        unlist(strsplit(x, split='\\[|\\]'))[2]
    }
})))
comps <- c('CCO-CYT', comps)

##genes & enzymes
gpr <- gpr[gpr$RXN %in% colnames(s.al),]
n2r <- NCUvector2rxns(unique(gpr$Genes), gpr) #spits out number of isozymes
length(unique(gpr$Genes))
length(unique(gpr$Enzyme))
#complexes start w/ 'CPLX'. (tried looking for enzymes associated w/ >1 gene, but complexes can be homo-mers.)
cmplx <- unique(grep('^CPLX', gpr$Enzyme, value=TRUE))
length(cmplx)
#number of unique genes in a complex
length(unique(gpr$Gene[gpr$Enzyme %in% cmplx]))

##rxns
#to get n rev rxns, can't use rxns.al$rev, since this computes reversibility from universal matrix
rev.model <- irrev2rev(s.al)
s.rev <- rev.model$sp
dim(s.rev)

##compartmentalized rxns
#to find number of rxns wholly in compartment x, count rxns that have no mets in any compartment other than x
colnames(s.al)[colSums(abs(s.al))==0] #nunca
rxns.comp <- list(); rxns.comp[comps] <- NA
rxns.comp[['CCO-CYT']] <- sum(colSums(abs(s.rev[grep('\\[CCO-', rownames(s.rev)),])) == 0)
for (i in 2:length(rxns.comp)){
    regexp <- paste('\\[', names(rxns.comp)[i], '\\]', sep='')
    rxns.comp[[i]] <- sum(colSums(abs(s.rev[-grep(regexp, rownames(s.rev)),])) == 0)
}
rxns.comp[['CCO-PM-FUNGI']] <- length(grep('[CCO-PM-FUNGI]', colnames(s.al), fixed=TRUE, value=TRUE))
rxns.comp[['CCO-VAC-MEM']] <- length(grep('[CCO-VAC-MEM]', colnames(s.al), fixed=TRUE, value=TRUE))

##exchange
#these are rxns that only give or only take from outside world
sum((colSums(sp.al==1)==1|colSums(sp.al==-1)==1) & colSums(sp.al==0)==nrow(sp.al)-1)

##trans rxns
#s.trans holds only rxns from 2 or more compartments
dim(s.trans <- s.rev[,!(cyt|mit|ex|glyox|nuc|vac)])
#keep all compartmental rxns except those from a given compartment, resulting in temporary S matrix w/ only other compartments (not inc. cyt), 
#then count number of columns that have a non-zero element.
sum(colSums(abs( s.trans[setdiff(grep('\\[CCO-', rownames(s.trans)), grep('\\[CCO-MIT\\]', rownames(s.trans))),] )) == 0)
sum(colSums(abs( s.trans[setdiff(grep('\\[CCO-', rownames(s.trans)), grep('\\[CCO-EXTRACELLULAR\\]', rownames(s.trans))),] )) == 0)
sum(colSums(abs( s.trans[setdiff(grep('\\[CCO-', rownames(s.trans)), grep('\\[CCO-GLYOXYSOME\\]', rownames(s.trans))),] )) == 0)

##mets
all(rowSums(abs(s.rev))>0)
#count number of unique mets (treating one met in diff compartments as one met)
length(unique(gsub('\\[CCO-.+', '', rownames(s.rev))))
#by compartment
mets.comp <- list(); mets.comp[comps] <- NA
mets.comp[['CCO-CYT']] <- nrow(s.rev)-length(grep('\\[CCO-', rownames(s.rev), value=TRUE))
for (i in 2:length(mets.comp)){
    regexp <- paste('\\[', names(mets.comp)[i], '\\]', sep='')
    mets.comp[[i]] <- length(grep(regexp, rownames(s.rev), value=TRUE))
}
#number of metabolites w/ structure (treating one met in diff compartments as one met)
smm.nc <- read.delim('/msc/neurospora/FBA/farm_data/Neurospora/nc10.smm', na='NIL')
smm <- smm.nc[smm.nc$FRAME %in% rownames(s.al) & !duplicated(smm.nc$FRAME),]
sum(smm$InChI!='')

##paths
nc.pwy2 <- nc.pwy[nc.pwy$rxn %in% colnames(s.al),]

##count number of paths
length(unique(nc.pwy2$PathwayID))
#count number of rxns in paths
length(unique(gsub('-L2R$|-R2L$', '', nc.pwy2$rxn)))

##pathway summary fig
setwd('Z:/reconstruct')
p <- read.delim('nrxns-of-level2-pwys.tsv', as.is=TRUE)
p$Level2Class <- gsub(paste(p$Level1Class, collapse='|'), '', p$Level2Class)
p <- p[p$Nrxns>2,]
p <- p[order(p[,1], -p$Nrxns),]
#color by level 1 class
rain <- rainbow(n=length(unique(p[,1]))) #c('black', 'white', 'grey')
col.v <- rain[as.numeric(as.factor(p[,1]))]

png('functional_classes.png')
par(mar=c(20, 4, 4, 2) + 0.1)
barplot(p$Nrxns, names.arg=p[,2], las=3, col=col.v)
legend(fill=rain, legend=levels(as.factor(p[,1])), x='topright')
dev.off()
