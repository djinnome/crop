##jmd
##9.22.11
##check_compartments.r

source('/msc/neurospora/FBA/farm/farm_config.r')
source('/msc/neurospora/FBA/farm/phenos2rxns.r')

##read
#rxns annot
rxns <- read.delim('rxn_annot.txt')
rxns <- rxns[colnames(s.sp),]
#rename, for below
s.al <- sg
rxns.al <- rxns

##get model rxns + compartmentalized rxns
model.rxns <- read.csv('model_rxns.csv')[,1]
g <- grep('\\[CCO-MIT\\]$|\\[CCO-GLYOXYSOME\\]$', colnames(s.sp), value=TRUE)
model.rxns <- union(model.rxns, g)
#subset s matrix
ss <- intersect(model.rxns, colnames(sg))
s.al <- s.al[,ss]
s.al <- s.al[rowSums(abs(s.al))>0,]
rxns.al <- rxns.al[ss,]
nmets <- nrow(s.al); nrxns <- ncol(s.al)

##get rid of no flux rxns
max.nv <- max.n.flux(s.al, allow.all.trans=FALSE, se='E')
no.v <- names(max.nv)[max.nv==0]

int <- intersect(g, no.v)
print(setdiff(int, grep('^TRANS-RXN|POLYMER-INST-', int , value=TRUE)))
