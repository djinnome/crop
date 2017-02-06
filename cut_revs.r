##jmd
##3.30.11
##cuts_rev.r
#constructs cutting planes to prevent island reversible rxns (cycles of length 2)
#gapless formulation: sv>=eps|s|beta is weak w/o binary constraints
#but, doesn't preclude back-to-back rev rxns = a met being in 2 cycles of length 2

#setwd('/msc/neurospora/FBA/farm')

#5.23.2011: ensure that this doesn't get bungled up by promiscuous cofactors like atp, nadh
#although i think these are already accounted for
#2.13.12: this only applies to sv>=0, but we're now using sv=0.
#2.27.12: each row only has 2 elements: one for the fwd & 1 for the reverse dir of rxn; also, all rownames are "1"

cut.revs <- function(sp, n.top.mets.rm=250, constrain.rxns=colnames(sp)){
 rxn.names <- colnames(sp)
 #don't want many combinatorial explosion of constraints for promiscuous cofactors
 rs <- rowSums(sp!=0); names(rs) <- rownames(sp)
 mets.nc <- names(rs)[rank(-rs)<n.top.mets.rm]

 #find mets that are in many rxns in nc.pwy
 #mets.nc <- rownames(s)[rowSums(s[,rxns$nc.pwy==1]>0)>=1]

 #check for rev rxns w/ mult mets not in mets.nc
 #there are >150 rev rxns w/ nnc>1
 #n.not.nc <- apply(s, MARGIN=2, FUN=function(x){ sum(!(rownames(s)[x<0] %in% mets.nc)) })
 #sum(n.not.nc>1 & rxns$rev==1)/2

 #get rev rxns
 base.rev.rxn <- gsub('-L2R|-R2L','', rxn.names)
 rev.rxns <- intersect(rxn.names[duplicated(base.rev.rxn)], constrain.rxns)
 
 ##rr is one of 2 rev rxns, chosen wlog
 #only restrict st one side needs independent 'in' source, b/c we allow export of any met
 const.mat <- Matrix(0, nrow=0, ncol=ncol(sp))
 for (rr in rev.rxns){
   #cat(rr, '/')
   #get rev rxns
   base.rr <- gsub('-L2R|-R2L','', rr)
   rxns.ind.tmp <- which(base.rev.rxn %in% base.rr)
   #ignore mets.nc, since these have beta's=1, so won't add any information
   mets.in <- setdiff(rownames(sp)[sp[,rr]<0], mets.nc)
   mets.out <- setdiff(rownames(sp)[sp[,rr]>0], mets.nc)
   n.eqns <- length(mets.in)*length(mets.out)
   if (n.eqns>0){
     #set-up submatrix to store constraints
     submat.tmp <- Matrix(0, nrow=n.eqns, ncol=ncol(sp), dimnames=list(paste(base.rr, 1:n.eqns, sep='_'), colnames(sp)))
     submat.tmp[,rxns.ind.tmp] <- 1/2
     #make matrix of pairs of in & out
     pairs.mat <- as.matrix(expand.grid(mets.in, mets.out, stringsAsFactors=FALSE))
     #for each pair of mets, find independent rxns which produce at least one of them
     for (pm.row in 1:nrow(pairs.mat)){
       mets.tmp <- pairs.mat[pm.row,]
       prod.rxns.tmp <- which(colSums(sp[mets.tmp,]>0)>0)
       prod.rxns.tmp <- setdiff(prod.rxns.tmp, rxns.ind.tmp)
       #set -1 for rxns where a met is a product
       if (length(prod.rxns.tmp)>0){ submat.tmp[pm.row, prod.rxns.tmp] <- -1 } #does this line slow it down??
     }#end for pm.row
     const.mat <- rBind(const.mat, submat.tmp)
   }#end if n.eqns>0
 }#end for rr
 return(const.mat)
}#end fcn
