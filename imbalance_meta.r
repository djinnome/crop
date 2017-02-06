##jmd
##7.21.11
##imbalance_meta.r
##get met.mat & rxns that generate mass
#representing mets as, eg (C 13)$(H 20)$(O 3)$, is annoying but clear, b/c many 'elements' are eg "GLY-tRNAs" or "Fructans"

library('CHNOSZ')

mm.ss <- get.met.mat(univ.smm=smm.meta, org.smm=smm.nc, s.rownames=rownames(s.sp), chem.form.col='CHEMICAL.FORMULA', colsums.mm.thresh=2)
write.csv(mm.ss, 'met_mat.csv')
#mm.ss <- as.matrix(read.csv('met_mat.csv', row.names=1))

##compute imbalance
#some mets are empty in all of CHNOPS
bal.mat <- as.matrix(t(t(mm.ss) %*% s.sp))
trans.rxns <- grep('-TRANS-RXN', colnames(s.sp), value=TRUE)
bal.mat[intersect(rownames(bal.mat), c(trans.rxns, bm$rxn)),] <- 0
bal.mat <- bal.mat[,colSums(abs(bal.mat))>0]
#need to write this for make_rxn_annot
write.csv(bal.mat, 'rxns_imbalance.csv')
#bal.mat <- read.csv('rxns_imbalance.csv', row.names=1)
#bal.mat <- bal.mat[,c('C', 'H', 'O', 'N', 'P', 'S')]

###check#########################################################################################
##check bal.mat
#rxns
#rs <- rowSums(abs(bal.mat)); names(rs) <- colnames(s.sp)
#mean(rs.nz <- rs!=0 & !(names(rs) %in% trans.rxns))
#rs[rs.nz][1:5]
#mets
#cs <- colSums(bal.mat!=0); names(cs) <- u.els
#mean(cs!=0) #cs[cs!=0]

##check a rxn
#rxn <- 'CIT-MET-RXN'
#(ss <- s.sp[s.sp[,rxn]!=0, rxn])
#(mm.ss <- met.mat[names(ss), colSums(abs(met.mat[names(ss),]))!=0])
#(bm.ss <- bal.mat[rxn, colnames(mm.ss)])

##check a met
#met <- 'POLYMER-INST-Peptidoglycans-C312-H512-N64-O152'
#met.mat[met,met.mat[met,]!=0]
#in.rxns <- which(s.sp[met,]!=0)
#bal.mat[in.rxns, bal.mat[in.rxns,]!=0]

##check an element
#elem <- 'O'
#bal.mat[bal.mat[,elem]!=0,1:10]
