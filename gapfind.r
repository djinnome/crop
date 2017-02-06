##jmd
##9.8.11
##gapfind.r

gapfind <- function(a, mets.oi=rownames(a)){
    a.intra <- a[-grep('\\[CCO-EXTRACELLULAR\\]$', rownames(a)),]
    root.noprod <- rownames(a.intra)[rowSums(a.intra!=0)>0 & rowSums(a.intra>0)==0]
    root.noconsume <- rownames(a.intra)[rowSums(a.intra!=0)>0 & rowSums(a.intra<0)==0]
    return(list(noprod.root=intersect(root.noprod, mets.oi), noconsume.root=intersect(root.noconsume, mets.oi)))
}

#get rid of Title Case classes & find min # meta rxns for gaplessness
#reverse BLAST
