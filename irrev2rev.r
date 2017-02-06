##jmd
##11.4.2011
##irrev2rev.r

irrev2rev <- function(sp, ub=matrix(Inf, nrow=1, ncol=ncol(sp), dimnames=list('1', colnames(sp)))[1,]){
  frames <- gsub('-L2R$|-R2L$', '', colnames(sp))
  rev.frames <- frames[duplicated(frames)]
  rev.rxns <- paste(rev.frames, '-L2R', sep='')
  
  sp.rev <- cBind(sp[,rev.rxns], sp[,!(frames %in% rev.frames)])
  lb <- c(-ub[rev.rxns], rep(0, ncol(sp.rev)-length(rev.rxns)))
  ub2 <- c(ub[rev.rxns], ub[!(frames %in% rev.frames)])
  names(lb) <- names(ub2) <- colnames(sp.rev)
  
  return(list(sp=sp.rev, lb=lb, ub=ub2))
}
