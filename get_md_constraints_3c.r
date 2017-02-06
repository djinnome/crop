##jmd
##2.24.12
##get_md_constraints_3c.r

#for accum, fill in rows w/ -s.intra*accum_v-eye*accum_slack_v+s.intra.norm*beta<=0
#for deplete, fill in rows w/ S.intra*accum_v-eye*accum_slack_v+s.intra.norm*beta<=0
#for both accum & deplete, fill in rows for: eye*accum_v-diagonal(md.f.ub)<=0
#dim(sp.intra)==dim(sp.intra.norm) & colnames(sp)==colnames(sp.intra)
get.md.constraints.3c <- function(sp.intra.3c, sp.intra.norm.3c, rnames.sp.intra, cnames.sp.intra, md.f.ub){
    cpds.accum.tmp <- paste('accum_', sp.intra.3c$cpd, sep='')
    cpds.deplete.tmp <- paste('deplete_', sp.intra.3c$cpd, sep='')
    
    nrxns <- length(cnames.sp.intra)
    nmets.intra <- length(rnames.sp.intra)
    
    md.3c <- rbind(
    data.frame(rxn=paste('accum_', sp.intra.3c$rxn, sep=''), cpd=cpds.accum.tmp, coeff=-sp.intra.3c$coeff),
    data.frame(rxn=paste('accum_slack_', rnames.sp.intra, sep=''), cpd=paste('accum_', rnames.sp.intra, sep=''), coeff=-1),
    data.frame(rxn=paste('beta_', sp.intra.3c$rxn, sep=''), cpd=cpds.accum.tmp, coeff=sp.intra.norm.3c$coeff),

    data.frame(rxn=paste('deplete_', sp.intra.3c$rxn, sep=''), cpd=cpds.deplete.tmp, coeff=sp.intra.3c$coeff),
    data.frame(rxn=paste('deplete_slack_', rnames.sp.intra, sep=''), cpd=paste('deplete_', rnames.sp.intra, sep=''), coeff=-1),
    data.frame(rxn=paste('beta_', sp.intra.3c$rxn, sep=''), cpd=cpds.deplete.tmp, coeff=sp.intra.norm.3c$coeff),
    
    data.frame(rxn=paste('accum_', cnames.sp.intra, sep=''), cpd=paste('accum_', cnames.sp.intra, sep=''), coeff=1),
    data.frame(rxn=paste('beta_', cnames.sp.intra, sep=''), cpd=paste('accum_', cnames.sp.intra, sep=''), coeff=-md.f.ub),
    
    data.frame(rxn=paste('deplete_', cnames.sp.intra, sep=''), cpd=paste('deplete_', cnames.sp.intra, sep=''), coeff=1),
    data.frame(rxn=paste('beta_', cnames.sp.intra, sep=''), cpd=paste('deplete_', cnames.sp.intra, sep=''), coeff=-md.f.ub)
    )
    
    md.cnames <- paste(rep(c('accum_', 'deplete_', 'accum_slack_', 'deplete_slack_'), times=c(nrxns, nrxns, nmets.intra, nmets.intra)),
    c(cnames.sp.intra, cnames.sp.intra, rnames.sp.intra, rnames.sp.intra), sep='')
    
    md.rnames <- c(paste(rep(c('accum_', 'deplete_'), times=c(nmets.intra, nmets.intra)), c(rnames.sp.intra, rnames.sp.intra), sep=''),
    paste(rep(c('accum_', 'deplete_'), times=c(nrxns, nrxns)), c(cnames.sp.intra, cnames.sp.intra), sep=''))

    return(list(amat.3c=md.3c, cnames=md.cnames, rnames=md.rnames))
}
