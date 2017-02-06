##jmd
##1.20.12
##metab_dilute_s.r
##prepare S for linear metabolite dilution FBA (LIMED-FBA)

#compartments represents suffixes
#by default, don't dilute extracellular cpds, since don't want to dilute trans rxns. can't prevent dilution of other trans rxns since eps acts on cpds not rxns.
#can make limed act on rxns by making an altered S on the rhs where column is all 0 if rxn shouldn't be diluted.
metab.dilute.s <- function(sp, eps=10**-4, nondilute.cpds=c('PROTON', 'WATER', 'OH', 'OXYGEN-MOLECULE', 'CARBON-DIOXIDE'),
nondilute.rxns=NULL, 
compartments.dilute=c('', paste('[CCO-', c('MIT', 'GLYOXYSOME', 'NUC-LUM', 'VAC-LUM'), ']', sep='')) ){
    #find indices of non-zero stoichiometric elements to dilute
    w <- which(sp!=0, arr.ind=TRUE)
    colnames(w) <- c('row', 'col')
    
    #get cpds to not dilute in any compartment
    #'' in rep is for cytosolic cpds
    nondilute.cpds.v <-  paste(nondilute.cpds, rep(compartments.dilute, each=length(nondilute.cpds)), sep='')
    nondilute.rows <- which(rownames(sp) %in% nondilute.cpds.v)
    
    #get rxns to not dilute
    nondilute.cols <- which(colnames(sp) %in% nondilute.rxns)
    
    #set elements to dilute
    w <- w[!(w[,'row'] %in% nondilute.rows) & !(w[,'col'] %in% nondilute.cols),]
    
    #get sum of rxns diluted per metabolite
    rs <- rowSums(sp[,setdiff(colnames(sp), nondilute.rxns)]!=0)
    names(rs) <- rownames(sp)
    #metab conc go up to 10^-2 (bennet et al, nat chem bio, 2009) and we use ub=10^3, so multiply by 10^-5
    #this sometimes gives numerical error, though, so use 10^-3 and verify KOs using FBA: limed-fba is tool to delete 'island cycles' and find unproducible intermediate mets, not for quant growth prediction.
    sp[w] <- sp[w]-eps/rs[w[,'row']]
    
    return(sp)
}
