##jmd
##7.15.11
##print_fluxes.r

#names of v are rxn names
print.flux <- function(v, ec.file='/msc/neurospora/FBA/farm_data/Neurospora/nc10.EC', out.file=NULL){
    ec <- read.delim(ec.file)
    flux.mat <- data.frame(Frame=NA, flux=v, rxn=names(v))
    #take non-zero fluxes
    flux.mat <- flux.mat[flux.mat$flux>0,]
    #mult R2l by -1
    flux.mat$flux[grep('-R2L$', flux.mat$rxn)] <- -flux.mat$flux[grep('-R2L$', flux.mat$rxn)]
    #match frames
    flux.mat$Frame <- ec$FRAME.ID[match(flux.mat$rxn, ec$rxn)]
    #order by frame
    flux.mat <- flux.mat[order(flux.mat$Frame),]
    #write
    if (!is.null(out.file)){ write.table(x=flux.mat, file=out.file, sep='\t', quote=FALSE, row.names=FALSE) }
    #return
    return(flux.mat)
}
