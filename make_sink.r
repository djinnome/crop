##jmd
##1.27.12
##make_sink.r
##make file of sink rxns

##this is too simple to account for vacuolar rxns, which need a diff name, e.g. HIS[CCO-VAC-LUM]-TRANS-RXN-R2L

#s should have a column named 'rxn'
make.sink <- function(s, cpd.colname='compound', 
extracell.regexp='\\[CCO-EXTRACELLULAR\\]$', rxn.suffix='-TRANS-RXN-R2L', 
write.file='/msc/neurospora/FBA/farm_data/sink.tsv', order.col.ind=1){
    ex.mets <- unique(grep(extracell.regexp, s[,cpd.colname], value=TRUE))
    exch.df <- data.frame(rxn=paste(gsub(extracell.regexp, '', x=ex.mets), rxn.suffix, sep='') , compound=ex.mets, coeff=-1)
    exch.df <- exch.df[order(exch.df[,order.col.ind]),]
    #write
    if (!is.null(write.file)){
        write.table(exch.df, write.file, sep='\t', quote=FALSE, row.names=FALSE)
    }
    return(exch.df)
}
