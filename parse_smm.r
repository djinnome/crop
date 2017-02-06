##jmd
##3.7.11
##parse_smm.r

#also deal w/: acp, coa?
#charge column mostly NA: oxidized/reduced e.g. Reduced-Quinones & Cytochromes-C-Reduced, charge e.g. FE+4?
#infer chem form from: fatty-acids e.g. cis-cis-D9-27-C46-2-ACPs, rxns w/ only 1 molecule w/o chem form?

##clean metacyc chem form in smm file
clean.meta.chem.form <- function(smm){
    #clean
    smm$CHEMICAL.FORMULA <- gsub('\\(|\\)', '', smm$CHEMICAL.FORMULA)
    #extract hidden chem.forms
    no.chem.ind <- which(smm$CHEM==''|is.na(smm$CHEM))
    #add '$' to end
    smm$CHEMICAL.FORMULA[no.chem.ind] <- sub('$', '$', smm$FRAME[no.chem.ind])
    #find hidden chem forms
    hidden.chem <- intersect(no.chem.ind, grep('-INST-.+-C[[:digit:]]', smm$CHEM))
    #strip, add space, & turn '-' into '$'
    smm$CHEMICAL.FORMULA[hidden.chem] <- gsub(pattern='([[:upper:]])([[:digit:]])', replacement='\\1 \\2', 
    gsub('-', '$', gsub(pattern='.+-(C[[:digit:]])', replacement='\\1', smm$CHEM[hidden.chem])))
    return(smm)
}  

#redox
#sm[grep('(^|-)ox(idized|$|-)', sm$FRAME, ignore.case=TRUE),]

##turn meta chem form into data.frame
chemForm2df <- function(chem.form){
    chem.form.pairs <- unlist(strsplit(chem.form, split='\\$'))
    elements <- unlist(strsplit(chem.form.pairs, split=' [[:digit:]]+$'))
    count <- sub(pattern='(.+) ([[:digit:]]+)$', replacement='\\2', chem.form.pairs)
    #when chem.form='1-3-alpha-D-Glucans' ie when it has no spaces, count gives the chem.form
    normal.met.cf.ind <- grep(pattern=' [[:digit:]]+$', x=chem.form.pairs)
    #can't use count[-grep(...)]=1 since when length(chem.form)=1, 
    #grep can come up empty and assigning 1 to count[character(0)]=count[numeric(0)] doesn't work
    count[setdiff(1:length(count), normal.met.cf.ind)] <- 1
    count <- as.numeric(count)
    dat <- data.frame(els=elements, count=count)
    #if elements are not unique
    if (any(duplicated(dat$els))){
        dat.v <- tapply(X=dat$count, INDEX=dat$els, FUN=sum)
        dat <- data.frame(els=names(dat.v), count=dat.v)
        rownames(dat) <- 1:length(dat.v)
    }
    return(dat)
}
