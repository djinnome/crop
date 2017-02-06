##jmd
##10.15.12
##sl_fromMlab_heat.r

# expect sl.mat to be symmetric, w/ rownames=colnames
getSLpartners <- function(sl.mat, genes){
    stopifnot(genes %in% rownames(sl.mat))
    submat <- sl.mat[genes,]
    partners <- colnames(submat)[colSums(submat>0)>0]
    return(partners)
}

source('farm_header.r')

##read
sl <- read.table('matlab_data/nonIsozymeSynthLethals.txt', as.is=TRUE, header=TRUE)
pwys.sl <- read.delim('matlab_data/nonIsozymeSynthLethalsPwys.tsv', as.is=TRUE)
# check that they're same, except for 1st 6 rows, which hold col-2 and bal
stopifnot(all(sl[,1:2]==pwys.sl[,1:2]))
# only need common pwys column of pwys.sl
sl <- pwys.sl[,1:3]

##add pairs
sl <- rbind(sl, c('NCU06532.5', 'NCU07334.5', ''), c('NCU03235.5', 'NCU04433.5', 'Sulfate Transport Rxn'))

## sl names
sl.names <- sort(unique(c(sl[,1], sl[,2])))
names(sl.names) <- sl.names
sl.names[intersect(rownames(ncu.gene.map), names(sl.names))] <- ncu.gene.map[intersect(rownames(ncu.gene.map), names(sl.names)), 'Symbols']
# create 2 more columns, with symbol if available, else ncu
sl <- cbind(name1=sl.names[sl[,1]], name2=sl.names[sl[,2]], sl)
rownames(sl) <- paste(sl[,1], ',', sl[,2], sep='')

### get mat #############################################################################################################
##make sl.mat for plotting
sl.mat <- matrix(0, nrow=length(sl.names), ncol=length(sl.names), dimnames=list(sl.names, sl.names))
#1: common; 2: diff
for (i in 1:nrow(sl.mat)){
    for (j in 1:ncol(sl.mat)){
        potential.rnames <- c(paste(sl.names[i], ',', sl.names[j], sep=''), paste(sl.names[j], ',', sl.names[i], sep=''))
        if (any(potential.rnames %in% rownames(sl))){
            if (sl[rownames(sl) %in% potential.rnames, 'Common.pwys']!=''){
                sl.mat[i,j] <- 1
            } else {
                sl.mat[i,j] <- 2
            }
        }
    }#end if j
}#end if i

## collapse
cmplx1.genes <- sl.names[intersect(gpr$Genes[gpr$Enz=='CPLX2T-21'], names(sl.names))]
cmplx3.genes <- sl.names[intersect(gpr$Genes[gpr$Enz=='CPLX2T-46'], names(sl.names))]
cmplx4.genes <- sl.names[intersect(gpr$Genes[gpr$Enz=='CPLX2T-51'], names(sl.names))]
cmplx.genes <- c(cmplx1.genes, cmplx3.genes, cmplx4.genes)
#cmplx 1 & 3 interact w/ suc
stopifnot(getSLpartners(sl.mat, c(cmplx1.genes, cmplx3.genes))=='suc')
# cmplx 4 should have same interactions for each of its subunits, these partners are rownames of cmplx4.submat
cmplx4.submat <- ss.mat(sl.mat, cmplx4.genes)
stopifnot(all(cmplx4.submat==2))

#re-build sl.mat
sl2 <- sl.mat[setdiff(rownames(sl.mat), cmplx.genes), setdiff(colnames(sl.mat), cmplx.genes)]
sl2 <- rbind(0, rbind(0, rbind(0, sl2)))
sl2 <- cbind(0, cbind(0, cbind(0, sl2)))
sl2[1:2, 'suc'] <- sl2['suc', 1:2] <- 2
sl2[3, rownames(cmplx4.submat)] <- sl2[rownames(cmplx4.submat), 3] <- 2
#re-name
rownames(sl2)[1:3] <- colnames(sl2)[1:3] <- paste('Complex', c('I','III','IV'))
sl2 <- sl2[order(rownames(sl2)), order(colnames(sl2))]

## set diag to -1
diag(sl2) <- -1

###plot##################################################################################################################
pdf('docs/Figures/SyntheticLethals/sl_heat_mlab.pdf')
hc <- hclust(dist(t(sl2)))
z <- sl2[hc$order, hc$order]
make.grid(z, color.v=c('gray', 'white', 'cyan', 'orange'), cex.axis=0.3, lwd=1)

z2 <- matrix(0, nrow=nrow(z), ncol=ncol(z), dimnames=dimnames(z))
z2[rbind(c('am', "en(am)-2"), c("en(am)-2", 'am'), c('pro-3', 'ota'), c('ota', 'pro-3'), 
c('pyr-1', 'uc-5'), c('uc-5', 'pyr-1'), c('cys-14', 'cys-13'), c('cys-13', 'cys-14'))] <- 1
my.rectangle(z2, equal.val=1, border='black', col=NULL, clear.prev.sq=FALSE, density=0)
dev.off()
