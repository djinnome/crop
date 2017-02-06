##jmd
##10.26.2011
##xylose_story.r

source('/msc/neurospora/FBA/farm/farm_header.r')

mm.colnames <- c('xylose-in', 'xylitol-out')
mm <- Matrix(0, nrow=nrow(s.al), ncol=length(mm.colnames))
dimnames(mm) <- list(rownames(s.al), mm.colnames)
mm['XYLOSE[CCO-EXTRACELLULAR]', 1] <- 1 #fba model won't secrete xylitol, b/c i think it's secreted by diffusion
mm['XYLITOL', 2] <- -1

sp.al <- s.sp[rownames(s.al), colnames(s.al)]
s2 <- cBind(s.al, mm)
#s2 <- cBind(s.sp[rownames(s.al), colnames(s.al)], mm)

ub <- rep(1000, ncol(s2)); names(ub) <- colnames(s2)

#ub["RXN-3341-L2R"] <- 0
#ub[c('XYLISOM-RXN-L2R', 'XYLISOM-RXN-R2L')] <- 0
#ub['L-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN-R2L'] <- 0
#ub[c('DALADALALIG-RXN-R2L', 'FORMATETHFLIG-RXN-R2L')] <- 0

#########################################################################################################################
##FBA for o2 vs atp
#########################################################################################################################
ub['OXYGEN-MOLECULE-TRANS-RXN-L2R'] <- 100
ff <- FBA(s2, fba.ub=ub, goal='atp-out')
print.flux(v=ff$xopt)

#########################################################################################################################
##xylose 2 etoh
#########################################################################################################################
##regulation of NADPH production
#ssm <- ss.mat(s.al, r='NADPH')
#ub[setdiff(names(ssm)[ssm>0], 'GLU6PDEHYDROG-RXN-L2R')] <- 0

#fbawmc=500; k=2; suc>=5
#(mat0 <- xylose2etoh(s2=s2, ub=ub, suc.ub=5, k=2, fbawmc=500))
o.thresh.v <- c(0, 10**(-9:2))
mat0 <- xylose2etoh(s2=s2, ub=ub, xr.nad.ub=1000, xylose.in.v=rep(10, length(o.thresh.v)), o.thresh.v=o.thresh.v, fba.obj='biomass'); 
mat0$out

mat0 <- xylose2etoh(s2=s2, ub=ub, xr.nad.ub=1000, k=100, fba.obj='ETOH-TRANS-RXN-R2L'); mat0$out

i=3
#pv=print.flux(v=mat0$flux[,i], out='xylose_flux_lowO2.tsv')
v <- mat0$flux[,i]
fom <- flux.of.met(a=s2, v=v, 'ATP')
fom <- fom[fom>0]
sort(fom)

#frames <- gsub('\\[CCO-.+|-L2R|-R2L', '', rownames(mat0$flux))
#v.df <- data.frame(frames, rxn.names=rownames(mat0$flux), mat0$flux)
#write.table(v.df, 'xylose_fluxes.tsv', row.names=FALSE, quote=FALSE, sep='\t')

lst <- list()
for (fbawmc in 100*seq(from=5, to=10, by=1)){
    for (k in c(1, 1.5, 2, 2.5, 3)){
        for (suc in 2:4){
            print(name <- paste('suc:', suc, ';fbawmc:', fbawmc, ';k:', k, sep=''))
            print(lst[[name]] <- xylose2etoh(s2=s2, ub=ub, suc.ub=suc, k=k, fbawmc=fbawmc)$output)
            print('')
        }
    }
}

filt.ind <- which(sapply(lst, FUN=function(x){
    all(c('0', '12.6') %in% rownames(x)[1:2], x['8.4', 'etoh.yield']>3*x['12.6', 'etoh.yield'], x['8.4', 'etoh.yield']>0.15, na.rm=TRUE)
}))
lst[filt.ind]

#xyl=20 & o2=0.75 -> etoh=29

##plot
#ouch!
mm <- c(0.0000, 0.2330, 0.2180, 0.1590, 0.0724)
names(mm) <- c("0", "3.3", "5.1", "8.4", "12.6")
m2 <- c(0.066, 0.25, 0.28, 0.34, 0)
names(m2) <- names(mm)
plot(x=as.numeric(names(m2)), y=m2, xlab='o2', ylab='etoh yield', type='l')
lines(x=as.numeric(names(mm)), y=mm, col='red')

##on pc
setwd('Z:/reconstruct')
x2e <- as.matrix(read.csv('xyl2etoh.csv'))
colnames(x2e) <- c('o2', 'etoh')
