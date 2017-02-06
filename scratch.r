##jmd
##10.31.11
##scratch.r

source('/msc/neurospora/FBA/farm/farm_header.r')
goals <- unique(bm$compound[bm$coeff<0])

####compare SL lists#####################################################################################################
setwd('C:/Documents and Settings/jdreyf/Desktop/farm_data/matlab_data')

sl1 <- read.table('synLethals_wNames.txt', header=TRUE, as.is=TRUE)
sl2 <- read.table('synLethals2.txt', header=TRUE, as.is=TRUE)

sl1.v <- paste(sl1[,3], sl1[,4], sep=';')
sl2.v <- paste(sl2[,1], sl2[,2], sep=';')

sl2.v <- sl2.v[!(sl2.v %in% sl1.v)]
sort(summary(as.factor(unlist(strsplit(sl2.v, split=';')))))

###growth data###########################################################################################################
#s.al['SUCROSE[CCO-EXTRACELLULAR]', 'SUCROSE-TRANS-RXN-L2R'] <- 0
#s.al['GLC[CCO-EXTRACELLULAR]', 'SUCROSE-TRANS-RXN-L2R'] <- 1
fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 3 #2.88 #0.81

n.atp <- 33
s.al['ATP', 'biomass'] <- -n.atp
s.al['ADP', 'biomass'] <- s.al['Pi', 'biomass'] <- n.atp

f <- FBA(s.al, fba.ub=fva.ub)

###oneprune##############################################################################################################
# nov 19, 2012
# current implementation of max.n.flux excludes purely extracellular rxns, eg invertase
v <- max.n.flux(s.al)
mr2 <- setdiff(model.rxns, names(v)[v==0])
s2 <- ss.mat(sg, c=mr2)
f.full <- FBA(s2, sense='E', fba.ub=fva.ub[colnames(s2)], fba.obj='FullBiomassComposition')

write.csv(mr2, 'model_rxns.csv', row.names=FALSE)

###supp-by-genes#########################################################################################################
setwd('C:/Documents and Settings/jdreyf/Desktop/farm_data')
smm <- read.delim('Neurospora/nc10.smm', na='NIL')
r.dat <- as.matrix(read.delim('docs/Figures/EssentialsAndSupplements/train_test_supps_x_genes.tsv'))
w <- which(matrix(r.dat %in% c(2:3, 5:6), nrow=nrow(r.dat)), arr.ind=TRUE)
#rm .R.. for pantothenate
met.names <- gsub('.', '-', fixed=TRUE, sub('.R..', '(R)-', fixed=TRUE, sub('^X', '', colnames(r.dat))))
#get frame name to match matlab output
met.frames <- gsub('\\[.+\\]$', '', smm[match(toupper(met.names), toupper(smm[,4])),1])
met.frames[is.na(met.frames)] <- c('MET and THR', 'PHE and TYR') #these 2 don't match, & they show up in this order
#output
pairs.mat <- cbind(rownames(r.dat)[w[,'row']], met.frames[w[,'col']])
write.table(pairs.mat, row.names=FALSE, quote=FALSE, sep='\t')

###fba accuracy##########################################################################################################
#model.rxns <- c(model.rxns, "RXN-7609-L2R", "TRANS-RXN2T-233-R2L[CCO-EXTRACELLULAR]") 
#"CARNITINE-O-ACETYLTRANSFERASE-RXN-L2R", "CARNITINE-O-ACETYLTRANSFERASE-RXN-R2L", )
s.ck <- ss.mat(sg, model.rxns)
sp.ck <- s.sp[rownames(s.ck), colnames(s.ck)]
sp.ck[c('NADH', 'FMN'), 'biomass'] <- -0.01
#add.trans <- Matrix(0, nrow=nrow(sp.ck), ncol=2, dimnames=list(rownames(sp.ck), c('carn-trans-rxn-L2R', 'carn-trans-rxn-R2L')))
#add.trans[c('CARNITINE', 'O-ACETYLCARNITINE'),] <- matrix(c(1,-1,-1,1), ncol=2)
#sp.ck <- cBind(sp.ck, add.trans)
s.ck <- sp.ck

#acu-5 (NCU06836.5)
ko.rxns <- NCUvector2rxns("NCU06836.5", gpr=gpr)
f <- FBA(s.ck, ko.rxns=ko.rxns)
w <- named.vec(1, name=rownames(s.ck))
w[c('ACET[CCO-MIT]')] <- 100
me <- min.export(s.ck[,setdiff(colnames(s.ck), ko.rxns)], w=w)

#ad-5
ko.rxns <- NCUvector2rxns("NCU06187.5", gpr=gpr)
s2 <- changeNuts(s.ck, add='ADENINE')
f <- FBA(s2, ko.rxns=ko.rxns, se="E")
w <- named.vec(1, name=rownames(s.ck))
w[c('AICAR', 'PHOSPHORIBOSYL-FORMAMIDO-CARBOXAMIDE', 'IMP[CCO-GLYOXYSOME]')] <- 100
w['GUANOSINE'] <- 0
me <- min.export(s2[,setdiff(colnames(s2), ko.rxns)], w=w)

#ad-8
ko.rxns <- NCUvector2rxns("NCU09789.5", gpr=gpr)
s2 <- changeNuts(s.ck, add='ADENINE')
f <- FBA(s2, ko.rxns=ko.rxns, se="E")
me <- min.export(s2[,setdiff(colnames(s2), ko.rxns)])

#pan-2 grows :(
ko.rxns <- NCUvector2rxns("NCU10048.5", gpr=gpr)
f <- FBA(s.ck, ko.rxns=ko.rxns, se="E", fba.ub=fva.ub)

###dups##################################################################################################################
rxns.v <- apply(sp.al, 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x, sep=':', collapse=';') })
rxns.dup <- rxns.v[duplicated(rxns.v)]
sort(rxns.v[rxns.v %in% rxns.dup])

del <- c('3.6.3.16-RXN-L2R[CCO-PM-FUNGI]', 'TYRAMINOTRANS-RXN-L2R', 'RXN2T-102-L2R[CCO-GLYOXYSOME]', 
'TRANS-RXN2T-18-GLC/PROTON//GLC/PROTON.23.-L2R[CCO-PM-FUNGI]', 
'TRANS-RXN2T-234-L-ALPHA-ALANINE/PROTON//L-ALPHA-ALANINE/PROTON.47.-L2R[CCO-EXTRACELLULAR]', 
'TRANS-RXN2T-68-L2R[CCO-PM-FUNGI]', 'TRANS-RXN-195-L2R[CCO-EXTRACELLULAR]', 
'TRANS-RXN2T-18-XYLOSE/PROTON//XYLOSE/PROTON.29.-L2R[CCO-PM-FUNGI]')

3.6.3.16-RXN-L2R[CCO-PM-FUNGI], TYRAMINOTRANS-RXN-L2R, RXN2T-102-L2R[CCO-GLYOXYSOME], 
TRANS-RXN2T-18-GLC/PROTON//GLC/PROTON.23.-L2R[CCO-PM-FUNGI], 
TRANS-RXN2T-234-L-ALPHA-ALANINE/PROTON//L-ALPHA-ALANINE/PROTON.47.-L2R[CCO-EXTRACELLULAR], 
TRANS-RXN2T-68-L2R[CCO-PM-FUNGI], TRANS-RXN-195-L2R[CCO-EXTRACELLULAR], 
TRANS-RXN2T-18-XYLOSE/PROTON//XYLOSE/PROTON.29.-L2R[CCO-PM-FUNGI]

###add transports########################################################################################################
trans <- read.csv('transport_added_22may12.csv')
#ID transport rxns
add.trans <- NULL
for (i in 1:nrow(trans)){
    add.trans <- c(add.trans, colnames(s.sp)[apply(s.sp[c(trans[i,1], paste(trans[i,1], '[CCO-EXTRACELLULAR]', sep='')),], 
    MARGIN=2, FUN=function(x){ all(x==c(1,-1)) })])
}
#add these, after deleting unwanted dups, to model.rxns
model.rxns <- c(model.rxns, trans$rxn)
write.csv(model.rxns, 'model_rxns.csv', row.names=FALSE)

###gene-by-supps##########################################################################################################
#ace-2,3,4
s2 <- changeNuts(sp.al, add.nuts=c('ACET', 'GLN'), rm.nuts='SUCROSE')
fva.ub <- named.vec(1000, name=colnames(s2)); fva.ub[paste(c('ACET', 'GLN'), '-TRANS-RXN-L2R', sep='')] <- 50
ff <- FBA(s2, ko="PYRUVDEH-RXN-L2R[CCO-MIT]", control=list(trace=0), fba.ub=fva.ub, min2=TRUE)
mit.rxns <- grep('CCO-MIT', names(ff$x)[ff$x>0], value=TRUE)
ss.mat(s2, mit.rxns)
#citrate not used, but suc & fum are imported thru TRANS-RXN2T-82/83


#acu-3
nn <- NCUvector2rxns(ncu='NCU04230.5', gpr=gpr)

s2a <- changeNuts(s.al, add.nuts=c('ADENINE', 'ACET'), rm.nuts='SUCROSE')
ff <- FBA(s2a, ko=nn, control=list(trace=0))

s2l <- changeNuts(s.al, add.nuts=c('LYS', 'ACET'), rm.nuts='SUCROSE')
wts <- named.vec(data=1, name=rownames(sg))
wts[rownames(s2l)] <- 0
m <- mma.bin(sg, wts=wts) #nothing
ff <- FBA(s2l, ko=nn, control=list(trace=0), goal=goals, ub.goals=0.01)

s2a <- changeNuts(s.al, add.nuts=c('ADENINE'))

##acu-5
nn <- NCUvector2rxns(ncu='NCU07659.5', gpr=gpr)
s2l <- changeNuts(sp.al, add.nuts=c('ACET'))
ff <- FBA(s2l, ko=nn, control=list(trace=0))

###biolog transport######################################################################################################
#get cpds for transport curation
pp=t(apply(pred, 1, FUN=function(x){
    y=unlist(strsplit(x=x['cpd'], split=','))
    as.numeric(c(any(y %in% rownames(s.al)),  any(paste(y, '[CCO-EXTRACELLULAR]', sep='') %in% rownames(s.al))
    )) 
}))
pp=cbind(pp, pred$pred)
colnames(pp) <- c('intra', 'extra', 'growth')

pp=pp[!duplicated(pred$cpd),]
pp=pp[rowSums(pp[,1:2])>0,]
pp=pp[pp[,2]==0,]
rownames(pp) <- gsub('(PM...)(.)$', replacement='\\10\\2', rownames(pp))

#look at non-imported extracellular biolog cpds: these are sugar polymers or di-amino acids
extra.cpds <- pb[rownames(pp)[pp[,2]==1], 'cpd']
e2 <- extra.cpds[!gsub('\\[CCO-.+', '', extra.cpds) %in% rownames(s.al)]
for (i in 1:length(e2)){ print(ss.rom(s.al, e2[i])) }

###arg-14################################################################################################################
f2 <- FBA(sp.al, sense='E', fba.ub=fva.ub, ko='N-ACETYLTRANSFER-RXN-L2R')
f2$x[c('GLUTAMATE-N-ACETYLTRANSFERASE-RXN-L2R', 'ACETYLGLUTKIN-RXN-L2R', 'N-ACETYLGLUTPREDUCT-RXN-R2L', 'ACETYLORNTRANSAM-RXN-R2L',
'ACETYLORNDEACET-RXN-L2R')]

mets <- c('ACETYL-GLU', 'N-ACETYL-GLUTAMYL-P', 'CPD-469', 'N-ALPHA-ACETYLORNITHINE')
for (i in 1:length(mets)){ 
    fom <- flux.of.met(sp.al, f2$x, mets[i])
    print(fom[order(fom)])
}

###make gene-by-supp heatmap#############################################################################################
setwd("C:/Documents and Settings/jdreyf/Desktop/farm_data")
sg <- as.matrix(read.delim('docs/Figures/EssentialsAndSupplements/train_test_supps_x_genes.tsv', row.names=1, as.is=TRUE))
z <- t(sg)

pdf('docs/Figures/train_test_supps_x_genes.pdf')
make.grid(z)

my.rec(equ=4, col='red')
my.rec(equ=5, col='blue')
my.rec(equ=6, col='purple')
dev.off()

###ox phos + tca#########################################################################################################
fva.ub <- rep(10**3, nrxns); names(fva.ub) <- colnames(s.al)
fva.ub['SUCROSE-TRANS-RXN-L2R'] <- 5
fva.ub['MAL-DEH-GLYOX-RXN-L2R'] <- 0
#fva.ub[grep('[CCO-GLYOXYSOME]', colnames(s.al), v=T, fixed=T)] <- 0
fva.ub['PYRUVATE-CARBOXYLASE-RXN-L2R[CCO-MIT]'] <- 0
#fva.ub['RXN2T-453-L2R[CCO-MIT]'] <- 0
se <- rep('E', nrow(s.al)); names(se) <- rownames(s.al); se['ATP'] <- 'G'
f <- FBA(s.al, fba.ub=fva.ub, se=se)

w=rep(1, nrow(s.al)); names(w) <- rownames(s.al); w[c('NADH[CCO-MIT]', 'PYRUVATE[CCO-MIT]')] <- 1000
me <- min.export(s.al, ub=fva.ub, v.min=0.01, w=w)

fg <- FBA(s.al, goal=goals, fba.ub=fva.ub, ub.goals=0.01)

rr <- nc.pwy$rxn[nc.pwy$PathwayID %in% c('PWY2T-2574', 'PWY2T-4302', 'PWY2T-3781', 'PWY-5690')]
r2 <- unlist(apply(as.matrix(rr), 1, FUN=function(x) grep(x, colnames(s.al), v=T)))
xx <- f$xopt[r2]
xx[order(-xx)]

r3 <- setdiff(grep('[CCO-MIT]', r2, v=T, f=T), 'FUMHYDR-RXN-L2R[CCO-MIT]')
f <- FBA(s.al, fba.ub=fva.ub, se='G', goal=)

flux.of.met(s.al, v=f$x, 'NADH[CCO-MIT]')

y <- numeric(nrow(s.sp)); names(y) <- rownames(s.sp)
y['ATP'] <- -1
wh <- which(apply(s.sp, 2, FUN=function(x) all(x==y)))

##wrong test supps#######################################################################################################
#ad-5
nn <- NCUvector2rxns(ncu='NCU02629.5', gpr=gpr)[,1]
s2a <- changeNuts(s.al, add.nuts=c('ADENINE'))
s2h <- changeNuts(s.al, add.nuts=c('HYPOXANTHINE'))

s2 <- s2h
se <- rep('E', nrow(s2)); names(se) <- rownames(s2); se['AICAR'] <- 'G'
ff <- FBA(s2, ko=nn, control=list(trace=0), se=se)

#ser-2
nn <- NCUvector2rxns(ncu='NCU01439.5', gpr=gpr)[1]
s2 <- changeNuts(s.al, add.nuts=c('SER'))
ff <- FBA(s2, ko=nn, control=list(trace=0), goal='SER')

#spe-2
nn <- NCUvector2rxns(ncu='NCU01083.5', gpr=gpr)[1]
s2 <- changeNuts(s.al, add.nuts=c('SPERMINE'))
ff <- FBA(s2, ko=nn, control=list(trace=0))

###ace-8#################################################################################################################
trace.mets(sp=s.al[,setdiff(colnames(s.al), c("PEPDEPHOS-RXN-R2L", '2.5.1.19-RXN-L2R', 'ANTHRANSYN-RXN-L2R'))], from='PHOSPHO-ENOL-PYRUVATE', to='PYRUVATE', pwy=nc.pwy)

###wrong test non-essentials#############################################################################################
d4 <- d3[d3$growth<10**-6,]
d4.rxns <- NCUvector2rxns(ncu=rownames(d4), gpr=gpr)
goals <- unique(bm$compound[bm$coeff<0])

#vanilla fba
sp.al <- s.sp[rownames(s.al), colnames(s.al)]
d4.ck <- check.ko(s.test=sp.al, ko.lst=d4.rxns, ub=fva.ub, sense='E', obs=1, quiet.fba=TRUE, ctrl=list(trace=0, method=1))$pred
names(d4.ck) <- names(d4.rxns)

#NCU00742 
f1 <- FBA(s.al, ko=d4.rxns[[1]])
f1b <- FBA(s.al, ko=d4.rxns[[1]][1])

#lys-7
f4 <- FBA(s.al, ko=d4.rxns[[4]], se='G')
f4g <- FBA(sp.al, ko=d4.rxns[[4]], goal=goals, quiet=FALSE, control=list(trace=0), se='G')

#ndk-1
f5g <- FBA(sp.al, ko=d4.rxns[[5]], goal=goals, quiet=FALSE, control=list(trace=0))

#pyr-2
f9g <- FBA(sp.al, ko=d4.rxns[[9]], goal=goals, quiet=FALSE, control=list(trace=0))

#NCU06348
f10 <- FBA(sp.al, ko=d4.rxns[[10]], se='G', goal=goals, quiet=FALSE, control=list(trace=0))

#pyr-1
f11 <- FBA(s.al, ko=d4.rxns[[11]][1:3], quiet=FALSE, control=list(trace=0))

#mtr
f12 <- FBA(s.al, ko=d4.rxns[[12]][-4], quiet=FALSE, control=list(trace=0))

#pab-1: grows on fba but not limed-fba
f13 <- FBA(sp.al, ko=d4.rxns[[13]], quiet=FALSE, control=list(trace=0))
#can't make tetrahydrofolate. 
#instead, recycles THF w/ 5,10-methylene-THF in PWY-2201, then turns tetrahydropteroyltri-L-glutamate (CPD-1301)
#from cys-methionine pwy (PWY2T-24) to tetrahydrofolate.
flux.of.met(a=sp.al, v=f13$xopt, met='CPD-1301')

#fks / do
f14 <- FBA(s.al, ko=d4.rxns[[14]])
f14g <- FBA(sp.al, ko=d4.rxns[[14]], goal=goals, quiet=FALSE, control=list(trace=0))

#pan-3
f16 <- FBA(s.al, ko=d4.rxns[[16]])
f16g <- FBA(sp.al, ko=d4.rxns[[16]], goal='PANTOTHENATE', quiet=FALSE, control=list(trace=0))

#mitochondrial carrier protein LEU5
f17 <- FBA(s.al, ko=d4.rxns[[17]])
f17g <- FBA(sp.al, ko=d4.rxns[[17]][1:2], goal='CO-A[CCO-MIT]', quiet=FALSE, se='G', control=list(trace=0))

#amine oxidase
f18 <- FBA(s.al, ko=d4.rxns[[18]])
f18g <- FBA(sp.al, ko=d4.rxns[[18]], goal='SPERMIDINE', quiet=FALSE, se='G', control=list(trace=0))

#rg-1
f19 <- FBA(s.al, ko=d4.rxns[[19]])
f19g <- FBA(sp.al, ko=d4.rxns[[19]], goal=goals, quiet=FALSE, se='G', control=list(trace=0))

###arg-4#################################################################################################################
arg4 <- "GLUTAMATE-N-ACETYLTRANSFERASE-RXN-L2R"
arg11 <- 'ACETYLORNDEACET-RXN-L2R'

f <- FBA(s.al, fba.ub=fva.ub)
v <- f$xopt
v[c(arg4, arg11)] #0.122, 0

f2 <- FBA(s.al, ko=arg4, fba.ub=fva.ub)
v2 <- f2$x
v2[c(arg4, arg11)] #0, 0.122

f3 <- FBA(s.al, ko=c(arg4, arg11)) #no growth

###vacuole###############################################################################################################
rr <- grep('[CCO-VAC-LUM]', model.rxns, fixed=TRUE, val=TRUE)
ff <- FBA(s.al, goal=rr)

###biolog################################################################################################################
model.rxns <- read.csv('model_rxns.csv')[,1]
s.al <- ss.mat(sg, c=intersect(model.rxns, colnames(sg)))
model.mets <- unique(gsub('\\[CCO-.+\\]$', '', rownames(s.al)))
    
pb <- predict.biolog(extra=TRUE, write.out=NULL)
pb <- pb[pb$cpd!='',]
pb <- pb[sapply(strsplit(gsub('\\[CCO-.+\\]$', '', pb$cpd), ','), FUN=function(x){ any(x %in% model.mets) }),]

pb.pm1 <- pb[pb$plate=='PM1',]

pm1-water-0.01-day3.tsv


###check model rxns that can't carry flux in vanilla FBA#################################################################
##one source of these is biosynthesis of cofactors which can't be used in FBA b/c they are recycled.
##it's ok to have some rxns in LIMED-FBA that are blocked in FBA, e.g. b/c of FBA recycles
model.rxns <- read.csv('model_rxns.csv')[,1]
s.al <- ss.mat(sg, c=model.rxns)

se <- rep('E', nrow(s.sp)); names(se) <- rownames(s.sp)
# 'MG-PROTOPORPHYRIN-MONOMETHYL-ESTER'
se[c('OXIDIZED-GLUTATHIONE', '2-HEXAPRENYL-6-METHOXY-14-BENZOQUINOL', 'OCTAPRENYL-METHYL-OH-METHOXY-BENZQ', '10-FORMYL-THF', 'HEME_D', 'CPD-675', 'SOLANESYL-PYROPHOSPHATE')] <- 'G' 

#f <- FBA(s.sp[rownames(s.al), colnames(s.al)], se=se[rownames(s.al)], goal="GLUTATHIONE-SYN-RXN-L2R")

#no flux rxns
max.nv <- max.n.flux(s.sp, se=se, allow.all.trans=TRUE)
max.nv2 <- max.n.flux(s.sp, se='G', allow.all.trans=TRUE)

(int <- intersect(model.rxns, names(max.nv)[max.nv<10**-6]))
intersect(model.rxns, names(max.nv2)[max.nv2<10**-6]) 

#(int2 <- intersect(int, nc.pwy$rxn))
int2 <- nc.pwy[nc.pwy$rxn %in% int, 3:1]
(int2 <- int2[order(int2$Pathway.Name),])

me <- min.export(s.sp, goals=int2$rxn[1])

i <- i+1
nfr <- int2$rxn[i]
ss.mat(s.sp, c=nfr)
me <- min.export(sp=ss.mat(s.sp, model.rxns), goal=nfr)

###check growth on vogel's using vanilla FBA#############################################################################
s2 <- s.sp[rownames(s.al), colnames(s.al)]
f <- FBA(s2, sense='E', fba.ub=fva.ub) #0.47; limed-fba gives 0.465

###parse meta chem.form##################################################################################################
cf <- smm.model$CHEM
cd=chemForm2df(cf)
w <- which(is.na(met.mat), arr.ind=TRUE)

w0=which(mm[1:8348,]!=met.mat[1:8348,], arr.ind=TRUE) #nunca!
smm[grep('CCO', rownames(smm.model)),][1:9,]

gmm=get.met.mat(global.smm=smm.meta, org.smm=smm.nc, s.rownames=rownames(s.sp), chem.form.col='CHEMICAL.FORMULA', colsums.mm.thresh=2)
all(gmm==mm.ss) #TRUE!

###cys-9 rescue##########################################################################################################
rad.rxns['NCU08352.5']
gg <- c('SO3[CCO-EXTRACELLULAR]', 'S2O3[CCO-EXTRACELLULAR]', 'MET[CCO-EXTRACELLULAR]', 'CYS[CCO-EXTRACELLULAR]')
se <- rep('E', nrow(s.al)); names(se) <- rownames(s.al)
se[gg] <- 'L'
se[se=='E'] <- 'G'

f <- FBA(s.al, ko=rad.rxns[['NCU08352.5']], goal=intersect(rownames(s.al), bm.goals$compound), sense=se)

f[c('GLUTCYSLIG-RXN', 'GLUTATHIONE-SYN-RXN', 'PRODISULFREDUCT-A-RXN', 'GLUTATHIONE-REDUCT-NADPH-RXN', 'RXN0-747', 'RIBONUCLEOSIDE-DIP-REDUCTII-RXN', 'RXN0-722', 'RXN2T-100')]

v0=c('GLUTCYSLIG-RXN', 'GLUTATHIONE-SYN-RXN', 'PRODISULFREDUCT-A-RXN', 'GLUTATHIONE-REDUCT-NADPH-RXN', 'RXN0-747', 'RIBONUCLEOSIDE-DIP-REDUCTII-RXN', 'RXN0-722', 'RXN2T-100')
v=paste(rep(v0, each=2), rep(c('-L2R', '-R2L'), times=length(v0)), sep='')
setdiff(v, colnames(s.al))
v=intersect(v, colnames(s.al))

f <- FBA(s.al, ko=rad.rxns[['NCU08352.5']], goal=v, sense=se)

###check perm vs ephemeral rxn stoich####################################################################################
s15e <- read.S('s15.perm.txt')
s15p <- read.S('s15.5.perm.txt')
map <- read.csv('ephemeral-permanent-map.csv', header=TRUE, as.is=TRUE)
mr <- read.csv('model_rxns.csv')[,1]

ss15e <- ss.mat(s15e, c=intersect(mr, colnames(ss15e)))
ss15p <- ss.mat(s15p, c=intersect(mr, colnames(ss15p)))

setdiff(map[,1], colnames(s15p))
length(intersect(map[,1], colnames(s15p)))
setdiff(rownames(ss15e), rownames(ss15p))

ss.mat(ss15e, r="ATP[CCO-MIT]")

all(ss15e==ss15p)

###check rxn bal#########################################################################################################
bal.mat <- read.csv('rxns_imbalance.csv', row.names=1)
bal.mat <- bal.mat[,c('C', 'O', 'N', 'P', 'S')]
unb.rxns <- rownames(bal.mat)[rowSums(bal.mat>0)>0]
(it=intersect(model.rxns, unb.rxns))

###make exchange rxns file###############################################################################################
ex.mets <- grep('\\[CCO-EXTRACELLULAR\\]$', rownames(s.al), value=TRUE)
exch.df <- data.frame(rxn=paste(gsub('\\[CCO-EXTRACELLULAR\\]$', '', x=ex.mets), 'TRANS-RXN-R2L', sep='-') , compound=ex.mets, coeff=-1)
write.table(exch.df, '/msc/neurospora/FBA/farm_data/Neurospora/nc10.export', sep='\t', quote=FALSE, row.names=FALSE)
#check
all(grep('-R2L$', vogel$rxn, value=TRUE) %in% exch.df$rxn) #TRUE

model.rxns <- read.csv('model_rxns.csv')[,1]
model.rxns <- union(model.rxns, exch.df$rxn)
write.csv(model.rxns, 'model_rxns.csv', row.names=FALSE)

###extracellular exports#################################################################################################
##n crassa can export from cyt->extracellular, but can't export from extracellular!
#added formate export to fix acetate
se2 <- rep('E', nrow(s.al)); names(se2) <- rownames(s.al)
c2=ckb[!is.na(ckb$pred) & ckb$pred<10**-6,]
w=rep(1, nrow(s.al)); names(w) <- rownames(s.al)
w['CPD-622'] <- 1000
w[grep('CCO-(MIT|GLYOXYSOME)', names(w), value=TRUE)] <- 1000
w[grep('CCO-EXTRACELLULAR', names(w), value=TRUE)] <- 10**-6

for (i in 1:4){ 
    cols <- setdiff(colnames(s.al), broad.rxns[[ c2$name[i] ]])
    me <- min.export(s.al[,cols], w=w, v.min=10**-6, ub=fva.ub[cols], ctrl=list(tilim=30, trace=0)) 
    cat('\n')
}

###ace-7#################################################################################################################
cn <- check.new.media(s.ck, rm='SUCROSE', add='XYLOSE[CCO-EXTRACELLULAR]', ctrl=list(trace=0), quiet=FALSE, ko=rad.rxns[['NCU09111.5']])

###dg0 consistency#######################################################################################################
#this only give 4 rxns which are part of pathways in correct direction
rr <- rxns.al[!is.na(rxns.al$dg) & !(rownames(rxns.al) %in% vogel$rxn) & (rxns.al$thermo<0.9 | rxns.al$dg0>50),]

###strip generics########################################################################################################
setdiff(grep('^[A-Z][a-z]', rownames(s.al), value=TRUE), grep('_bm$', rownames(s.al), value=TRUE))
ssm <- ss.mat(s.al, r=grep('^Fatty-Acids|^Nucleotides', rownames(s.al), v=T))
model.rxns <- setdiff(model.rxns, colnames(ssm))

###rm wrong lysine pwys##################################################################################################
#pwy-5327
pwys <- c('P163-PWY','LYSDEGII-PWY','LYSINE-DEG1-PWY','PWY-5280','PWY-5283','PWY-5298','PWY-5311','PWY-5314','PWY-5324','PWY0-461')
ncp <- nc.pwy[nc.pwy$PathwayID %in% pwys,]
x <- intersect(ncp$rxn, colnames(s.al))
y <- c('1.5.1.8-RXN','1.5.1.9-RXN','2-AMINOADIPATE-AMINOTRANSFERASE-RXN','2-KETO-ADIPATE-DEHYDROG-RXN','RXN-8173','RXN-8162','L-LYSINE-AMINOTRANSFERASE-RXN','ACETYL-COA-ACETYLTRANSFER-RXN','ACETATEKIN-RXN','LYSINE-23-AMINOMUTASE-RXN')
(x <- x[!(gsub('-L2R|-R2L', '', x) %in% y)])

###imbalance#############################################################################################################
bal.mat <- read.csv('rxns_imbalance.csv', row.names=1)
bmat <- bal.mat[,c('C', 'H', 'O', 'N', 'P', 'S')]
bmat <- bmat[colnames(s.al),]
bmat <- bmat[rowSums(abs(bmat))>0,]
bmat$path <- NA
for (i in 1:nrow(bmat)){
    bmat$path[i] <- nc.pwy$Pathway.Name[match(rownames(bmat)[i], nc.pwy$rxn)]
}
#filter rxns that use thioredoxin, ferrodoxin, cyt c, trna, etc.
smm <- read.csv('smm.csv'); rownames(smm) <- smm$FRAME; smm <- smm[rownames(s.al),]
mm.ss <- as.matrix(read.csv('met_mat.csv', row.names=1)); mm <- mm.ss[rownames(s.al),]
(cpds.noelems <- rownames(mm)[rowSums(abs(mm[,c('C', 'H', 'O', 'N', 'P', 'S','CL', 'K', 'MG', 'FE', 'NA.')]))==0])
bmat <- bmat[setdiff(rownames(bmat), colnames(ss.mat(s.al, r=cpds.noelems))),]

###export dilution#######################################################################################################
#add exports st model doesn't rely on 'export dilution' of OH and HCO3.
add = c('RXN0-5224-R2L', 'RXN0-5224-R2L[CCO-MIT]', "RXN2T-70-R2L[CCO-MIT]")
model.rxns <- union(model.rxns, add)
#write.csv(model.rxns, 'model_rxns.csv', row.names=FALSE)

###limed-fba ub's########################################################################################################
f=FBA(s.ck, ko=ko.lst[[9]], fba.ub=rep(1000, ncol(s.ck)))
rhs=(s.ck %*% f$x)[,1]; names(rhs)=rownames(s.ck); rhs=rhs[rhs>0]; sort(rhs)

###acet, round 3#########################################################################################################
spa <- changeNuts(s.al, rm='SUCROSE', add='ACET[CCO-EXTRACELLULAR]')
f <- FBA(spa, se='E') #growth, but no growth w/ se='E'
rhs=(spa %*% f$x)[,1]; names(rhs)=rownames(spa); rhs=rhs[rhs>0]; sort(rhs)
me <- min.export(spa) #needs to export HCO3

n <- numeric(nrow(spa)); names(n) <- rownames(spa)
ssm <- ss.mat(s.sp, 'RXN0-5224-R2L')
n[paste(names(ssm), '[CCO-MIT]', sep='')] <- ssm
spa <- cBind(spa, n)

###no.v##################################################################################################################
ap <- apply(as.matrix(1:21), 1, FUN=function(i){
    check.new.media(s.ck[,!(colnames(s.ck) %in% no.v[i])], rm='SUCROSE', add='ACET[CCO-EXTRACELLULAR]', ctrl=list(trace=0), quiet=FALSE)$obj
})

##min export
ssm <- read.csv('smm.csv', row.names=1)
smm$MASS[is.na(smm$MASS)] <- 10**5
ssp <- s.sp[,ss]; ssp <- ssp[rowSums(abs(ssp))>0,]
ssp2 <- ssp[-grep('CCO-EXTRACELLULAR', rownames(ssp)),]
#w <- rep(1, nrow(ssp2)); names(w) <- rownames(ssp2); w[intersect(names(w), rownames(smm))] <- smm[intersect(names(w), rownames(smm)), 'MASS']
#me <- min.export(ssp2, w=w, goals="FADSYN-RXN-L2R")

model.rxns <- setdiff(model.rxns, 'ACETYL-COA-HYDROLASE-RXN-L2R[CCO-GLYOXYSOME]')
model.rxns <- setdiff(model.rxns, no.v)

###no limed##############################################################################################################
ssp <- s.sp[,ss]; ssp <- ssp[rowSums(abs(ssp))>0,]
ssp[c('NADH', 'FMN'), 'biomass']=-0.001 #'ACETYL-GLU'
#me <- min.export(ssp)
#broad
ckb <- check.ko(s.test=ssp, ko.lst=broad.rxns, sense='E')
(ckb2 <- ckb[!is.na(ckb$pred) & ckb$pred==0,])
#rad
ck <- check.ko(s.test=ssp, ko.lst=rad.rxns, sense='E', print.v=FALSE, annot.df=rad0[names(rad.rxns),], ub=rep(10**3, ncol(ssp)))
(ck2 <- ck[ck$pred >= 10**-6,]); (rad0[ck2$name,])

###model.rxns in metacyc#################################################################################################
##src make_s.r & farm_script.r

#get rxn strings
mr=setdiff(intersect(model.rxns, colnames(s.sp)), c(nc$rxn, bm$rxn, vogel$rxn))
rxns.v <- apply(s.sp[,union(model.rxns, nc$rxn)], 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x, sep=':', collapse=';') })

#get metacyc model rxns
#(mr.rxns <- gsub('NADP', 'NAD', rxns.v[mr]))
(mr.rxns <- rxns.v[mr][rxns.v[mr]!=''])
(int=intersect(mr.rxns, rxns.v[nc$rxn]))
(mr.rxns <- mr.rxns[mr.rxns %in% int])
(nc.rxns <- rxns.v[unique(nc$rxn)][ rxns.v[unique(nc$rxn)] %in% int ])
#check
all(sort(mr.rxns)==sort(nc.rxns))

rck.rxns=union(setdiff(rck.rxns, names(mr.rxns)), names(nc.rxns))

(ss1=ss.mat(sg, names(sort(mr.rxns))))
(ss2=ss.mat(sg, names(sort(nc.rxns))))
all(ss1==ss2)

model.rxns=union(setdiff(model.rxns, c(RXN-12084-L2R, RXN-12084-R2L, RXN0-276-L2R, RXN-10862-L2R, RXN-11726-R2L, THYMIDINE-TRIPHOSPHATASE-RXN-L2R)), 
c(ESTRADIOL-17-BETA-DEHYDROGENASE-RXN-NADP/CPD-352//NADPH/ESTRONE/PROTON.35.-L2R, ESTRADIOL-17-BETA-DEHYDROGENASE-RXN-NADP/CPD-352//NADPH/ESTRONE/PROTON.35.-R2L, 
RXN-2962-S-HYDROXYMETHYLGLUTATHIONE/NAD//CPD-548/NADH/PROTON.52.-L2R, TRANSENOYLCOARED-RXN-BUTYRYL-COA/NADP//CROTONYL-COA/NADPH/PROTON.44.-R2L, 
NUCLEOSIDE-DIPHOSPHATASE-RXN-ADP/WATER//Pi/AMP/PROTON.25.-L2R, NUCLEOSIDE-TRIPHOSPHATASE-RXN-TTP/WATER//Pi/TDP/PROTON.25.-L2R))

print(paste("rid rxns = c(", paste(mr[10:14], collapse="\', \'"), ")", sep="\'"))
paste(names(nc.rxns), collapse=', ')

#RXN-11726-R2L is needed for growth on acet, even tho 'TRANSENOYLCOARED-RXN-BUTYRYL-COA/NADP//CROTONYL-COA/NADPH/PROTON.44.-R2L' should replace it

#rck.rxns=union(setdiff(rck.rxns, "NADPH-PEROXIDASE-RXN-L2R"), "NADH-PEROXIDASE-RXN-L2R")
#rck.rxns=union(setdiff(rck.rxns, "RXN-8672-L2R"), c('D-AMINO-ACID-OXIDASE-RXN-D-ALANINE/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/PYRUVATE/PROTON.75.-L2R[CCO-MIT]',
#'TRANS-RXN2T-49-CPD-218//CPD-218.17.-L2R[CCO-MIT]'))
#model.rxns=union(setdiff(model.rxns, c(RXN-8672-L2R, NADPH-PEROXIDASE-RXN-L2R)), c(D-AMINO-ACID-OXIDASE-RXN-D-ALANINE/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/PYRUVATE/PROTON.75.-L2R[CCO-MIT],
#TRANS-RXN2T-49-CPD-218//CPD-218.17.-L2R[CCO-MIT], NADH-PEROXIDASE-RXN-L2R))

model.rxns=union(model.rxns, c('RXN-11820-TETRACOSANOATE/HYDROGEN-PEROXIDE//CPD-8472/WATER.49.-L2R', 'RXN-11820-TETRACOSANOATE/HYDROGEN-PEROXIDE//CPD-8472/WATER.49.-R2L'))

#export errors:
#'RXN-4464' instances have same cpd on both sides
#carbpsyn-rxn doesn't show up in cytosol

###fix new trans rxns####################################################################################################
nr=c('RXN-9615-L2R[CCO-MIT]', 'TRANS-RXN2T-32-R2L[CCO-GLYOXYSOME]', 'TRANS-RXN2T-32-L2R[CCO-GLYOXYSOME]', 'TRANS-RXN2T-32-R2L[CCO-MIT]', 'TRANS-RXN2T-32-L2R[CCO-MIT]',  'TRANS-RXN2T-33-L2R[CCO-MIT]', 'TRANS-RXN2T-33-R2L[CCO-GLYOXYSOME]', 'TRANS-RXN2T-33-L2R[CCO-GLYOXYSOME]', 'RXN2T-68-L2R[CCO-MIT]', 'TRANS-RXN2T-82-R2L[CCO-MIT]', 'TRANS-RXN2T-83-R2L[CCO-MIT]', 'TRANS-RXN2T-233-L2R[CCO-EXTRACELLULAR]', 'TRANS-RXN2T-236-L2R[CCO-MIT]', 'TRANS-RXN2T-241-L2R[CCO-GLYOXYSOME]')

for (i in 1:length(nr)){
    cat(i, grep(nr[i], colnames(s.sp), v=T, f=T), '\n \n')
}

(ssm=ss.mat(s.sp, nr))

###arg mutants###########################################################################################################
rr=rad.rxns[['NCU07732.5']]
se=rep('E', nrow(s.al)); names(se)=rownames(s.al); se['L-ORNITHINE[CCO-MIT]']='L'
f=FBA(s.sp[rownames(s.al), colnames(s.al)], ko=rr,  goal=unique(bm$comp[bm$coeff<0]), se=se)
f=FBA(s.al, ko=rr,  goal='CARBAMOYL-P'), se=se)

#done: hco3, citrulline out, 
#add "TRANS-RXN2T-234-GLN/PROTON//GLN/PROTON.23.-L2R[CCO-MIT]", "TRANS-RXN2T-49-GLT//GLT.9.-R2L[CCO-PM-FUNGI]", "TRANS-RXN2T-193-L2R[CCO-MIT]"
#move ORNCARBAMTRANSFER-RXN to cco-mit in model.rxns
grep('ORNCARBAMTRANSFER-RXN', colnames(s.sp), v=T)

model.rxns <- read.csv('model_rxns.csv')[,1]
model.rxns=setdiff(union(model.rxns, c('CARBPSYN-RXN-L2R[CCO-MIT]',"ORNCARBAMTRANSFER-RXN-L2R[CCO-MIT]", "TRANS-RXN2T-234-GLN/PROTON//GLN/PROTON.23.-L2R[CCO-MIT]", "TRANS-RXN2T-49-GLT//GLT.9.-R2L[CCO-MIT]", "TRANS-RXN2T-193-L2R[CCO-MIT]", "TRANS-RXN2T-193-R2L[CCO-MIT]")), c("ORNCARBAMTRANSFER-RXN-L2R", "ORNCARBAMTRANSFER-RXN-R2L"))

f$x[grep('ORNCARBAMTRANSFER-RXN', colnames(s.al), v=T)]

check.rxn.inputs(s.al, 'CARBPSYN-RXN-L2R[CCO-MIT]')
check.rxn.inputs(s.al, 'ORNCARBAMTRANSFER-RXN-L2R[CCO-MIT]')
f2=FBA(s.al, se=se)

me=min.export(s.al, goals='CARBAMOYL-P[CCO-MIT]')

ff=FBA(s.al, ko=rad.rxns[['NCU07732.5']], goal='L-CITRULLINE')
tm=trace.mets(
(fom=flux.of.met(s.al, ff$x[1:ncol(s.al)], 'L-CITRULLINE'))

###growth on ace, round 2################################################################################################
seg=rep('G', nrow(sg)); names(seg)=rownames(sg)
seg[c('ACET[CCO-EXTRACELLULAR]')]='L'
#seg['HCO3[CCO-MIT]']='G'
#seg['ATP']='L'
tp=rxns$thermo.prob; tp[is.na(tp)]=1
pen=1.1-rxns$prob * tp; names(pen)=rownames(rxns)
pen[colnames(s.al)]=0; pen[is.na(pen)]=1.1
pen[c("MALSYN-RXN-L2R", "CITSYN-RXN-L2R")]=10
f.ub=rep(5*10**4, ncol(sg)); names(f.ub)=colnames(sg)
f.ub['SUCROSE-TRANS-RXN-L2R']=0
f.ub[c('PYRUVFORMLY-RXN-R2L','FRUCTOSE-6-PHOSPHATE-PHOSPHOKETOLASE-RXN-R2L')]=0

f=FBA(s.al, fba.ub=f.ub[colnames(s.al)], se=seg[rownames(s.al)], ko=rad.rxns[['NCU09873.5']])
(fom=flux.of.met(s.al, f$x, '2-PG'))
(tm=trace.mets(s.al, f$x, from='ACET', to='2-PG', pwy=nc.pwy,
rm.mets.grep='^((A|C|G|T|U)(M|D|T)P|^NAD(P|)(H|)|FAD(H.*|)|FMN|CO-A|.+redoxin(s|)|Cytochromes.*|AMMONIA|NITR(A|I)TE|GLT|GLN|Pi|PPI|PROTON|WATER|CARBON-DIOXIDE|HYDROGEN-PEROXIDE|OXYGEN-MOLECULE)(\\[|$)'))

me=min.export(s.al[,!(colnames(s.al) %in% 'SUCROSE-TRANS-RXN-L2R')], se=seg[rownames(s.al)])

#model.rxns=setdiff(union(model.rxns, "TRANS-RXN2T-220-L2R[CCO-GLYOXYSOME]"), c('TRANS-RXN2T-233-L2R[CCO-GLYOXYSOME]', 'TRANS-RXN2T-233-R2L[CCO-GLYOXYSOME]'))
model.rxns=setdiff(union(model.rxns, 'TRIOSEPISOMERIZATION-RXN-L2R'), c('PYRUVFORMLY-RXN-R2L','FRUCTOSE-6-PHOSPHATE-PHOSPHOKETOLASE-RXN-R2L'))

model.rxns=setdiff(union(model.rxns, TRIOSEPISOMERIZATION-RXN-L2R), c(PYRUVFORMLY-RXN-R2L,FRUCTOSE-6-PHOSPHATE-PHOSPHOKETOLASE-RXN-R2L))

(mp=min.pen(sp=sg, pen=pen, f.=f.ub, se=seg))
(mpm <- min.pen.milp(sp=sg, pen.v=pen, rxns0=rownames(mp), n.add=1, f.=f.ub, se=seg, ctrl=list(tilim=180)))

##ace supp###############################################################################################################
(tt=grep('TRANS-RXN2T-191', colnames(s.al), v=T))
(ts=grep('TRANS-RXN2T-156', colnames(s.al), v=T))

# se=rep('E', nrow(s.al)); names(se)=rownames(s.al)
# se[c('ACET')]='L'
# se['FADH2[CCO-GLYOXYSOME]']='G'
f=FBA(s2, ko=c(rad.rxns[['NCU06482.5']]), fba.ub=c(fva.ub, rep(Inf, 1), 0))
f2=FBA(s2, ko=c(rad.rxns[['NCU06482.5']]), fba.ub=c(fva.ub, rep(Inf, 1), 100))
f3=FBA(s2, ko=c(rad.rxns[['NCU06482.5']]), goal=unique(bm$comp[bm$coeff<0]))
ff=FBA(s2, ko=c(rad.rxns[['NCU06482.5']]), goal='CIT')

seg=rep('E', nrow(sg)); names(seg)=rownames(sg)
seg[c('ACET[CCO-EXTRACELLULAR]')]='L'
pen=1.1-rxns$prob; names(pen)=rownames(rxns)
pen[colnames(s.al)]=0; pen[c(tt, 'CITSYN-RXN-L2R', 'ISOCIT-CLEAV-RXN-R2L', "PYRUVDEH-RXN-L2R[CCO-MIT]", "OXOGLUTARATE-DEHYDROGENASE-NADP+-RXN-R2L", "GLUTRNAREDUCT-RXN-L2R")]=10**6
(mp=min.pen(sg, pen=pen, se=seg))

#r=rownames(mp)
ss.mat(s.sp, r)
rxns[r,]

model.rxns <- read.csv('model_rxns.csv')[,1]
model.rxns=setdiff(model.rxns, c('RXN-9951-R2L', 'RXN-9951-L2R', 'ISOCITDEH-RXN-L2R', 'RXN-8642-L2R', 'RXN-8642-R2L'))

cit.sp=Matrix(0, nrow=nrow(s.al), ncol=2); rownames(cit.sp)=rownames(s.al); 
cit.sp[c('CIT','CIT[CCO-GLYOXYSOME]'),1]=c(1,-1)
# cit.sp[c('CIT','CIT[CCO-MIT]'),2]=c(-1,1)
# cit.sp[c('2-KETOGLUTARATE','2-KETOGLUTARATE[CCO-MIT]'),3]=c(1,-1)
cit.sp['ACET',2]=1
s2=cBind(s.al, cit.sp); colnames(s2)[(ncol(s2)-ncol(cit.sp)+1):ncol(s2)]=c('TRANS-CIT-RXN-L2R', 'ACET-IN')

# f=FBA(s.al, ko=c('RXN66-3-L2R','RXN0-3962-L2R', rad.rxns[['NCU06482.5']]), goals='ACET', min2=TRUE)
fom=flux.of.met(s.al, f$x[colnames(s.al)], 'CIT'); sort(fom[fom>0])
(tm=trace.mets(s.al, v=f$x[colnames(s.al)], f='SUCROSE[CCO-EXTRACELLULAR]', t='ACET', pwy=nc.pwy))
rxns[tm$rxn[tm$rxn!=''],]

#supp

####sensitivity to bm coeffs#############################################################################################
mma=min.met.adj(s.al, 

f <- FBA(s.al, sense='G', fba.ub=fva.ub)
lam=f$extra$lam; names(lam)=rownames(s.al)
lam[order(lam)][1:19]

###check filtered model rxns############################################################################################# 
rr=rxns[intersect(colnames(s.al), filt.rxns$rxn), c('Pathways', 'Reason')]
rr=rr[rr$Reason %in% c('RXNS-WITH-INSTANCELESS-CLASSES', 'POLYMERIZATION-PWY-RXNS', 'INSTANTIABLE-GENERIC-RXNS', 'CANNOT-BALANCE-RXNS'),]
(rr=rr[order(rr$Reason),])

model.rxns=setdiff(model.rxns, 'R601-RXN-L2R[CCO-MIT]')
f$x['R601-RXN-L2R[CCO-MIT]']

###missed rxns###########################################################################################################
rr=unique(rxns$rxn[rxns$nc==1 & rxns$gpr==1 & rxns$prob>0.9 & !is.na(rxns$dg0) & rxns$dg0<5])
r2=setdiff(rr, c(rxns.al$rxn, filt.rxns$rxn))

#redundant rxns
rxns.v <- apply(s.sp[,union(model.rxns, r2)], 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x, sep=':', collapse=';') })
r3 <- r2[!(rxns.v[r2] %in% rxns.v[model.rxns])]

#no.v
ss2 <- intersect(union(model.rxns, r3), colnames(sg))
s2 <- sg[,ss2]; s2 <- s2[rowSums(abs(s2))>0,]
max.nv <- max.n.flux(s2, se='E', allow.all.trans=TRUE)
no.v <- names(max.nv)[max.nv==0]
(r4 <- setdiff(r3, no.v))

#NOTE some rxns in this set b/c: use NADP instead of NAD; are associated w/ EC number, but not gene, based on gene annotation;

rck.rxns=union(rck.rxns, c("L-AMINO-ACID-OXIDASE-RXN-GLY/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/GLYOX/PROTON.66.-L2R", "TRANS-RXN0-277-L2R[CCO-EXTRACELLULAR]"))

added L-AMINO-ACID-OXIDASE-RXN-GLY/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/GLYOX/PROTON.66.-L2R & TRANS-RXN0-277-L2R[CCO-EXTRACELLULAR]
###gapfind###############################################################################################################
gf=gapfind(s.al)
np=gf$noprod.root
np=np[order(nchar(np))]
(np.mit=grep('\\[CCO-MIT\\]$', np, val=T))
(trans.np.mit=grep('TRANS-', colnames(ss.mat(s.al, r=np.mit)), v=T))
gf=gapfind(s.al[,!(colnames(s.al) %in% trans.np.mit)])

###add export rxns#######################################################################################################
s2=s.al[,model.rxns]
r.ox.thio=r.1063=numeric(nrow(s2))
names(r.ox.thio)=names(r.1063)=rownames(s2)
r.ox.thio['Ox-Thioredoxin']=r.1063['CPD-1063']=-1
s3=cBind(s2, r.ox.thio, r.1063)
ff=FBA(s3)
ck <- check.ko(s.test=s3, ko.lst=broad.rxns, sense='E'); table(ck$obs, ck$pred>10^-6)

ko.rxns.lst=list(pyrdeh="PYRUVDEH-RXN-L2R[CCO-MIT-LUM]")

(aao=grep('.-AMINO-ACID-OXIDASE-RXN.+-R2L$', colnames(s.sp), v=T))
write.csv(farm.rid, 'farm_rid_rxns.csv', row.names=FALSE)

s2=s.ck[,!(colnames(s.al) %in% c('CO-A-GOAL', ntp, fac))]
ff=FBA(s2)

(ntp=grep('NUCLEOSIDE-TRIPHOSPHATASE-RXN.+-R2L', colnames(s.ck), v=T))
(fac=grep('RXN-9917.+-R2L', colnames(s.ck), v=T))
fac=RXN-9917-LINOLENIC_ACID/ATP/CO-A//CPD2T-28/PPI/AMP.42.-R2L

mm[intersect(names(mm), bm$comp),]

f=FBA(s.al, ko=c('ACETOACETATE--COA-LIGASE-RXN-R2L', 'RXN2T-62-CPD1G-567/WATER//ARACHIDIC_ACID/CO-A/PROTON.44.-R2L'))
ss=setdiff(ss, c('ACETOACETATE--COA-LIGASE-RXN-R2L', 'RXN2T-62-CPD1G-567/WATER//ARACHIDIC_ACID/CO-A/PROTON.44.-R2L', 'PHOSACETYLTRANS-RXN-L2R', 'TRANS-RXN0-277-L2R[CCO-EXTRACELLULAR]'))
#could use 'ACETATEKIN-RXN-R2L' instead of 'PHOSACETYLTRANS-RXN-L2R'; together they make ATP+acetate+coa from ADP+phosphat+acetyl-coa, with a combined dG<0
ss=setdiff(ss, c('ACETOACETATE--COA-LIGASE-RXN-R2L', 'RXN2T-62-CPD1G-567/WATER//ARACHIDIC_ACID/CO-A/PROTON.44.-R2L', 'ACETATEKIN-RXN-R2L', 'TRANS-RXN0-277-L2R[CCO-EXTRACELLULAR]'))

ss=setdiff(ss, c('ACETOACETATE--COA-LIGASE-RXN-R2L', 'RXN2T-62-CPD1G-567/WATER//ARACHIDIC_ACID/CO-A/PROTON.44.-R2L', 'TRANS-RXN0-277-L2R[CCO-EXTRACELLULAR]', 'ACETYL-COA-HYDROLASE-RXN-R2L'))
#if set to al.lb[c("beta_ACETATEKIN-RXN-R2L", "beta_PHOSACETYLTRANS-RXN-L2R")]=1, get ACETYL-COA-HYDROLASE-RXN-R2L, which solves am & ace-8 synergistically

dg0=rxns.al$dg0; dg0[is.na(dg0)]=10
pen.v[rownames(rxns.al)[dg0>9]]=dg0[dg0>9]
pen.v[c('ATPSYN-RXN-L2R[CCO-EXTRACELLULAR]')]

thresh <- 10**-9
rr <- rxns.al[beta>thresh,]
rxns1 <- union(rownames(rr), names(pen.v)[pen.v==0])
print(paste('Testing if new rxn set gives FBA growth w/ threshold', thresh))
fba <- FBA(s.al[,rxns1], control=list(trace=0), sense='E')
cat('Number of bad Gibbs rxns to be added:', length(gibbs.add <- intersect(rownames(rxns.al)[rxns.al$dg0>20], rxns1)), '\n')
gibbs.rm=setdiff(rownames(rxns.al)[dg0>20], gibbs.add)

nr=setdiff(names(pen.v)[rxns$nc.pwy==1 & pen.v==0], model.rxns)
rck.rxns=union(rck.rxns, nr)
#pwy rxns which sink mass are ok
nr2=names(pen.v)[rxns$nc.pwy==1 & pen.v==0.1 & names(pen.v) %in% rownames(filt.rxns)[filt.rxns$Reason=='CANNOT-BALANCE-RXNS'] & !(names(pen.v) %in% unb.rxns)]
s.ck['CPD-12930', 'biomass']=-0.1

d=f2$x[f2$x>0 & f$x==0]

 rr=c("TRANS-RXN2T-94-Pi//Pi.7.-L2R[CCO-PM-FUNGI]",
"TRANS-RXN2T-92-SUCROSE//MALTOSE.17.-L2R[CCO-PM-FUNGI]",
"TRANS-RXN2T-51-LINEAR-ALPHA-GLUCAN//LINEAR-ALPHA-GLUCAN.41.-L2R[CCO-PM-FUNGI]",
"TRANS-RXN2T-74-L2R[CCO-PM-FUNGI]",
"TRANS-RXN2T-232-R2L[CCO-MIT]",
"TRANS-RXN2T-232-R2L[CCO-EXTRACELLULAR]")

tr=grep('^TRANS-RXN.+//', colnames(s.al), v=T); pen.v[tr]=1
tr92=grep('^TRANS-RXN2T-92-', tr, v=T); pen.v[tr92]=1000

grep('^TRANS-RXN2T-232-', colnames(s.al), value=TRUE),
al.ub["beta_RXN-10862-R2L"]=0

g=grep('\\[CCO-MIT\\]$|\\[CCO-GLYOXYSOME\\]$', colnames(s.sp), value=TRUE)
model.rxns=union(model.rxns, c(vogel$rxn, g))

int=intersect(g, no.v)
setdiff(grep('CCO-GLYOXYSOME', int , value=TRUE), grep('^TRANS-RXN|POLYMER-INST-', int , value=TRUE))

grep("DIHYDROOROTATE-DEHYDROGENASE-RXN-", colnames(s.sp), value=TRUE)

ssm=ss.mat(s.al, "MALSYN-RXN-L2R[CCO-GLYOXYSOME]")
f=FBA(s.al, goals=names(ssm)[ssm<0])
f$x[paste(names(ssm)[ssm<0], 'GOAL-RXN-L2R', sep='-')]

me=min.export(s.al, w=rxns.al$mass)

se=rep('E', nrow(s.al)); names(se)=rownames(s.al)
se[nm0[c(1:3)]]='G'
f2=FBA(s.al, sense=se)

(aao.l2r=grep('.-AMINO-ACID-OXIDASE-RXN-.+-L2R', colnames(s.sp), value=TRUE))
"D-AMINO-ACID-OXIDASE-RXN-D-ALANINE/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/PYRUVATE/PROTON.75.-L2R"

ssm=ss.mat(s.sp, r="CPD-218")
rxns[names(ssm), 'nc']

c('TRANS-RXN2T-49-CPD-218//CPD-218.17.-L2R[CCO-PM-FUNGI]', 'D-AMINO-ACID-OXIDASE-RXN-CPD-218/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/CPD-479/PROTON.72.-L2R')

###check growth on other media###########################################################################################
se <- rep('E', nrow(s.al)); names(se) <- rownames(s.al)
se['ACET[CCO-EXTRACELLULAR]']='L'
f=FBA(s.al, ko='SUCROSE-TRANS-RXN-L2R', se=se)
cn=check.new.media(s.al, rm='SUCROSE', add=c('ACET[CCO-EXTRACELLULAR]'), goal="ISOCIT-CLEAV-RXN-L2R[CCO-GLYOXYSOME]", se='G')

cn=check.new.media(s.al, rm='SUCROSE', add=c('ACET[CCO-EXTRACELLULAR]', 'ATP', 'ERGOSTEROL', 'FRUCTOSE-16-DIPHOSPHATE', 'TETRACOSANOATE'), se='G')
cn=check.new.media(s.al, rm='SUCROSE', add=c('ACET[CCO-EXTRACELLULAR]', 'ATP', 'DELTA3-ISOPENTENYL-PP', 'PALMITATE[CCO-GLYOXYSOME]', 'GLYCERALD'), se='G')

s.all=sg
colnames(s.all)[colnames(s.all)=='SUCROSE-TRANS-RXN-L2R']=colnames(s.al)[colnames(s.al)=='SUCROSE-TRANS-RXN-L2R']='ACET-TRANS-RXN-L2R'
s.all[c('ACET[CCO-EXTRACELLULAR]', 'SUCROSE[CCO-EXTRACELLULAR]'), 'ACET-TRANS-RXN-L2R']=s.al[c('ACET[CCO-EXTRACELLULAR]', 'SUCROSE[CCO-EXTRACELLULAR]'), 'ACET-TRANS-RXN-L2R']=1:0
#f=FBA(s.al)
rxns.all <- read.delim('rxn_annot.txt')
rownames(rxns.all)[rownames(rxns.all)=='SUCROSE-TRANS-RXN-L2R']='ACET-TRANS-RXN-L2R'
vogel$rxn[vogel$rxn=='SUCROSE-TRANS-RXN-L2R']='ACET-TRANS-RXN-L2R'
model.rxns=union(model.rxns, 'ACET-TRANS-RXN-L2R')
# model.rxns=union(model.rxns, c('2TRANSKETO-RXN-R2L', 'ADENYL-KIN-RXN-R2L', 'FRUCTOSE-6-PHOSPHATE-PHOSPHOKETOLASE-RXN-R2L', 'RXN0-5408-R2L', 'RXN3O-127-R2L', 'SARCOX-RXN-R2L'))
rownames(rxns.al)[rownames(rxns.al)=='SUCROSE-TRANS-RXN-L2R']='ACET-TRANS-RXN-L2R'

ko.rxns.lst=rad.rxns['NCU09873.5'] #acu-6
rm(s.all)

k=rownames(rxns.al)[((rxns.al$gpr==1 & rxns.al$nc.pwy==1) | rxns.al$prob>0.9) & (is.na(rxns.al$dg0) | rxns.al$dg0<6)]
keep.rxns=union(keep.rxns, k)
al.lb[intersect(paste('beta', unique(keep.rxns), sep='_'), names(al.lb))] <- 1

fva.ub['ACET-TRANS-RXN-L2R']
ko=rid2
f=FBA(s.al, fba.ub=fva.ub)
f2=FBA(s.al, ko=c(ko.rxns.lst[[1]]), fba.ub=fva.ub)
f3=FBA(s.al, ko=c(ko.rxns.lst[[1]], 'MALIC-NADP-RXN-L2R', 'MALIC-NAD-RXN-L2R'), fba.ub=fva.ub)
(tm <- trace.mets(s.al, from='ACET', to='PYRUVATE', v=f3$x, pwy=nc.pwy))

model.rxns=union("ACETATE--COA-LIGASE-RXN-L2R", setdiff(model.rxns, c('MALIC-NADP-RXN-L2R', 'MALIC-NAD-RXN-L2R', "TRANS-RXN2T-220-L2R[CCO-GLYOXYSOME]")))

model.rxns <- read.csv('model_rxns.csv')[,1]
model.rxns=union("ACETATE--COA-LIGASE-RXN-L2R", setdiff(model.rxns, c("TRANS-RXN2T-220-L2R[CCO-GLYOXYSOME]", 'TRANS-RXN2T-220-R2L[CCO-GLYOXYSOME]',
"ACETATE--COA-LIGASE-RXN-L2R[CCO-GLYOXYSOME]", 'TRANS-RXN2T-5-ADP//ADP.9.-R2L[CCO-MIT]', 'TRANS-RXN2T-5-ADP//ADP.9.-L2R[CCO-MIT]', 'TRANS-RXN2T-5-ATP//ATP.9.-R2L[CCO-MIT]')))

'new model.rxns. solved acu-6 by moving ACETATE--COA-LIGASE-RXN from glyox to cytosol & removing glyox coa transport & mit adp trans.' 

s.al[c('PROTON', 'PROTON[CCO-MIT]'), 'TRANS-RXN2T-5-ATP//ATP.9.-L2R[CCO-MIT]']=c(-1,1)

ssm=ss.mat(s.sp[,colnames(s.al)], r='CO-A[CCO-GLYOXYSOME]')
ssm[ssm>0]
ss.mat(s.sp[,colnames(s.al)], names(ssm)[ssm>0])

#adenine phosphates
am=colSums(s.al[paste('A', c('M','D','T'), 'P[CCO-MIT]', sep=''),]); names(am)=colnames(s.al)
(am=am[abs(am)>0.9])
ss.mat(s.al, names(am))
trace.mets(s.al, from='AMP[CCO-MIT]', to='ATP[CCO-MIT]')

al.ub[paste('beta_', c('1.4.3.19-RXN-L2R', 'METHIONINE--GLYOXYLATE-TRANSAMINASE-RXN-R2L'), sep='')]=0

cn=check.new.media(s.al, rm.nuts='SUCROSE', add.nuts='ACET[CCO-EXTRACELLULAR]', goals='OXYGEN-MOLECULE', se='G')
check.rxn.inputs(sp, "ACETATE--COA-LIGASE-RXN-L2R[CCO-GLYOXYSOME]")
f=FBA(sp, goals=c('ATP', 'CO-A'))
(mma=min.met.adj(sp))
(mb=mma.bin(sp, M=5*10**3, wts=mm[rownames(sp),2], ctrl=list(tilim=30)))

cn=check.new.media(s.al, rm='SUCROSE', add='GLUCOSE', 'ATP'), se='G')

se <- rep('E', nrow(s.al)); names(se) <- rownames(s.al)
se['CPD-218[CCO-EXTRACELLULAR]']='L'
f=FBA(s.al, ko='SULFATE-TRANS-RXN-L2R', se=se)
cn=check.new.media(s.al, rm='SULFATE', add='CPD-218[CCO-EXTRACELLULAR]')

se <- rep('E', nrow(s.al)); names(se) <- rownames(s.al)
se['BUTYRIC_ACID[CCO-EXTRACELLULAR]']='L'
f=FBA(s.al, ko='SUCROSE-TRANS-RXN-L2R', se=se)
flux.of.met(s.al, f$x, 'BUTYRIC_ACID[CCO-GLYOXYSOME]')

f=FBA(s.al, goal='ACETOACETYL-COA[CCO-GLYOXYSOME]')
v=f$x; v=v[!(names(v) %in% 'ACETOACETYL-COA[CCO-GLYOXYSOME]-GOAL-RXN-L2R')]
flux.of.met(s.al, v, 'ACETOACETYL-COA[CCO-GLYOXYSOME]')

sp2=s.al
sp2[c('SUCROSE[CCO-EXTRACELLULAR]', 'BUTYRIC_ACID[CCO-EXTRACELLULAR]'), 'SUCROSE-TRANS-RXN-L2R']=0:1
mma=min.met.adj(sp2)
mb=mma.bin(sp2, mets0=names(mma))

max.nv <- max.n.flux(sp2, se='E', allow.all.trans=FALSE)
no.v <- names(max.nv)[max.nv==0]
(ssm=ss.mat(sp2, r='BUTYRIC_ACID'))
names(ssm)[names(ssm) %in% no.v]

me=min.export(s.al, rxns.min="ACYLCOASYN-RXN-BUTYRIC_ACID/ATP/CO-A//BUTYRYL-COA/PPI/AMP.43.-L2R[CCO-GLYOXYSOME]")

cn=check.new.media(s.al, add='BUTYRIC_ACID[CCO-EXTRACELLULAR]', goal='ACETOACETYL-COA[CCO-GLYOXYSOME]')

model.rxns=union(model.rxns, c('RXN2T-72-L2R', "CYSTATHIONASE-RXN-L2R", "CYSTATHIONASE-RXN-L2R", "RXN2T-74-L2R"))
rm.frames=c('HOMOCYSTEINE-S-METHYLTRANSFERASE-RXN', 'MMUM-RXN', 'R15-RXN')
rr=unlist(apply(as.matrix(rm.frames), 1, FUN=function(x){ grep(x, model.rxns, v=T) }))
f=FBA(s.al)
f2=FBA(s.al, ko=rr[1:4], se='G')
f$x[rr]
(sr=nc.pwy$rxn[nc.pwy$PathwayID %in% c('PWY-5141', 'PWY-6151')])

model.rxns=setdiff(model.rxns, rr)

cn=check.new.media(s.al, rm='SULFATE', add='MET[CCO-EXTRACELLULAR]')
cn=check.new.media(s.al[,!(colnames(s.al) %in% rad.rxns[['NCU06558.5']])], rm='SULFATE', add='MET[CCO-EXTRACELLULAR]')

mets=c('NCU02430.5', 'NCU07001.5', 'NCU07987.5')
(gg=gpr[gpr$RXN %in% gpr$RXN[gpr$Genes %in% 'NCU05752.5'],])
rr=intersect(gg$RXN, colnames(s.al))
#'RXN-10724-L2R'
f2=FBA(s.al, ko=rr)
trace.mets(s.al, f2$x, from='SUCROSE[CCO-EXTRACELLULAR]', to='QUINOLINATE')

###metacyc rxns in model#################################################################################################
ssm=ss.mat(s.sp, grep('GLUCOSE-1-DEHYDROGENASE-', colnames(s.sp), v=T))
colnames(ssm)[6:9]
model.rxns=union(colnames(ssm)[6:9], setdiff(model.rxns, colnames(ssm)[1:4]))

ssm=ss.mat(s.sp, grep('DIHYDROOROTATE-DEHYDROGENASE-', colnames(s.sp), v=T))

###check limed-fba diff
f=FBA(s.al, ko=rad.rxns[['NCU08313.5']], goal=unique(bm$comp[bm$coeff<0]), q=F)
f=FBA(s.sp[,colnames(s.al)], ko=rad.rxns[['NCU08313.5']], goal=unique(bm$comp[bm$coeff<0]), q=F)

###ace muts
f=FBA(s.al, ko=c("PEPDEPHOS-RXN-R2L", 'ALANINE-AMINOTRANSFERASE-RXN-L2R', 'RXN-6902-R2L', '2.5.1.19-RXN-L2R'), fba.ub=fva.ub)
(tm=trace.mets(s.al, f$x, from='SUCROSE[CCO-EXTRACELLULAR]', to='PYRUVATE', pwy=nc.pwy))
f1=FBA(s.al, ko=c('ALANINE-AMINOTRANSFERASE-RXN-L2R', 'RXN-6902-R2L'), goal=unique(bm$comp[bm$coeff<0])) #'2.5.1.19-RXN-L2R'
f1$x[tm$rxns]

model.rxns=union(grep('GAPOXNPHOSPHN-RXN', colnames(s.sp), v=T), setdiff(model.rxns, grep('GAPDHSYNEC-RXN', colnames(s.al), v=T)))

##broad set#############################################################################################################
bk=names(br)
i=4
f=FBA(s.al, ko=broad.rxns[[ bk[i] ]], goals=unique(bm$comp[bm$coeff<0]), eps=10**-3, se='G')
f2=FBA(sg, ko=broad.rxns[[ bk[i] ]], goals=unique(bm$comp[bm$coeff<0]), eps=10**-3, se='E')
#setdiff(names(f2$x)[f2$x>10**-6], colnames(s.al))

new.rxns=NULL
pv=pen.v
pv[ broad.rxns[[ bk[i] ]] ]=10**6
mp <- min.pen(sp=s.all, goals=c('biomass'), pen.v=pv)
#rxns[rownames(mp),]
(mpm <- min.pen.milp(sp=s.all, pen.v=pv, rxns0=rownames(mp), goals=c('biomass'), n.add=0, ctrl=list(tilim=180)))
ss.mat(s.sp, rownames(mpm))
rxns[rownames(mpm),]

new.rxns=c(new.rxns, rownames(mpm))

###h balance#############################################################################################################
bal.mat <- read.csv('rxns_imbalance.csv', row.names=1)
bal.mat <- bal.mat[,c('H', 'C', 'O', 'N', 'P', 'S')]
unb.rxns <- rownames(bal.mat)[rowSums(bal.mat>0)>0]
i=intersect(colnames(s.al), unb.rxns)
s.al['PROTON', i]=s.al['PROTON', i]-bal.mat[i,'H']-10**-3
f=FBA(s.al, fba.ub=fva.ub)

###check exports#########################################################################################################
g=grep('CCO-EXTRA', rownames(s.al), v=T)
unique(unlist(apply(as.matrix(g), 1, FUN=function(x){ fom=flux.of.met(s.al, f$x, x); if (length(fom)>0) cat(x, ':', mean(fom), '\n') })))

(fom=flux.of.met(s.al, f$x, 'CARBON-DIOXIDE[CCO-EXTRACELLULAR]'))
ss.mat(s.al, names(fom))

###supplements###########################################################################################################
cn=check.new.media(s.al, add=c('ACET[CCO-EXTRACELLULAR]', "ACETYL-COA[CCO-MIT]"), ko=rad.rxns[['NCU06482.5']], goal=bm.goals$comp)
cn=check.new.media(s.al, add=c('ACET[CCO-EXTRACELLULAR]', 'ATP', 'CO-A'), ko=rad.rxns[['NCU06482.5']], goal='ACETYL-COA')

cn=check.new.media(s.al, add=c('L-CITRULLINE'), ko=rad.rxns[['NCU07732.5']])

###vogels################################################################################################################
model.rxns=union(setdiff(model.rxns, c("RXN2T-70-L2R")), "TRANS-RXN2T-94-OH//OH.7.-L2R[CCO-PM-FUNGI]")

###full biomass##########################################################################################################
ss.mat(s.al, 'FullBiomassComposition')
f=FBA(s.al, goal=c('ERGOSTEROL', 'CPD-12930', 'CPD1F-129'))
ff=FBA(s.al, fba.obj='FullBiomassComposition')

###biomass###############################################################################################################
bm[bm$coeff>0,]
bm.neg=bm[bm$coeff<0,]
bm.neg[bm.neg$comp %in% bm.neg$comp[duplicated(bm.neg$comp)],]

s.al['CPD2T-25', 'Lipid']=s.al['CPD1F-129', 'Lipid']=s.al['ERGOSTEROL', 'Lipid']=s.al['CPD-12930', 'CellWall']=0
s.al['FADH2', 'biomass']=-1
bm.goals.v=setdiff(union('FADH2', bm.goals$comp), c('CPD2T-25', 'CPD1F-129', 'ERGOSTEROL', 'CPD-12930'))
f0=FBA(s.al, eps=10**-5)
s.al['WATER', 'biomass']=0

i=4
(kos=intersect(br[[i]], colnames(s.al)))
f0$x[kos][f0$x[kos]>0]
f=FBA(s.al, ko=kos, goal=bm.goals.v, eps=10**-5)

###add rxns to model.rxns for growth on D-met############################################################################
cn=check.new.media(sp=s.al, rm.nuts='SULFATE', add.nuts='CPD-218[CCO-EXTRACELLULAR]')

frames=c('D-AMINO-ACID-OXIDASE-RXN', 'NADH-PEROXIDASE-RXN', 'TRANS-RXN2T-218', 'R15-RXN', 'TRANS-RXN2T-49', 'RXN-11811', 'RXN-9615')
apply(as.matrix(frames), 1, FUN=function(y) grep(pattern=y, x=colnames(s.sp), value=TRUE))

new.rxns=c("D-AMINO-ACID-OXIDASE-RXN-CPD-218/OXYGEN-MOLECULE/WATER//AMMONIA/HYDROGEN-PEROXIDE/CPD-479/PROTON.72.-L2R[CCO-MIT]",
"NADH-PEROXIDASE-RXN-L2R[CCO-MIT]", 
"TRANS-RXN2T-218-L2R[CCO-MIT]", "TRANS-RXN2T-218-R2L[CCO-MIT]",
"R15-RXN-CPD-479/L-ALPHA-ALANINE//MET/PYRUVATE.38.-L2R[CCO-MIT]",
"TRANS-RXN2T-49-CPD-218//CPD-218.17.-L2R[CCO-MIT]", "TRANS-RXN2T-49-L-ALPHA-ALANINE//L-ALPHA-ALANINE.33.-L2R[CCO-MIT]", 
"TRANS-RXN2T-49-MET//MET.9.-R2L[CCO-MIT]",
"RXN-11811-L2R[CCO-MIT]", "RXN-11811-R2L[CCO-MIT]",  
"RXN-9615-L2R[CCO-MIT]", "RXN-9615-R2L[CCO-MIT]")

model.rxns=union(model.rxns, new.rxns)

###paths.tsv2summary#####################################################################################################
##plot summary of rxns per level 3 class
#jz will do this
paths <- read.delim('Neurospora/pwy-class-subclass.tsv')
#level 1, since treat "Pathways" as level 0
paths1 <- paths[paths$ParentId=='Pathways', 'ChildId']
paths2 <- paths[paths$ParentId %in% paths1, 'ChildId']
paths3 <- paths[paths$ParentId %in% paths2, 1:2]
paths4 <- paths[paths$ParentId %in% paths3$ChildId, 1:2]
paths5 <- paths[paths$ParentId %in% paths4$ChildId, 1:2]
paths6 <- paths[paths$ParentId %in% paths5$ChildId, 1:2]
#there are no paths7

#jz named columns 'ChildId' in paths but 'ChildID' in paths.i
paths.i <- read.delim('Neurospora/pwy-class-instances.tsv')
(i3 <- paths.i[paths.i$ChildID %in% nc.pwy2$PathwayID & paths.i$ParentID %in% paths3$ChildId,])
(i4 <- paths.i[paths.i$ChildID %in% nc.pwy2$PathwayID & paths.i$ParentID %in% paths4$ChildId,])
(i5 <- paths.i[paths.i$ChildID %in% nc.pwy2$PathwayID & paths.i$ParentID %in% paths5$ChildId,])
(i6 <- paths.i[paths.i$ChildID %in% nc.pwy2$PathwayID & paths.i$ParentID %in% paths6$ChildId,])

all(nc.pwy2$PathwayID %in% paths.i$ChildID) #TRUE

###undo strip fuckup#####################################################################################################
sp=s.sp
#get generics
g=grep('//', colnames(s.al), v=T)
rv <- apply(sp[,!(colnames(sp) %in% g)], 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x/min(abs(x)), sep=':', collapse=';') })
rvg <- apply(sp[,g], 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x/min(abs(x)), sep=':', collapse=';') })
#dup rxns
rv1 <- rv[rv %in% rvg]
rvg0 <- rvg[rvg %in% rv1]

rxns.v <- c(rv1, rvg0)
rxns.annot <- rxns[names(rxns.v),]
rxns.v.o <- rxns.v[order(rxns.annot$nc, !(names(rxns.v) %in% grep('//', names(rxns.v), value=TRUE)), rxns.annot$gpr, rxns.annot$nc.pwy, rxns.annot$prob, rxns.annot$rev, -nchar(names(rxns.v)), decreasing=TRUE)]
rxns.out <- rxns.v.o[!duplicated(rxns.v.o)]
sum(names(rvg0) %in% names(rxns.out))

model.rxns <- setdiff(mr, names(rvg0))
model.rxns <- union(model.rxns, names(rxns.out))

#check defuckup
rr=rxns.al[rxns.al$nc==0,]
rn <- apply(s.sp[,rownames(rxns.al)[rxns.al$nc==1]], 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x/min(abs(x)), sep=':', collapse=';') })
rr <- apply(s.sp[,intersect(model.rxns, rownames(rxns.al)[rxns.al$nc==0])], 2, FUN=function(x){ x <- x[x!=0]; paste(names(x), x/min(abs(x)), sep=':', collapse=';') })
sum(rr %in% rn) #0

c(RXN-10965-R2L, RXN-10862-R2L)
