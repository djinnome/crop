##jmd
##3.7.11
##check_fba.r

source('fba.r')

##run fba
#fba grows on vogel's, dies when rm source of {c,n,s, or p}
#growth insensitive to ub=1000 or ub=10^9
check.fba <- function(a, ko.rxns.lst=NULL, sense='E'){
	pred <- rep(NA, length(ko.rxns.lst))
	for (ko.rxns.ind in 1:length(ko.rxns.lst)){
    		ko.rxns <- ko.rxns.lst[[ko.rxns.ind]]
		if (sum(ko.rxns %in% colnames(a))>0){ 
			pred[ko.rxns.ind] <- fba(a=a, ko.rxns=ko.rxns, sense=sense)$lp$obj 
		}
	}
	return(pred)
}

##run dual
#growth = fba, except when ko.rxns=sugar
check.dual <- function(a, ko.rxns.lst){
	pred.dual <- rep(NA, length(ko.rxns.lst))
	for (ko.rxns.ind in 1:length(ko.rxns.lst)){
    		ko.rxns <- ko.rxns.lst[[ko.rxns.ind]]
		if (sum(ko.rxns %in% colnames(a))>0){
			pred.dual[ko.rxns.ind] <- fba.dual(a=s.sp, ko.rxns=ko.rxns)$lp$obj 
		}
	}
	return(pred.dual)
}

###check#################################################################################################################
##check a met
#met <- 'POLYMER-INST-Peptidoglycans-C312-H512-N64-O152'
#in.rxns <- s.sp[met,]!=0
#(ss <- s.sp[s.sp[,in.rxns]!=0, in.rxns])
