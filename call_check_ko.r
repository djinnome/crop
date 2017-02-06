##jmd
##1.23.12
##call_check_ko.r
##fcn to call check.ko & present results

##each element of the list ko.rxns represents a gene
call.check.ko <- function(sp, ko.rxns, ko.names=names(ko.rxns), set.name, obs=0, eps=10**-6, ...){
    cat('\n Checking', set.name, 'KOs \n')
    
    ck <- check.ko(s.test=sp, ko.lst=ko.rxns, obs=obs, ...)
    rownames(ck) <- ko.names
    
    ko.rxns.not.in.model <- sapply(ko.rxns, FUN=function(x){ length(intersect(x, colnames(sp)))==0 })
    if (sum(ko.rxns.not.in.model)>0){ cat(sum(ko.rxns.not.in.model), 'genes are not in model. \n') }
    
    cat('summary of check.ko growth:\n')
    ck <- ck[!ko.rxns.not.in.model,]
    print(summary(as.factor(ck$pred>eps)))
    
    cat('Incorrect', set.name, 'predictions:', rownames(ck)[!is.na(ck$pred) & (ck$pred>eps)!=ck$obs], '\n')
    cat('NA', set.name, 'predictions:', rownames(ck)[is.na(ck$pred)], '\n')
    cat('There are', sum(is.na(ck$pred)|(ck$pred>eps)!=ck$obs), 'predictions that are wrong or NA.\n')
    
    return(ck)
}
