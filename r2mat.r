##jmd
##6.19.12
##r2mat.r

r2mat <- function(rfile, loc="C:/Documents and Settings/jdreyf/Desktop/farm"){
    rcode <- scan(paste(loc, '/', rfile, sep=''), what='character', sep='@')
        
    mcode <- r2mat.regexp(rcode)
    
    mfile <- gsub('r$', 'm', rfile, ignore.case=TRUE)
    write(mcode, paste(loc, '/matlab/', mfile, sep=''))
    return(paste('written to ', loc, '/matlab/', mfile, sep=''))
}

r2mat.regexp <- function(rcode){
mcode <- 
    gsub('#', '%', fixed=TRUE,
    gsub('[', '(', fixed=TRUE, 
    gsub(']', ')', fixed=TRUE, 
    gsub('!', '~', fixed=TRUE, 
    gsub('<-', '=', fixed=TRUE,
    gsub('$', '.', fixed=TRUE,
    gsub('\\[\\[(.+)\\]\\]', '\\{\\1\\}',
    gsub('.', '_', fixed=TRUE,
    gsub('$', ';',
    gsub('{', '', fixed=TRUE,
    gsub('}', 'end', fixed=TRUE,
    
    gsub('NA', 'NaN',
    gsub('NULL', '{}', fixed=TRUE,  
    gsub('for \\((.+) in (.+)\\).*\\{', 'for \\1=\\2',
    gsub('while \\((.+) in (.+)\\).*\\{', 'for \\1=\\2',
    gsub('(.+) <- function\\((.+)\\)\\{', 'function \\[ret\\]=\\1\\(\\2\\)',
    
    #handle vec/matrix constructions, incl. some from 'Matrix' pkg
    #ignore rbind for now, since it requires gsub mult. commas into semicolons
    gsub('c\\((.+)\\)', '\\[ \\1 \\]', 
    gsub('cbind\\((.+)\\)', '\\[ \\1 \\]', 
    gsub('cBind\\((.+)\\)', 'sparse\\(\\[ \\1 \\]\\)', 
    
    gsub('matrix\\(0,', 'zeros\\(', ignore.case=TRUE,
    gsub('matrix\\(1,', 'ones\\(', ignore.case=TRUE,
    gsub('diag\\(', 'eye\\(', ignore.case=TRUE,
    gsub('nrow\\((.+)\\)', 'size\\(\\1,1\\)',    
    gsub('ncol\\((.+)\\)', 'size\\(\\1,2\\)',
    
    gsub('dim(', 'size(', fixed=TRUE,    
    gsub('which(', 'find(', fixed=TRUE,
    gsub('grep\\((.+),(.+)\\)', 'regexp\\(\\2,\\1\\)',
    gsub('list(', 'struct(', fixed=TRUE,
    gsub('stop(', 'error(', fixed=TRUE,
    gsub('as.numeric(', 'str2double(', fixed=TRUE,
    gsub('print(', 'disp(', fixed=TRUE,
    gsub('setwd\\((.+)\\)', 'cd \\1',
    gsub('next', 'continue',
    gsub('is.na(', 'is.nan(', fixed=TRUE,
    gsub('paste\\((.+)\\)', '\\[\\1\\]',
    
    # when want all rows or columns of a matrix/df
    # there could be a space in the R code
    gsub('\\[( ),', '\\(:,',
    gsub(',( )\\]', ',:\\)',
    
    rcode
    ))))))))))))))))))))))))))))))))))
    return(mcode)
}

## to do
# df=struct w cell array & dimnames
# named.vec=struct w 1d cell array or matrix & names
# rbind(a,b)=[a;b]

# should manipulate mm[rownames,1:2], but I don't know the name of the object in the structure that holds the row names
# fcn parameter passing, eg
#   m=matrix(data=0, nrow=4, ncol=2, dimnames=list(1:2,1:4)) &
#   new.mat <- function(data=0, nrow=1, ncol=1, dimnames){
# apply
# matrix operations: * to .*, %*% to *
# fcn return
