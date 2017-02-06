##jeremy zucker
##7.20.2011
##read_fba.r

require('Matrix')

read.flux <- function(flux.txt) {
  flux.mat <- read.delim(flux.txt)
  v <- flux.mat$flux
  names(v) <- flux.mat$rxn
  return(v)
}

read.S <-  function(S.txt) {
  s <- read.delim(S.txt, as.is=FALSE, na='NIL')
  s.sp <- sparseMatrix(i=as.numeric(s$compound), j=as.numeric(s$rxn), x=as.real(s$coeff))
  dimnames(s.sp) <- list(levels(s$compound), levels(s$rxn))
  return(s.sp)
}
