##jmd
##4.1.11
##cut_mets.r
##make matrix: beta_in <= sum(beta_out) & vice versa, which enforces gaplessnes on metabolites

require('Matrix')

##return only out.minus.sum.in if sense='G'
cut.mets <- function(s, n.top.mets.rm=75, sense='G'){
 ##don't want many constraints for promiscuous cofactors
 rs <- rowSums(s!=0)
 names(rs) <- rownames(s)
 s <- s[rank(-rs)>n.top.mets.rm,]
 
 ##split by sign
 s.plus <- s>0
 s.minus <- s<0

 ##make expanded matrix w/ 1 beta.in (or beta.out) coeff per line
 s.plus.e <- sparseMatrix(i=1:sum(s.plus), j=which.col(s>0), x=1, dims=c(sum(s.plus), ncol(s)))
 s.minus.e <- sparseMatrix(i=1:sum(s.minus), j=which.col(s<0), x=1, dims=c(sum(s.minus), ncol(s)))

 ##make sum matrix w/ rows matching expanded s *of opposite sign*
 #a 1 for each positive (or negative) element of s
 sigma.plus <- s.plus[rep(1:nrow(s.plus), times=rowSums(s.minus)),]
 sigma.minus <- s.minus[rep(1:nrow(s.minus), times=rowSums(s.plus)),]

 ##combine
 in.minus.sum.out <- s.plus.e-sigma.minus
 out.minus.sum.in <- s.minus.e-sigma.plus

 if (sense=='G'){ ret <- out.minus.sum.in } else { ret <- rBind(in.minus.sum.out, out.minus.sum.in) }
 if (any(duplicated(rownames(ret)))){ rownames(ret) <- paste(1:nrow(ret), rownames(ret), sep='_') }
 
 return(ret)
}
