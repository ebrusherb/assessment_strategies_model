sub2ind <- function(v,sub){
	P = prod(v)
	m = array(1:P,dim=v)
	ind = m[matrix(sub,nrow=1)]
	return(ind)
}