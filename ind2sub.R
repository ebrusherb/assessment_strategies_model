ind2sub <-function(v,ind){
	# v=c(2,3)
	# ind=5
	num_dim = length(v)
	sub <- array(0,dim=c(1,num_dim))
	P = cumprod(v)
	for(l in num_dim:2){ #work backwards to fill in from the last dimension
		ind = (ind-1)%%P[l]+1 # index in current dimension
		div = P[l-1] #size of current dimension
		sub[l] = floor((ind-1)/div)+1 #number of times ind-1 goes into product of higher dimensions + 1
		}
	sub[1] = (ind-1)%%div+1 #left over after ind-1 goes into number or rows
	return(sub)
}