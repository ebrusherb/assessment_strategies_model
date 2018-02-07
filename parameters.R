source('sub2ind.R')

# group size, category width, memory window, rho
parameters = data.frame(N=c(),rho=c(),w=c(),delta=c(),l=c())

N_vals = c(25,50,100)
rho_vals = c(seq(0.5,0.9,by=0.2),0.99)
wind_vals = c(250,seq(500,2500,by=500),Inf,5000)
delta_vals = c(seq(0,1,by=0.25))
learn_rate_vals = c(0.05,0.1,0.25,0.5)

xN = length(N_vals)
xrho = length(rho_vals)
xwind = length(wind_vals)
xdelta = length(delta_vals)
xlearn = length(learn_rate_vals)

for(l in 1:xrho){
	for(k in 1:xwind){
		for(j in 1:xdelta){		
			for(m in 1:xlearn){
				for(i in 1:xN){
					ind = sub2ind(c(xN,xdelta,xwind,xrho,xlearn),c(i,j,k,l,m))
					parameters = rbind(parameters,data.frame(N=N_vals[i],rho=rho_vals[l],w=wind_vals[k],delta=delta_vals[j],l=learn_rate_vals[m]))
				}
			}
		}
	}
}

# 
breaks = seq(1,dim(parameters)[1],by=xN*xdelta*xlearn)
num_breaks = length(breaks)
chunk = list()
if(num_breaks>1){
	for(i in 1:(num_breaks-1)){
		chunk[[i]] = breaks[i]:(breaks[i+1]-1)
	}
}
chunk[[num_breaks]] = breaks[num_breaks]:dim(parameters)[1]

first = lapply(chunk,function(x) x[1])
first = unlist(first)
