library(lattice)
library(pracma)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(digest)
library(RColorBrewer)
setwd('/Users/eleanorbrush/Dropbox/evo_badgesVSrecognition/code_extended_model')
source('get_legend.R')
load('noise=0.1/all_parameters.Rdata')
source('sub2ind.R')
source('ind2sub.R')

parameters_cat_to_plot = parameters_cat
parameters_cat_to_plot$w[parameters_cat$w==Inf] = 3000

parameters_rule_to_plot = parameters_rule
parameters_rule_to_plot$w[parameters_rule$w==Inf] = 3000

Tfights = fixed_parameters$Tfights
down_sample = fixed_parameters$down_sample
time = seq(1,fixed_parameters$Tfights,by=fixed_parameters$down_sample)
xN_cat = length(unique(parameters_cat$N))
xdelta_cat = length(unique(parameters_cat$delta))
xl_cat = length(unique(parameters_cat$l))
xw_cat = length(unique(parameters_cat$w))
xrho_cat = length(unique(parameters_cat$rho))
dim_cat = data.frame(N=xN_cat,rho=xrho_cat,w=xw_cat,delta=xdelta_cat,l=xl_cat)
xN_rule = length(unique(parameters_rule$N))
xdelta_rule = length(unique(parameters_rule$delta))
xl_rule = length(unique(parameters_rule$l))
xw_rule = length(unique(parameters_rule$w))
xrho_rule = length(unique(parameters_rule$rho))
dim_rule = data.frame(N=xN_rule,rho=xrho_rule,w=xw_rule,delta=xdelta_rule,l=xl_rule)

find <- function(N=NA,rho=NA,w=NA,delta=NA,l=NA,c=NA,p=parameters_cat){
	params = unlist(as.list(environment())[1:6])
	provided = which(!is.na(params))
	vec = 1:dim(p)[1]
	for(i in provided){
		now = names(params)[i]
		vec = intersect(vec,which(p[[now]]==params[i]))
	}
	return(vec)
}

find_df <- function(params,p=parameters_cat){
	to_return = NULL
	for(j in 1:dim(params)[1]){	
		provided = which(!is.na(params[j,]))
		vec = 1:dim(p)[1]
		for(i in provided){
			now = names(params)[i]
			vec = intersect(vec,which(p[[now]]==params[[now]][j]))
		}
	if(sum(dim(p)==dim(parameters_cat))==2){
		if(sum(p==parameters_cat)==prod(dim(p)) && length(vec)<prod(dim_cat)/prod(dim_cat[provided])){
			vec = c(vec,rep(NA,prod(dim_cat)/prod(dim_cat[provided])-length(vec)))
		}
	}
	if(sum(dim(p)==dim(parameters_rule))==2){
		if(sum(p==parameters_rule)==prod(dim(p)) && length(vec)<prod(dim_rule)/prod(dim_rule[provided])){
			vec = c(vec,rep(NA,prod(dim_rule)/prod(dim_rule[provided])-length(vec)))
		}
	}
	to_return = rbind(to_return,vec)
	}
	return(to_return)
}

myround <-function(x,digits){
	round(x*10^digits)/10^digits
}

moving_average <- function(v,k=5,keep_end=TRUE){
	v2 = array(NA,dim=c(length(v),1))
	if(keep_end){
		M = length(v)
	}else{
		M = length(v)-floor((k-1)/2)
	}
	for(i in 1:M){
		v2[i] = mean(v[i+(max(-floor(k/2),0):min(floor((k-1)/2),length(v)))],na.rm=TRUE)
	}
	return(v2)
}

moving_var <- function(v,k=5,keep_end=TRUE){
	v2 = array(NA,dim=c(length(v),1))
	if(keep_end){
		M = length(v)
	}else{
		M = length(v)-floor((k-1)/2)
	}
	for(i in 1:M){
		v2[i] = var(v[i+(max(-floor(k/2),0):min(floor((k-1)/2),length(v)))],na.rm=TRUE)
	}
	return(v2)
}

find_changes <- function(s){
	intersect(which(diff(s)!=0),intersect(which(c(diff(s)[2:(length(s)-1)],0)==0),which(c(0,diff(s)[1:(length(s)-2)])==0)))+1
}

find_changes_immediate <- function(s){
	which(diff(s)!=0)+1
}

optimization<- function(w){	
	cost_cat = apply(error_cat,2,function(v) 2*w*v+(1-w)*time/30001)
	
	best_cat = array(NA,dim(parameters_cat)[1])
	start = 3
	start2 = 10
	start3 = 20
	noisy  = array(NA,dim(parameters_cat)[1])
	for(i in 1:dim(parameters_cat)[1]){
		d = c(rep(NA,start-1),diff(cost_cat[start:length(time),i]))
		k = 70
		d_avg = moving_average(d,k)
		# d_var = moving_var(d,10)
		d_var = moving_var(error_cat[start:length(time),i],10)
		max_var = max(d_var[start3:length(d_var)],na.rm=TRUE)
		M = diff(range(error_cat[start2:length(time),i],na.rm=TRUE))
		# if((diff(range(error_cat[start2:dim(cost_cat)[1],i],na.rm=TRUE)))>0.15){
		if(M>0.075 && max_var<0.00006*M){
			offset = 4
			s = sign(myround(d_avg,4.2))
			# }else if((diff(range(error_cat[start2:dim(cost_cat)[1],i],na.rm=TRUE)))>0.1){
			noisy[i] = FALSE			
			}else if((M>0.1 && max_var<0.0004*M && w<1) || max_var<0.00018*M){				
			offset = 4
			s = sign(myround(d_avg,3.3))
			noisy[i] = TRUE			
			}else{
			offset = 8
			s = sign(myround(d_avg,3.1))	
			noisy[i] = TRUE
			}
		# if((noisy && !(parameters_cat$N[i]==25 && parameters_cat$w[i]>=2500)) || w==1){
			s[s==0] = 1		
		# }else{
			# s[s==0] = -1
		# }
		z = apply(matrix(1:length(s),ncol=1),1,function(x) identical(s[(max(0,x-offset):min(length(s),x+offset))],rep(s[x],offset*2+1)))
		if(length(s[z])>1){
			f = which(z)[find_changes(s[z])]-offset
		}else{
			f = NULL
		}	
		if(length(f)==0 || length(f)>3){
			if(sum(z)==0 || max(s[z],na.rm=TRUE)>(-1)){
				best_cat[i] = start2
			} else{
				best_cat[i] = length(time)
			}
		}else {
			if(max(s[f])==-1){
				best_cat[i] = length(time)
			}else{
				best_cat[i] = min(f[which(s[f]>(-1))])
			}
		}
	}
	
	optimized_cat = apply(matrix(c(best_cat,1:dim(parameters_cat)[1]),nrow=2,byrow=TRUE),2,function(x) cost_cat[x[1],x[2]])
	
	# best_cat[intersect(which(best_cat==start2),which(noisy))] = 1
	
	cost_rule = apply(error_rule,2,function(v) 2*w*v+(1-w)*time/30001)
	
	best_rule = array(NA,dim(parameters_rule)[1])
	start = 3
	start2 = 10
	start3 = 20
	for(i in 1:dim(parameters_rule)[1]){
		d = c(rep(NA,start-1),diff(cost_rule[start:length(time),i]))
		k = 5
		d_avg = moving_average(d,k,keep_end=FALSE)[1:(length(d)-floor((k-1)/2))]	
		d_var = moving_var(error_rule[start:length(time),i],10)	
		offset = 5
		s = sign(myround(d_avg,3))	
		if(w<0.8){
			offset = 6
			s = sign(myround(d_avg,3.2))
		}else if(max(d_var[start3:length(d_var)],na.rm=TRUE)<4e-05){				
			offset = 6
			s = sign(myround(d_avg,3.1)) ##changed from 3.1, which was changed from 3.15
		}else{
			offset = 10
			s = sign(myround(d_avg,3))	
		}
		z = apply(matrix(1:length(s),ncol=1),1,function(x) identical(s[(max(0,x-offset):min(length(s),x+offset))],rep(s[x],offset*2+1)))
		if(length(s[z])>1){
			f = which(z)[find_changes_immediate(s[z])]-offset
		}else{
			f = NULL
		}	
		if(length(f)==0 || length(f)>3){
			if(sum(z)==0 || max(s[z],na.rm=TRUE)>(-1)){
				best_rule[i] = start2
			} else{
				best_rule[i] = length(time)
			}
		}else {
			if(max(s[f])==-1){
				best_rule[i] = length(time)
			}else{
				best_rule[i] = min(f[which(s[f]>(-1))])
			}
		}
	}
	
	optimized_rule = apply(matrix(c(best_rule,1:dim(parameters_rule)[1]),nrow=2,byrow=TRUE),2,function(x) cost_rule[x[1],x[2]])

	# best_rule[best_rule==start2] = 1
	return(list(best_cat=best_cat,optimized_cat=optimized_cat,best_rule=best_rule,optimized_rule=optimized_rule))
}

optimization_cat<- function(error,w){
	if(length(error)>length(time)){
	N = error[length(time)+1]
	wind = error[length(time)+2]
	error = error[1:length(time)]	
	}else{
		N = 25
		wind = Inf
	}
	cost = 2*w*error+(1-w)*time/30001
		
	start = 3
	start2 = 10
	start3 = 20
	noisy  = NA
	best_cat = NA

	d = c(rep(NA,start-1),diff(cost[start:length(cost)]))
	k = 70
	d_avg = moving_average(d,k)
	# d_var = moving_var(d,10)
	d_var = moving_var(error[start:length(error)],10)
	max_var = max(d_var[start3:length(d_var)],na.rm=TRUE)
	M = diff(range(error[start2:length(error)],na.rm=TRUE))
	# if((diff(range(error_cat[start2:dim(cost_cat)[1],i],na.rm=TRUE)))>0.15){
	if(M>0.075 && max_var<0.00006*M){
		offset = 4
		s = sign(myround(d_avg,4.2))
		# }else if((diff(range(error_cat[start2:dim(cost_cat)[1],i],na.rm=TRUE)))>0.1){
		noisy = FALSE		
		}else if((M>0.1 && max_var<0.0004*M && w<1) || max_var<0.00018*M){	###changed to 0.9 from 1		
		offset = 4
		s = sign(myround(d_avg,3.3))
		noisy = TRUE		
		}else{
		offset = 8
		s = sign(myround(d_avg,3.1))	
		noisy = TRUE		
		}
	# if((noisy && !(N==25 && wind>=2500)) || w>=0.8){
		s[s==0] = 1		
	# }else{
		# s[s==0] = -1
	# }			
	z = apply(matrix(1:length(s),ncol=1),1,function(x) identical(s[(max(0,x-offset):min(length(s),x+offset))],rep(s[x],offset*2+1)))
	if(length(s[z])>1){
		f = which(z)[find_changes(s[z])]-offset
	}else{
		f = NULL
	}	
	if(length(f)==0 || length(f)>3){
		if(sum(z)==0 || max(s[z],na.rm=TRUE)>(-1)){
			best_cat = start2
		} else{
			best_cat = length(error)
		}
	}else {
		if(max(s[f])==-1){
			best_cat = length(error)
		}else{
			best_cat = min(f[which(s[f]>(-1))])
		}
	}	
	
	optimized_cat = cost[best_cat]
	
	#if(best_cat==start2 && noisy){best_cat=1} 
	
	return(list(best_cat=best_cat,optimized_cat=optimized_cat))
}	



optimization_rule<- function(error,w){		
	cost = 2*w*error+(1-w)*time/30001
	
	start = 3
	start2 = 10
	start3 = 20
	best_rule = NA
	d = c(rep(NA,start-1),diff(cost[start:length(cost)]))
	k = 5
	d_avg = moving_average(d,k,keep_end=FALSE)[1:(length(d)-floor((k-1)/2))]	
	d_var = moving_var(error[start:length(error)],10)	
	offset = 5
	s = sign(myround(d_avg,3))	
	if(w<0.8){
			offset = 6
			s = sign(myround(d_avg,3.2))
		}else if(max(d_var[start3:length(d_var)],na.rm=TRUE)<4e-05){				
			offset = 6
			s = sign(myround(d_avg,3.1)) ##changed from 3.1, which was changed from 3.15
		}else{
			offset = 10
			s = sign(myround(d_avg,3))	
		}
	z = apply(matrix(1:length(s),ncol=1),1,function(x) identical(s[(max(0,x-offset):min(length(s),x+offset))],rep(s[x],offset*2+1)))
	if(length(s[z])>1){
		f = which(z)[find_changes_immediate(s[z])]-offset
	}else{
		f = NULL
	}	
	if(length(f)==0 || length(f)>3){
		if(sum(z)==0 || max(s[z],na.rm=TRUE)>(-1)){
			best_rule = start2
		} else{
			best_rule = length(time)
		}
	}else {
		if(max(s[f])==-1){
			best_rule = length(time)
		}else{
			best_rule = min(f[which(s[f]>(-1))])
		}
	}
	
	optimized_rule = cost[best_rule]

	# best_rule[best_rule==start2] = 1
	
	return(list(best_rule=best_rule,optimized_rule=optimized_rule))
}

####### plotting parameters
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

lwd=2
squish = log
# unsquish = function(x){x^2}
unsquish = exp

#### palettes 

set1cols=brewer.pal(5,'Set1')
# warm_pal = brewer.pal(8,'Spectral')[-5]
red_pal = rev(brewer.pal(9,'Reds'))
blue_pal = rev(brewer.pal(9,'Blues')[2:9])
two_blues = brewer.pal(9,'Blues')[c(4,6)]
div_pal = brewer.pal(11,'RdBu')
red_blue_purp = c(brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[c(4,6,7)],brewer.pal(9,'Greys')[5])
purp_pal = brewer.pal(9,'Purples')[c(5,7,8)]
green_pal = brewer.pal(9,'Greens')[c(3,5,7,8)]


#### optimization over nearly all parameters
examples_optimized = data.frame(N=rep(c(25,50,100),each=(xw_cat)*xrho_cat),rho=rep(rep(c(0.5,0.9,0.99),each=xw_cat),times=xN_cat),w=rep(c(250,500,1000,1500,2000,2500,Inf),times=xN_cat*xrho_cat))

# w_vec = seq(0.5,1,by=0.1)
# f_cat = as.vector(find_df(examples_optimized,p=parameters_cat))
# f_rule = as.vector(find_df(examples_optimized,p=parameters_rule))

# optimal_cognition_cat = data.frame(c=c(),N=c(),rho=c(),w=c(),delta=c(),l=c(),cost=c(),time=c(),error=c())
# optimal_cognition_rule = data.frame(c=c(),N=c(),rho=c(),w=c(),delta=c(),l=c(),cost=c(),time=c(),error=c())

# for(j in 1:length(w_vec)){
	# w = w_vec[j]
	# opt = apply(rbind(error_cat[,f_cat],parameters_cat$N[f_cat],parameters_cat$w[f_cat]),2,optimization_cat,w=w)
	# optimized_cat = unlist(lapply(opt,function(x) x$optimized_cat))
	# best_cat = unlist(lapply(opt,function(x) x$best_cat))
	# opt = apply(error_rule[,f_rule],2,optimization_rule,w=w)
	# optimized_rule = unlist(lapply(opt,function(x) x$optimized_rule))
	# best_rule = unlist(lapply(opt,function(x) x$best_rule))
	# for(i in 1:dim(examples_optimized)[1]){
		# f = find(N=examples_optimized$N[i],rho=examples_optimized$rho[i],w=examples_optimized$w[i],p=parameters_cat[f_cat,])
		# best = f[which.min(optimized_cat[f])]
		# optimal_cognition_cat = rbind(optimal_cognition_cat,data.frame(c=w,N=examples_optimized$N[i],rho=examples_optimized$rho[i],w=examples_optimized$w[i],delta=parameters_cat$delta[f_cat[best]],l=parameters_cat$l[f_cat[best]],cost=optimized_cat[best],time=best_cat[best],error=error_cat[best_cat[best],f_cat[best]]))
		# f = find(N=examples_optimized$N[i],rho=examples_optimized$rho[i],w=examples_optimized$w[i],p=parameters_rule[f_rule,])
		# if(!isempty(f)){
			# optimal_cognition_rule = rbind(optimal_cognition_rule,data.frame(c=w,N=examples_optimized$N[i],rho=examples_optimized$rho[i],w=examples_optimized$w[i],delta=0,l=0.05,cost=optimized_rule[f],time=best_rule[f],error=error_rule[best_rule[f],f_rule[f]]))
		# }
	# }
# }

# optimal_cognition_gen = data.frame(c=c(),N=c(),rho=c(),w=c(),delta=c(),l=c(),cost=c(),time=c(),error=c())

# for(j in 1:length(w_vec)){
	# w = w_vec[j]
	# for(i in 1:dim(examples_optimized)[1]){		
		# f = find(N=examples_optimized$N[i],rho=examples_optimized$rho[i],p=parameters_genetic)
		# optimal_cognition_gen = rbind(optimal_cognition_gen,data.frame(c=w,N=examples_optimized$N[i],rho=examples_optimized$rho[i],w=examples_optimized$w[i],cost=w*2*error_genetic[1,f],time=1,error=error_genetic[1,f]))
	# }
# }

# optimal_cognition_rule$w[which(optimal_cognition_rule$w==250)] = 0
# optimal_cognition_cat$w[which(optimal_cognition_cat$w==250)] = 0
# optimal_cognition_cat$w[which(optimal_cognition_cat$w==Inf)] = 3000
# optimal_cognition_gen$w[which(optimal_cognition_gen$w==250)] = 0
# optimal_cognition_gen$w[which(optimal_cognition_gen$w==Inf)] = 3000

# optimal_cognition_rule_expanded = data.frame(c=c(),N=c(),rho=c(),w=c(),delta=c(),l=c(),cost=c(),time=c(),error=c())
# hold = optimal_cognition_rule[find(w=2500,p=optimal_cognition_rule),]
# hold$w = 3000
# for(i in 1:length(hold$w)){
	# optimal_cognition_rule_expanded = rbind(optimal_cognition_rule_expanded,optimal_cognition_rule[(1:6)+(i-1)*6,],hold[i,])
# }
# optimal_cognition_rule = optimal_cognition_rule_expanded
# rm(optimal_cognition_rule_expanded)

# # save(optimal_cognition_rule,optimal_cognition_cat,optimal_cognition_gen,file='optimal_cognition.Rdata')

load('optimal_cognition.Rdata')
w_vec = unique(optimal_cognition_cat$c)

#### error and stopping time as a function of c

examples = data.frame(N=c(50,25,25),rho=c(0.5,0.9,0.99),w=c(2500,2500,250),delta=c(0,0,0),l=c(0.25,0.25,0.25))
examples = data.frame(N=c(25,25),rho=c(0.99,0.9),w=c(2500,2500),delta=c(0,0),l=c(0.5,0.5))
examples_delta = examples
examples_delta$delta = 0.25
num_ex = dim(examples)[1]
rule_relevant = intersect(which(names(examples)!='delta'),which(names(examples)!='l'))
gen_relevant = rule_relevant[which(names(examples)[rule_relevant]!='w')]

w_vec = seq(0,1,by=0.1)
f_recog = as.vector(find_df(examples[1,],p=parameters_cat))
f_cat = as.vector(find_df(examples_delta,p=parameters_cat))
f_rule = as.vector(find_df(examples[,rule_relevant],p=parameters_rule))

opt_c_recog = array(NA,length(w_vec))
opt_c_cat = array(NA,dim=c(dim(examples)[1],length(w_vec)))
opt_c_rule = array(NA,dim=c(dim(examples)[1],length(w_vec)))
best_c_recog = array(NA,length(w_vec))
best_c_cat = array(NA,dim=c(dim(examples)[1],length(w_vec)))
best_c_rule = array(NA,dim=c(dim(examples)[1],length(w_vec)))

for(j in 1:length(w_vec)){
	w = w_vec[j]
	opt = optimization_cat(c(error_cat[,f_recog],parameters_cat$N[f_recog],parameters_cat$w[f_recog]),w)
	opt_c_recog[j] = error_cat[opt$best_cat,f_recog]
	best_c_recog[j] = opt$best_cat
	opt = apply(rbind(error_cat[,f_cat],parameters_cat$N[f_cat],parameters_cat$w[f_cat]),2,optimization_cat,w=w)
	best_cat = unlist(lapply(opt,function(x) x$best_cat))
	opt_c_cat[,j] = apply(matrix(1:length(f_cat),nrow=1),2,function(i){error_cat[best_cat[i],f_cat[i]]})
	best_c_cat[,j] = best_cat
	opt = apply(error_rule[,f_rule],2,optimization_rule,w=w)
	best_rule = unlist(lapply(opt,function(x) x$best_rule))	
	opt_c_rule[,j] = apply(matrix(1:length(f_rule),nrow=1),2,function(i){error_rule[best_rule[i],f_rule[i]]})
	best_c_rule[,j] = best_rule
}


#### run optimization 
w=1
# opt_no_time_costs = optimization(w)

w=0.9
# opt_small_time_costs = optimization(w)
load('/Users/eleanorbrush/Dropbox/evo_badgesVSrecognition/code_extended_model/opt_data.Rdata')

best_cat = opt_no_time_costs$best_cat
optimized_cat = opt_no_time_costs$optimized_cat/2
best_rule = opt_no_time_costs$best_rule
optimized_rule = opt_no_time_costs$optimized_rule/2
best_gen = rep(squish(10),dim(parameters_genetic)[1])
optimized_gen = error_genetic[1,]

##### learning curves based on w, delta, l

marg = c(0.4,0.4,0.07,0.05)
omarg = c(0.4,0.3,0.7,0.7)

width = 6.85
height =4


N = 25
rho = 0.99
delta = 0
l = 0.05
w = 2500

ylim= c(0,0.5)
xlim=c(-0.1,30000)

pdf('/Users/eleanorbrush/Desktop/speed_accuracy_tradeoff.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)

layout(matrix(c(1,2,3,4),nrow=2,byrow=T))

for(w in c(2500,500)){
	f = find(N=N,rho=rho,delta=delta,w=w)
plot(time,time,ylim=ylim,xlim=xlim,t='n',xlab='',ylab='',xaxt='n',yaxt='n',xaxs='i',yaxs='i')
axis(1,mgp=c(3,0.5,0))
axis(2,mgp=c(3,0.5,0))
f_rule = find(N=N,rho=rho,w=w,p=parameters_rule)
lines(time,error_rule[,f_rule],t='l',col=two_blues[1],lwd=lwd)
points(time[best_rule[f_rule]],optimized_rule[f_rule],col=two_blues[1],lwd=lwd,pch=19)
f_gen = find(N=N,rho=rho,p=parameters_genetic)
lines(time,error_genetic[,f_gen],t='l',col=two_blues[2],lwd=lwd)
# points(500,error_genetic[1,f_gen],col=two_blues[2],lwd=lwd,pch=19)
for(i in 1:length(f)){
	lines(time,error_cat[,f[i]],t='l',col=red_pal[i+2],lwd=lwd)
	points(time[best_cat[f[i]]],optimized_cat[f[i]],lwd=lwd,col=red_pal[i+2],pch=19)
}
# mtext('Time, t',1,line=1.7,at=mean(time),cex=largefontsize/smallfontsize)
mtext(paste('w = ',w,sep=''),3,line=0.1,cex=largefontsize/smallfontsize)
w_leg = c()
for(i in 1:length(f)){	
	w_leg = c(w_leg,bquote(r*.(paste(' = ',unique(parameters_cat$l[f])[i],sep=''))))
	}
w_leg = sapply(w_leg,as.expression)
if(w==2500){legend(18000,0.52,legend=w_leg,col=red_pal[(1:length(f))+2],lty=1,lwd=lwd,bty='n')
	legend(6000,0.52,legend=c('Learned rule','Innate rule'),col=two_blues,lty=1,lwd=lwd,bty='n')
	mtext(expression(paste('Error, ',epsilon,sep='')),2,line=1.3,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext('(t)',2,line=1.5,at=mean(ylim)+0.12,cex=1)
	}
}

for(w in c(2500,500)){
f = find(N=N,rho=rho,l=l,w=w)[-4]
plot(time,time,ylim=ylim,xlim=xlim,t='n',xlab='',ylab='',xaxt='n',yaxt='n',xaxs='i',yaxs='i')
axis(1,mgp=c(3,0.5,0))
axis(2,mgp=c(3,0.5,0))
f_rule = find(N=N,rho=rho,w=w,p=parameters_rule)
lines(time,error_rule[,f_rule],t='l',col=two_blues[1],lwd=lwd)
points(time[best_rule[f_rule]],optimized_rule[f_rule],col=two_blues[1],lwd=lwd,pch=19)
f_gen = find(N=N,rho=rho,p=parameters_genetic)
lines(time,error_genetic[,f_gen],t='l',col=two_blues[2],lwd=lwd)
# points(500,error_genetic[1,f_gen],col=two_blues[2],lwd=lwd,pch=19)
for(i in 1:length(f)){
	lines(time,error_cat[,f[i]],t='l',col=c(red_pal[3],purp_pal)[i],lwd=lwd)
	points(time[best_cat[f[i]]],optimized_cat[f[i]],lwd=lwd,col=c(red_pal[3],purp_pal)[i],pch=19)
}
mtext(expression(paste('Error, ',epsilon,sep='')),2,line=1.3,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext('(t)',2,line=1.5,at=mean(ylim)+0.12,cex=1)
mtext('Time, t',1,line=1.7,at=mean(time),cex=largefontsize/smallfontsize)
w_leg = c()
for(i in 1:length(f)){	
	w_leg = c(w_leg,bquote(delta*.(paste(' = ',unique(parameters_cat$delta[f])[i],sep=''))))
	}
w_leg = sapply(w_leg,as.expression)
if(w==2500){legend(18000,0.52,legend=w_leg,col=c(red_pal[3],purp_pal),lty=1,lwd=lwd,bty='n')}}

dev.off()


### effect of parameters


marg = c(0.3,0.3,0.04,0.05)
omarg = c(0.7,0.65,0.65,0.05)

width = 6.85
height =4

legend_labels=c()

for(i in 1:num_ex){
	# rest_of_label = paste('=',examples$rho[i],', N=',examples$N[i],', w=',examples$w[i],sep='')
	rest_of_label = paste(' = ',examples$rho[i],sep='')
	legend_labels = c(legend_labels,bquote(rho*.(rest_of_label)))
	}
legend_labels=c('Recognition','Categorization','Learned rule','Innate rule',legend_labels)
legend_labels = sapply(legend_labels,as.expression)

# parameters_to_plot=c('rho','N','w','l','delta')
# labels = c(expression(paste('Correlation, ',rho,sep='')),'Group size, N','Memory window, w','Updating rate, r',expression(paste('Category width, ',delta,sep='')))

parameters_to_plot=c('rho','N','w')
labels = c(expression(paste('Correlation, ',rho,sep='')),'Group size, N','Memory window, w','Cost of errors, c')

time_ylim = squish(c(400,5001))

pdf('/Users/eleanorbrush/Desktop/parameters_exploration.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)

layout(matrix(1:8,nrow=2))

for(var in parameters_to_plot){
	examples_var = examples
	examples_var[[var]] = NA
	f_cat = find_df(examples_var)
	f_rule = find_df(examples_var[,rule_relevant],p=parameters_rule)
	f_gen = find_df(examples_var[,gen_relevant],p=parameters_genetic)
	examples_var = examples_delta
	examples_var[[var]] = NA
	f_cat2 = find_df(examples_var)	
	
	plot(parameters_cat_to_plot[[var]],squish(time[best_cat]),t='n',xlab='',ylab='',ylim=time_ylim,yaxt='n',xaxt='n')
	if(var=='w'){
	axis(1,at=seq(500,2500,by=1000),labels=seq(500,2500,by=1000),mgp=c(3,0.5,0))
	axis(1,at=3000,labels='All',mgp=c(3,0.5,0))
	}else{
		axis(1,mgp=c(3,0.5,0),at=unique(parameters_cat[[var]]),labels=unique(parameters_cat[[var]]))}
	ticks = c(1,1000,4000,9000,16000,25000)
	ticks = c(500,1000,2000,4000,16000)
	axis(2,at=squish(ticks),labels=ticks,mgp=c(3,0.5,0))
	for(i in 1:dim(f_rule)[1]){
		if(is.element(var,names(parameters_rule)[rule_relevant])){
			points(parameters_rule_to_plot[[var]][f_rule[i,]],squish(time[best_rule[f_rule[i,]]]),t='o',col=red_blue_purp[2],lty=i,lwd=lwd,pch=19)
		}else{
			points(unique(parameters_cat_to_plot[[var]]),rep(squish(time[best_rule[f_rule[i,1]]]),length(unique(parameters_cat[[var]]))),t='l',col=red_blue_purp[2],lty=i,lwd=lwd)
		}		
	}
	for(i in 1:dim(f_gen)[1]){
		if(is.element(var,names(parameters_genetic)[gen_relevant])){
			points(parameters_genetic[[var]][f_gen[i,]],best_gen[f_gen[i,]],t='o',col=red_blue_purp[3],lty=i,lwd=lwd,pch=19)
		}else{
			points(unique(parameters_cat_to_plot[[var]]),rep(best_gen[f_gen[i,1]],length(unique(parameters_cat[[var]]))),t='l',col=red_blue_purp[3],lty=i,lwd=lwd)
		}		
	}	
	if(var!='delta'){
		# for(i in 1:num_ex){
		i = 1
		points(parameters_cat_to_plot[[var]][f_cat[i,]],squish(time[best_cat[f_cat[i,]]]),t='o',col=red_blue_purp[1],lwd=lwd,lty=i,pch=19)
	# }
	}	
	for(i in 1:num_ex){
		points(parameters_cat_to_plot[[var]][f_cat2[i,]],squish(time[best_cat[f_cat2[i,]]]),t='o',col=purp_pal[2],lwd=lwd,lty=i,pch=19)
	}	
	if(var=='rho'){
	mtext(expression(paste('Stopping time, ',tau,sep='')),2,line=1.4,cex=largefontsize/smallfontsize)}
	plot(parameters_cat_to_plot[[var]],optimized_cat,t='n',xlab='',ylab='',xaxt='n',yaxt='n',ylim=c(0,0.35))
	if(var=='w'){
	axis(1,at=seq(500,2500,by=1000),labels=seq(500,2500,by=1000),mgp=c(3,0.5,0))
	axis(1,at=3000,labels='All',mgp=c(3,0.5,0))
	}else{
		axis(1,mgp=c(3,0.5,0),at=unique(parameters_cat[[var]]),labels=unique(parameters_cat[[var]]))}	
	axis(2,mgp=c(3,0.5,0),at=seq(0,0.4,by=0.1),labels=seq(0,0.4,by=0.1))
	for(i in 1:dim(f_rule)[1]){
		if(is.element(var,names(parameters_rule)[rule_relevant])){
			points(parameters_rule_to_plot[[var]][f_rule[i,]],optimized_rule[f_rule[i,]],t='o',col=red_blue_purp[2],lty=i,lwd=lwd,pch=19)
		}else{
			points(unique(parameters_cat_to_plot[[var]]),rep(optimized_rule[f_rule[i,1]],length(unique(parameters_cat[[var]]))),t='l',col=red_blue_purp[2],lty=i,lwd=lwd)
		}		
	}
	for(i in 1:dim(f_gen)[1]){
		if(is.element(var,names(parameters_genetic)[gen_relevant])){
			points(parameters_genetic[[var]][f_gen[i,]],optimized_gen[f_gen[i,]],t='o',col=red_blue_purp[3],lty=i,lwd=lwd,pch=19)
		}else{
			points(unique(parameters_cat_to_plot[[var]]),rep(optimized_gen[f_gen[i,1]],length(unique(parameters_cat[[var]]))),t='l',col=red_blue_purp[3],lty=i,lwd=lwd)
		}		
	}
	if(var!='delta'){
		# for(i in 1:num_ex){
		i = 1
		points(parameters_cat_to_plot[[var]][f_cat[i,]],optimized_cat[f_cat[i,]],t='o',col=red_blue_purp[1],lwd=lwd,lty=i,pch=19)
	# }
	}
	for(i in 1:num_ex){
		points(parameters_cat_to_plot[[var]][f_cat2[i,]],optimized_cat[f_cat2[i,]],t='o',col=purp_pal[2],lwd=lwd,lty=i,pch=19)
	}	
	xlim = range(unique(parameters_cat_to_plot[[var]]))
	mtext(labels[which(parameters_to_plot==var)],1,line=2,at=0.55*xlim[1]+0.45*xlim[2],cex=largefontsize/smallfontsize)
	ylim=range(optimized_cat,na.rm=TRUE)
	if(var=='rho'){
	mtext(expression(paste('Error, ',epsilon(tau),sep='')),2,line=1.4,cex=largefontsize/smallfontsize)}
		
}	

f_gen = find_df(examples[,gen_relevant],p=parameters_genetic)
w_to_plot = 5:11
plot(w_vec[w_to_plot],w_vec[w_to_plot],t='n',xlab='',ylab='',yaxt='n',ylim=time_ylim)		
axis(2,at=squish(ticks),labels=ticks,mgp=c(3,0.5,0))
for(i in 1:dim(f_rule)[1]){
		points(w_vec[w_to_plot],squish(time[best_c_rule[i,w_to_plot]]),t='o',col=red_blue_purp[2],lty=i,lwd=lwd,pch=19)			
}
for(i in 1:dim(f_gen)[1]){				points(w_vec[w_to_plot],rep(squish(time[best_gen[f_gen[i]]]),length(w_to_plot)),t='l',col=red_blue_purp[3],lty=i,lwd=lwd,pch=19)	
}	
	i = 1
	points(w_vec[w_to_plot],squish(time[best_c_recog[w_to_plot]]),t='o',col=red_blue_purp[1],lwd=lwd,lty=i,pch=19)	
examples_var[[var]] = NA
for(i in 1:num_ex){
	points(w_vec[w_to_plot],squish(time[best_c_cat[i,w_to_plot]]),t='o',col=purp_pal[2],lwd=lwd,lty=i,pch=19)
}
legend(0.35,y=8.7,legend=legend_labels[1:6],lty=c(rep(1,4),1,2),lwd=lwd,col=c(red_blue_purp[1],purp_pal[2],red_blue_purp[2:3],rep('black',num_ex)),bty='n',seg.len=2)
# legend(0.8,y=8.7,legend=legend_labels[5:6],lty=1:num_ex,lwd=lwd,col=rep('black',num_ex),bty='n',seg.len=2)
	

plot(w_vec[w_to_plot],w_vec[w_to_plot],t='n',xlab='',ylab='',xaxt='n',yaxt='n',ylim=c(0,0.35))	
axis(1,mgp=c(3,0.5,0))
axis(2,mgp=c(3,0.5,0),at=seq(0,0.4,by=0.1),labels=seq(0,0.4,by=0.1))
for(i in 1:dim(f_rule)[1]){
		points(w_vec[w_to_plot],opt_c_rule[i,w_to_plot],t='o',col=red_blue_purp[2],lty=i,lwd=lwd,pch=19)			
}
for(i in 1:dim(f_gen)[1]){				points(w_vec[w_to_plot],rep(optimized_gen[f_gen[i]],length(w_to_plot)),t='l',col=red_blue_purp[3],lty=i,lwd=lwd,pch=19)	
}	
	i = 1
	points(w_vec[w_to_plot],opt_c_recog[w_to_plot],t='o',col=red_blue_purp[1],lwd=lwd,lty=i,pch=19)	
examples_var[[var]] = NA
for(i in 1:num_ex){
	points(w_vec[w_to_plot],opt_c_cat[i,w_to_plot],t='o',col=purp_pal[2],lwd=lwd,lty=i,pch=19)
}	
xlim = range(w_vec[w_to_plot])
mtext(labels[4],1,line=2,at=0.55*xlim[1]+0.45*xlim[2],cex=largefontsize/smallfontsize)


dev.off()


######## parameter interaction heat maps


contour_error = list()
contour_time = list()

N = 25
rho = 0.99

f = find(N=N,rho=rho,l=0.1)
error_delta_w = data.frame(delta=parameters_cat$delta[f],w=parameters_cat$w[f],error=optimized_cat[f])
time_delta_w = data.frame(delta=parameters_cat$delta[f],w=parameters_cat$w[f],time=best_cat[f])
# error_delta_w=error_delta_w[-which(error_delta_w$w==250),]
error_delta_w$w[error_delta_w$w==250]=0
error_delta_w$w[which(error_delta_w$w==Inf)]=3000
# time_delta_w=time_delta_w[-which(time_delta_w$w==250),]
time_delta_w$w[time_delta_w$w==250]=0
time_delta_w$w[which(time_delta_w$w==Inf)]=3000
time_delta_w$time = squish(time[time_delta_w$time])

f = find(N=N,rho=rho,w=1000)
error_delta_l = data.frame(delta=parameters_cat$delta[f],l=parameters_cat$l[f],error=optimized_cat[f])
time_delta_l = data.frame(delta=parameters_cat$delta[f],l=parameters_cat$l[f],time=best_cat[f])
error_delta_l$l = replace(error_delta_l$l,1:20,c(0.05,0.125,0.275,0.5))
error_delta_l$width=replace(error_delta_l$l,1:20,values=c(0.05,0.1,0.2,0.25))
time_delta_l$l = replace(time_delta_l$l,1:20,c(0.05,0.125,0.275,0.5))
time_delta_l$width=replace(time_delta_l$l,1:20,values=c(0.05,0.1,0.2,0.25))
time_delta_l$time = squish(time[time_delta_l$time])

f = find(N=N,rho=rho,delta=0)
error_w_l = data.frame(w=parameters_cat$w[f],l=parameters_cat$l[f],error=optimized_cat[f])
time_w_l = data.frame(w=parameters_cat$w[f],l=parameters_cat$l[f],time=best_cat[f])
error_w_l$l = replace(error_w_l$l,1:28,c(0.05,0.125,0.275,0.5))
error_w_l$width=replace(error_w_l$l,1:28,values=c(0.05,0.1,0.2,0.25))
# error_w_l=error_delta_w[-which(error_delta_w$w==250),]
error_w_l$w[error_w_l$w==250]=0
error_w_l$w[which(error_w_l$w==Inf)]=3000
time_w_l$l = replace(time_w_l$l,1:28,c(0.05,0.125,0.275,0.5))
time_w_l$width=replace(time_w_l$l,1:28,values=c(0.05,0.1,0.2,0.25))
# time_w_l=time_delta_w[-which(time_w_l$w==250),]
time_w_l$w[time_w_l$w==250]=0
time_w_l$w[which(time_w_l$w==Inf)]=3000
time_w_l$time = squish(time[time_w_l$time])

M_error = range(c(range(error_delta_w$error,na.rm=TRUE),range(error_delta_l$error,na.rm=TRUE),range(error_w_l$error,na.rm=TRUE)))
M_time = range(c(range(time_delta_w$time,na.rm=TRUE),range(time_delta_l$time,na.rm=TRUE),range(time_w_l$time,na.rm=TRUE)))
error_cuts = seq(M_error[1],M_error[2],length.out=11)
time_cuts = seq(M_time[1],M_time[2],length.out=11)

error_fake = data.frame(z =seq(M_error[1],M_error[2],length.out = 22), x = 1:22,y = 1:22)
time_fake = data.frame(z =seq(M_time[1],M_time[2],length.out = 22), x = 1:22,y = 1:22)

marg = c(0.4,0,0.04,0.1)
big_marg = c(0.2,0.01,0.04,0.045)

contour_error[[1]] = ggplot(error_delta_w,aes(y=w,x=delta,z=error)) +
	geom_tile(aes(fill=error))+
	scale_fill_gradientn(colours=rev(red_pal),limits=M_error,guide='colorbar')+
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(big_marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25))+
	ylab('Memory window, w') + xlab(expression(paste('Category width, ',delta,sep='')))

contour_error[[2]] = ggplot(error_delta_l,aes(y=l,x=delta,z=error)) +
	geom_tile(aes(fill=error,height=width))+
	scale_fill_gradientn(colours=rev(red_pal),limits=M_error,guide='colorbar')+
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(big_marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_y_continuous(expand = c(0,0),breaks=unique(error_delta_l$l),labels=unique(parameters_cat$l)) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25))+
	ylab('Updating rate, r') + xlab(expression(paste('Category width, ',delta,sep='')))
	
contour_error[[3]] = ggplot(error_w_l,aes(y=w,x=l,z=error)) +
	geom_tile(aes(fill=error,width=width))+
	scale_fill_gradientn(colours=rev(red_pal),limits=M_error,guide='colorbar')+
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(big_marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_x_continuous(expand = c(0,0),breaks=unique(error_delta_l$l)+c(-0.01,0,0,0),labels=unique(parameters_cat$l)) + scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) +
	ylab('Memory window, w')+xlab('Updating rate, r') 

contour_time[[1]] = ggplot(time_delta_w,aes(y=w,x=delta,z=time)) +
	geom_tile(aes(fill=time))+
	scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,guide='colorbar')+
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(big_marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25))+
	ylab('Memory window, w') + xlab(expression(paste('Category width, ',delta,sep='')))

contour_time[[2]] = ggplot(time_delta_l,aes(y=l,x=delta,z=time)) +
	geom_tile(aes(fill=time,height=width))+
	scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,guide='colorbar',na.value=brewer.pal(9,'YlGnBu')[9])+
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(big_marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_y_continuous(expand = c(0,0),breaks=unique(time_delta_l$l),labels=unique(parameters_cat$l)) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25))+
	ylab('Updating rate, r') + xlab(expression(paste('Category width, ',delta,sep='')))
	
contour_time[[3]] = ggplot(time_w_l,aes(y=w,x=l,z=time)) +
	geom_tile(aes(fill=time,width=width))+
	scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,guide='colorbar')+
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(big_marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_x_continuous(expand = c(0,0),breaks=unique(error_delta_l$l)+c(-0.01,0,0,0),labels=unique(parameters_cat$l)) + scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) +
	ylab('Memory window, w')+xlab('Updating rate, r') 	
	
error_legend = ggplot(error_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(red_pal),limits=M_error,breaks = seq(round(M_error[1],1),round(M_error[2],1),by=0.1), labels = seq(round(M_error[1],1),round(M_error[2],1),by=0.1))+labs(fill=expression(paste('Error, ',epsilon,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

error_legend = get_legend(error_legend)

time_legend = ggplot(time_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,breaks = squish(c(1000,4000,16000)),labels=c(1000,4000,16000))+labs(fill=expression(paste('Time, ',tau,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), plot.title=element_text(size=smallfontsize),legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

time_legend = get_legend(time_legend)

graphics.off()

width = 6.85
height =4.5

pdf('/Users/eleanorbrush/Desktop/parameters_interactions_full.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(contour_time[[1]],contour_time[[2]],contour_time[[3]],time_legend,contour_error[[1]],contour_error[[2]],contour_error[[3]],error_legend,ncol=4,widths=c(0.4,0.4,0.4,0.16))	

dev.off()

###compare learned rule of thumb and recognition

marg = c(0.5,0.48,0.04,0.1)
omarg = c(0.25,0.33,0.3,0.0)

contour_diff = list()

M_diff = c(-0.37,0.37)

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		f_gen_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_gen)
			diff_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],c=optimal_cognition_cat$c[f_cat_now],diff=optimal_cognition_rule$cost[f_rule_now]-optimal_cognition_cat$cost[f_cat_now])		
		if(i==3 && j==2){
			diff_df$diff[abs(diff_df$diff)<0.04]=0
		}
		diff_df$contour = rep(NA,dim(diff_df)[1])
		hold = 2750
		for(q in 1:length(w_vec)){
			f = which(diff_df$c==w_vec[q])
			f2 = which(sign(diff_df$diff[f])!=sign(diff_df$diff[f])[1])
			f2 = c(f2[1],f2[length(f2)])
			f2 = diff_df$w[f2]
			if(!isempty(which(diff_df$diff[f]<=0))){
				f2= c(f2[1]-250,f2[2]+250)
				diff_df$contour[f] = c(rep(min(hold,f2[2]),length(f)-1),f2[1])
				hold = f2[1]
			}
		}
		contour_diff[[I]] = ggplot(diff_df,aes(x=c,y=w,z=diff)) +
			geom_tile(aes(fill=diff))+
			scale_fill_gradientn(colours=rev(div_pal),limits=M_diff,guide='colorbar')+			
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
			if(!isempty(which(!is.na(diff_df$contour)))){
				if(length(which(!is.na(diff_df$contour)))>1){
				# contour_diff[[I]] = contour_diff[[I]] +geom_line(aes(x=c-0.05,y=contour),colour='black')
				}
			}
	}
}

diff_fake = data.frame(z =seq(M_diff[1],M_diff[2],length.out = 22), x = 1:22,y = 1:22)

diff_legend = ggplot(diff_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(div_pal),limits=M_diff,breaks = seq(round(M_diff[1],1),round(M_diff[2],1),by=0.1), labels = seq(round(M_diff[1],1),round(M_diff[2],1),by=0.1))+labs(fill='Difference')+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

diff_legend = get_legend(diff_legend)

graphics.off()

width = 6.85
height= 4 

blank <- grid.rect(gp=gpar(col="white"))

pdf('/Users/eleanorbrush/Desktop/recog_vs_learned_rule.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_diff[[4]],contour_diff[[5]],contour_diff[[6]],diff_legend,blank,contour_diff[[7]],contour_diff[[8]],contour_diff[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.15))

dev.off()

####best assessment strategy

contour_diff = list()
hold = list()

M_diff = c(1,5)

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		f_gen_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_gen)
			diff_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],c=optimal_cognition_cat$c[f_cat_now],diff=optimal_cognition_rule$cost[f_rule_now]-optimal_cognition_cat$cost[f_cat_now])	
		if(i==3 && j==2){
			flip = 0.05
		}else{flip = 0.05}
		for(k in 1:length(f_cat_now)){
			if(abs(optimal_cognition_rule$cost[f_rule_now[k]]-optimal_cognition_gen$cost[f_gen_now[k]])<flip){
				if(abs(mean(optimal_cognition_rule$cost[f_rule_now[k]],optimal_cognition_gen$cost[f_gen_now[k]])-optimal_cognition_cat$cost[f_cat_now[[k]]])<flip){
					diff_df$diff[k] = 5
				}else if(mean(optimal_cognition_rule$cost[f_rule_now[k]],optimal_cognition_gen$cost[f_gen_now[k]])<optimal_cognition_cat$cost[f_cat_now[[k]]]){
						diff_df$diff[k] = 4
					}else{
						diff_df$diff[k] = 1
					}
				}else{
					diff_df$diff[k] = which.min(c(optimal_cognition_cat$cost[f_cat_now[[k]]],optimal_cognition_rule$cost[f_rule_now[[k]]],optimal_cognition_gen$cost[f_gen_now[[k]]]))
				}
			}							
		contour_diff[[I]] = ggplot(diff_df,aes(x=c,y=w,z=diff)) +
			geom_tile(aes(fill=diff))+
			scale_fill_gradientn(colours=(red_blue_purp),limits=M_diff,guide='colorbar')+		
			# scale_fill_manual(values=hold[[I]])+			
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.1),limits=c(0.55,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')			
	}
}

diff_fake = data.frame(z =cut(c(1:6),breaks=0:6), x = 1:6,y = 1:6)

diff_legend = ggplot(diff_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_manual(values=(c(red_blue_purp[1],purp_pal[2],red_blue_purp[2:5])),labels =c('Recog','Categ','Learned rule','Innate rule','Either rule','Recog / rule'))+labs(fill='')+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

diff_legend = get_legend(diff_legend)

graphics.off()

width = 6.85
height= 5.5

blank <- grid.rect(gp=gpar(col="white"))

pdf('/Users/eleanorbrush/Desktop/best_type_of_learning.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_diff[[1]],contour_diff[[2]],contour_diff[[3]],diff_legend,blank,contour_diff[[4]],contour_diff[[5]],contour_diff[[6]],blank,blank,contour_diff[[7]],contour_diff[[8]],contour_diff[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.23))

dev.off()

####optimized strategies

marg = c(0.5,0.48,0.04,0.1)
omarg = c(0.25,0.33,0.3,0.0)

M_delta = c(0,0.25)

contour_delta = list()

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		delta_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],c=optimal_cognition_cat$c[f_cat_now],delta=optimal_cognition_cat$delta[f_cat_now])		
		contour_delta[[I]] = ggplot(delta_df,aes(x=c,y=w,z=delta)) +
			geom_tile(aes(fill=delta))+
			scale_fill_gradientn(colours=c(red_blue_purp[1],purp_pal[2]),limits=M_delta,guide='colorbar')+
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
	}
}

delta_fake = data.frame(z =cut(M_delta,breaks = 2), x = 1:2,y = 1:2)

# delta_legend = ggplot(delta_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(red_pal),limits=M_delta,breaks = seq(0,M_delta[2],by=0.25), labels = seq(0,M_delta[2],by=0.25))+labs(fill='')+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

delta_legend = ggplot(delta_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=factor(z))) + scale_fill_manual(values=c(purp_pal[2],red_blue_purp[1]),labels=c(0.25,0))+labs(fill=expression(paste('Cat. width, ',delta,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

delta_legend = get_legend(delta_legend)

graphics.off()

M_l = c(0.05,0.5)

contour_l = list()

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		l_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],c=optimal_cognition_cat$c[f_cat_now],l=cut(optimal_cognition_cat$l[f_cat_now],breaks=c(0.025,0.075,0.15,0.375,0.55)))		
		contour_l[[I]] = ggplot(l_df,aes(x=c,y=w,z=l)) +
			geom_tile(aes(fill=l))+
			scale_fill_manual(values=green_pal,limits=levels(l_df$l))+
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
	}
}

l_fake = data.frame(z =cut(seq(M_l[1],M_l[2],length.out=4),breaks = 4), x = 1:4,y = 1:4)

l_legend = ggplot(l_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=factor(z))) + scale_fill_manual(values=rev(green_pal),labels=rev(unique(parameters_cat$l)))+labs(fill='Rate, r')+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

l_legend = get_legend(l_legend)

graphics.off()

M_error = c(0,0.30)

contour_error = list()

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		error_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],c=optimal_cognition_cat$c[f_cat_now],error=optimal_cognition_cat$error[f_cat_now])		
		contour_error[[I]] = ggplot(error_df,aes(x=c,y=w,z=error)) +
			geom_tile(aes(fill=error))+
			scale_fill_gradientn(colours=rev(red_pal),limits=M_error,guide='colorbar')+
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
	}
}

error_fake = data.frame(z =seq(M_error[1],M_error[2],length.out = 22), x = 1:22,y = 1:22)

error_legend = ggplot(error_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(red_pal),limits=M_error,breaks = seq(round(M_error[1],1),round(M_error[2],1),by=0.1), labels = seq(round(M_error[1],1),round(M_error[2],1),by=0.1))+labs(fill=expression(paste('Error, ',epsilon,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

error_legend = get_legend(error_legend)

graphics.off()

M_time = squish(time[c(10,601)])

contour_time = list()

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		time_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],c=optimal_cognition_cat$c[f_cat_now],time=squish(time[optimal_cognition_cat$time[f_cat_now]]))		
		contour_time[[I]] = ggplot(time_df,aes(x=c,y=w,z=time)) +
			geom_tile(aes(fill=time))+
			scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,guide='colorbar')+
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
	}
}

time_fake = data.frame(z =seq(M_time[1],M_time[2],length.out = 22), x = 1:22,y = 1:22)

time_legend = ggplot(time_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,breaks = squish(c(1000,4000,16000)),labels=c(1000,4000,16000))+labs(fill=expression(paste('Time, ',tau,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"),legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

time_legend = get_legend(time_legend)

graphics.off()

blank <- grid.rect(gp=gpar(col="white"))

width = 6.85
height= 6

pdf('/Users/eleanorbrush/Desktop/delta_heat_maps.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_delta[[1]],contour_delta[[2]],contour_delta[[3]],delta_legend,blank,contour_delta[[4]],contour_delta[[5]],contour_delta[[6]],blank,blank,contour_delta[[7]],contour_delta[[8]],contour_delta[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.17))

dev.off() 

width = 6.85
height= 6

pdf('/Users/eleanorbrush/Desktop/l_heat_maps.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_l[[1]],contour_l[[2]],contour_l[[3]],l_legend,blank,contour_l[[4]],contour_l[[5]],contour_l[[6]],blank,blank,contour_l[[7]],contour_l[[8]],contour_l[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.17))

dev.off()

width = 6.85
height= 6

pdf('/Users/eleanorbrush/Desktop/time_heat_maps.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_time[[1]],contour_time[[2]],contour_time[[3]],time_legend,blank,contour_time[[4]],contour_time[[5]],contour_time[[6]],blank,blank,contour_time[[7]],contour_time[[8]],contour_time[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.17))

dev.off()

width = 6.85
height= 6

pdf('/Users/eleanorbrush/Desktop/error_heat_maps.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_error[[1]],contour_error[[2]],contour_error[[3]],error_legend,blank,contour_error[[4]],contour_error[[5]],contour_error[[6]],blank,blank,contour_error[[7]],contour_error[[8]],contour_error[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.19))

dev.off()

width = 6.85
height= 4

pdf('/Users/eleanorbrush/Desktop/strategies_heat_maps.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_delta[[4]],contour_delta[[5]],contour_delta[[6]],delta_legend,blank,contour_delta[[7]],contour_delta[[8]],contour_delta[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.17))

dev.off()

#### learned rule heat maps
M_time = squish(time[c(10,601)])

contour_time = list()

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		time_df = data.frame(w=optimal_cognition_rule$w[f_rule_now],c=optimal_cognition_rule$c[f_rule_now],time=squish(time[optimal_cognition_rule$time[f_rule_now]]))		
		contour_time[[I]] = ggplot(time_df,aes(x=c,y=w,z=time)) +
			geom_tile(aes(fill=time))+
			scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,guide='colorbar')+
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
	}
}

time_fake = data.frame(z =seq(M_time[1],M_time[2],length.out = 22), x = 1:22,y = 1:22)

time_legend = ggplot(time_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(blue_pal),limits=M_time,breaks = squish(c(1000,4000,16000)),labels=c(1000,4000,16000))+labs(fill=expression(paste('Time, ',tau,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

time_legend = get_legend(time_legend)

graphics.off()


width = 6.85
height= 4

pdf('/Users/eleanorbrush/Desktop/time_heat_maps_rule.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_time[[4]],contour_time[[5]],contour_time[[6]],time_legend,blank,contour_time[[7]],contour_time[[8]],contour_time[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.17))

dev.off()

M_error = c(0,0.30)

contour_error = list()

for(i in 1:3){
	for(j in 1:3){
		I = sub2ind(c(xN_cat,xrho_cat),c(i,j))
		f_cat_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
		f_rule_now = find(N=unique(examples_optimized$N)[i],rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
		error_df = data.frame(w=optimal_cognition_rule$w[f_rule_now],c=optimal_cognition_rule$c[f_rule_now],error=optimal_cognition_rule$error[f_rule_now])		
		contour_error[[I]] = ggplot(error_df,aes(x=c,y=w,z=error)) +
			geom_tile(aes(fill=error))+
			scale_fill_gradientn(colours=rev(red_pal),limits=M_error,guide='colorbar')+
			theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
			scale_y_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_x_continuous(expand = c(0,0),breaks=seq(0,1,by=0.25),limits=c(0.45,1.05))+
			ylab('Memory window, w') + xlab('Cost of errors, c')
	}
}

error_fake = data.frame(z =seq(M_error[1],M_error[2],length.out = 22), x = 1:22,y = 1:22)

error_legend = ggplot(error_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=rev(red_pal),limits=M_error,breaks = seq(round(M_error[1],1),round(M_error[2],1),by=0.1), labels = seq(round(M_error[1],1),round(M_error[2],1),by=0.1))+labs(fill=expression(paste('Error, ',epsilon,sep='')))+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

error_legend = get_legend(error_legend)

graphics.off()

width = 6.85
height= 4

pdf('/Users/eleanorbrush/Desktop/error_heat_maps_rule.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(blank,contour_error[[4]],contour_error[[5]],contour_error[[6]],error_legend,blank,contour_error[[7]],contour_error[[8]],contour_error[[9]],ncol=5,widths=c(0.05,0.4,0.4,0.4,0.17))

dev.off()

####

width = 6.85
height= 3

marg = c(0.5,0.5,0.2,0.1)
omarg = c(0.1,0.33,0.4,0.1)

pdf('/Users/eleanorbrush/Desktop/stopping_time.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)

layout(matrix(1:3,nrow=1))

w_vals = c(1,0.9,0.6)

f = find(N=25,rho=0.9,w=2500,delta=0,l=0.05)

error = error_cat[,f]

for(i in 1:3){
	w = w_vals[i]
	opt_now = optimization_cat(error,w)
	cost = 2*w*error+(1-w)*time/30001
	plot(time,cost,ylim=c(0,0.5),t='l',lwd=lwd,xlab='',ylab='',cex=largefontsize/smallfontsize,xaxt='n')
	axis(1,seq(0,30000,by=10000))
	mtext(paste('c = ',w,sep=''),cex=largefontsize/smallfontsize,side=3,line=0.5)
	mtext('Time, t',side=1,line=2.3,cex=largefontsize/smallfontsize)
	mtext('Total cost, C(t)',side=2,line=2.2,cex=largefontsize/smallfontsize)
	points(time[opt_now$best_cat],opt_now$optimized_cat,col=red_blue_purp[1],lwd=lwd,pch=19)
}

dev.off()