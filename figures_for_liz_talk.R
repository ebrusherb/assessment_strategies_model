library(lattice)
library(pracma)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
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
load('opt_data.Rdata')

best_cat = opt_no_time_costs$best_cat
optimized_cat = opt_no_time_costs$optimized_cat/2
best_rule = opt_no_time_costs$best_rule
optimized_rule = opt_no_time_costs$optimized_rule/2
best_gen = rep(squish(10),dim(parameters_genetic)[1])
optimized_gen = error_genetic[1,]


### effect of parameters


marg = c(0.3,0.3,0.04,0.05)
omarg = c(0.7,0.65,0.65,0.05)

width = 5
height =2

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

pdf('/Users/eleanorbrush/Desktop/parameters_exploration_full.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)

layout(matrix(1:3,nrow=1))

for(var in parameters_to_plot){
	examples_var = examples
	examples_var[[var]] = NA
	f_cat = find_df(examples_var)
	f_rule = find_df(examples_var[,rule_relevant],p=parameters_rule)
	f_gen = find_df(examples_var[,gen_relevant],p=parameters_genetic)
	examples_var = examples_delta
	examples_var[[var]] = NA
	f_cat2 = find_df(examples_var)	
		plot(parameters_cat_to_plot[[var]],optimized_cat,t='n',xlab='',ylab='',xaxt='n',yaxt='n',ylim=c(0,0.35))
	if(var=='w'){
	axis(1,at=seq(500,2500,by=1000),labels=seq(500,2500,by=1000),mgp=c(3,0.5,0))
	axis(1,at=3000,labels='All',mgp=c(3,0.5,0))
	}else{
		axis(1,mgp=c(3,0.5,0),at=unique(parameters_cat[[var]]),labels=unique(parameters_cat[[var]]))}	
	axis(2,mgp=c(3,0.5,0),at=seq(0,0.4,by=0.1),labels=seq(0,0.4,by=0.1))
	for(i in 1:1){
		if(is.element(var,names(parameters_rule)[rule_relevant])){
			points(parameters_rule_to_plot[[var]][f_rule[i,]],optimized_rule[f_rule[i,]],t='o',col=red_blue_purp[2],lty=i,lwd=lwd,pch=19)
		}else{
			points(unique(parameters_cat_to_plot[[var]]),rep(optimized_rule[f_rule[i,1]],length(unique(parameters_cat[[var]]))),t='l',col=red_blue_purp[2],lty=i,lwd=lwd)
		}		
	}
	for(i in 1:1){
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
	for(i in 1:1){
		points(parameters_cat_to_plot[[var]][f_cat2[i,]],optimized_cat[f_cat2[i,]],t='o',col=purp_pal[2],lwd=lwd,lty=i,pch=19)
	}	
	xlim = range(unique(parameters_cat_to_plot[[var]]))
	mtext(labels[which(parameters_to_plot==var)],1,line=2,at=0.55*xlim[1]+0.45*xlim[2],cex=largefontsize/smallfontsize)
	ylim=range(optimized_cat,na.rm=TRUE)
	if(var=='rho'){
	mtext(expression(paste('Error, ',epsilon(tau),sep='')),2,line=1.4,cex=largefontsize/smallfontsize)}
		
}	

legend(530,y=0.37,legend=legend_labels[1:4],lty=c(rep(1,4)),lwd=lwd,col=c(red_blue_purp[1],purp_pal[2],red_blue_purp[2:3]),bty='n',seg.len=2)

dev.off()


pdf('/Users/eleanorbrush/Desktop/parameters_exploration_just_categorization.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)

layout(matrix(1:3,nrow=1))

for(var in parameters_to_plot){
	examples_var = examples
	examples_var[[var]] = NA
	f_cat = find_df(examples_var)
	f_rule = find_df(examples_var[,rule_relevant],p=parameters_rule)
	f_gen = find_df(examples_var[,gen_relevant],p=parameters_genetic)
	examples_var = examples_delta
	examples_var[[var]] = NA
	f_cat2 = find_df(examples_var)	
		plot(parameters_cat_to_plot[[var]],optimized_cat,t='n',xlab='',ylab='',xaxt='n',yaxt='n',ylim=c(0,0.35))
	if(var=='w'){
	axis(1,at=seq(500,2500,by=1000),labels=seq(500,2500,by=1000),mgp=c(3,0.5,0))
	axis(1,at=3000,labels='All',mgp=c(3,0.5,0))
	}else{
		axis(1,mgp=c(3,0.5,0),at=unique(parameters_cat[[var]]),labels=unique(parameters_cat[[var]]))}	
	axis(2,mgp=c(3,0.5,0),at=seq(0,0.4,by=0.1),labels=seq(0,0.4,by=0.1))
	
	if(var!='delta'){
		# for(i in 1:num_ex){
		i = 1
		points(parameters_cat_to_plot[[var]][f_cat[i,]],optimized_cat[f_cat[i,]],t='o',col=red_blue_purp[1],lwd=lwd,lty=i,pch=19)
	# }
	}
	for(i in 1:1){
		points(parameters_cat_to_plot[[var]][f_cat2[i,]],optimized_cat[f_cat2[i,]],t='o',col=purp_pal[2],lwd=lwd,lty=i,pch=19)
	}	
	xlim = range(unique(parameters_cat_to_plot[[var]]))
	mtext(labels[which(parameters_to_plot==var)],1,line=2,at=0.55*xlim[1]+0.45*xlim[2],cex=largefontsize/smallfontsize)
	ylim=range(optimized_cat,na.rm=TRUE)
	if(var=='rho'){
	mtext(expression(paste('Error, ',epsilon(tau),sep='')),2,line=1.4,cex=largefontsize/smallfontsize)}
		
}	

legend(530,y=0.37,legend=legend_labels[1:2],lty=c(rep(1,2)),lwd=lwd,col=c(red_blue_purp[1],purp_pal[2]),bty='n',seg.len=2)

dev.off()




###compare categorization and recognition

marg = c(0.5,0.48,0.04,0.1)
omarg = c(0.25,0.33,0.3,0.0)

rho = 0.99
l = 0.5
		
f_recog_now = find(rho=rho,l=l,delta=0)
f_cat_now = find(rho=rho,l=l,delta=0.25)
diff_df = data.frame(w=parameters_cat$w[f_cat_now],N=parameters_cat$N[f_cat_now],diff=opt_no_time_costs$optimized_cat[f_cat_now]-opt_no_time_costs$optimized_cat[f_recog_now])		
diff_df$w[diff_df$w==Inf]=3000
diff_df$w[diff_df$w==250]=0
diff_df$height = diff_df$N
diff_df$height = replace(diff_df$height,1:21,c(15,35,70))
M_diff = c(-max(abs(diff_df$diff)),max(abs(diff_df$diff)))
contour_diff = ggplot(diff_df,aes(x=w,y=N,z=diff,height=height)) +
	geom_tile(aes(fill=diff))+
	scale_fill_gradientn(colours=c(rev(brewer.pal(5,'Purples')),'white',brewer.pal(5,'Reds')),limits=M_diff,guide='colorbar')+			
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_x_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_y_continuous(expand = c(0,0),breaks=c(25,50,100),labels=c(25,50,100))+
	xlab('Memory window, w') + ylab('Group size, N')

diff_fake = data.frame(z =seq(M_diff[1],M_diff[2],length.out = 22), x = 1:22,y = 1:22)

diff_legend = ggplot(diff_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_gradientn(colours=c(rev(brewer.pal(5,'Purples')),'white',brewer.pal(5,'Reds')),limits=M_diff,breaks = seq(round(M_diff[1],1),round(M_diff[2],1),by=0.1), labels = seq(round(M_diff[1],1),round(M_diff[2],1),by=0.1))+labs(fill='Difference')+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

diff_legend = get_legend(diff_legend)

graphics.off()

width = 3.5
height= 2.7

blank <- grid.rect(gp=gpar(col="white"))

pdf('/Users/eleanorbrush/Desktop/recog_vs_categorization.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(contour_diff,diff_legend,ncol=2,widths=c(1,0.3))

dev.off()
graphics.off()

####best assessment strategy

c = 1
M_diff = c(1,5)

f_cat_now = find(c=c,rho=unique(examples_optimized$rho)[j],p=optimal_cognition_cat)
f_cat_now = f_cat_now[order(optimal_cognition_cat[f_cat_now,]$w)]
f_rule_now = find(c=c,rho=unique(examples_optimized$rho)[j],p=optimal_cognition_rule)
f_rule_now = f_rule_now[order(optimal_cognition_rule[f_rule_now,]$w)]
f_gen_now = find(c=c,rho=unique(examples_optimized$rho)[j],p=optimal_cognition_gen)
f_gen_now = f_gen_now[order(optimal_cognition_gen[f_gen_now,]$w)]
	best_df = data.frame(w=optimal_cognition_cat$w[f_cat_now],N=optimal_cognition_cat$N[f_cat_now],diff=optimal_cognition_rule$cost[f_rule_now]-optimal_cognition_cat$cost[f_cat_now])	
flip = 0.05
for(k in 1:length(f_cat_now)){
	if(abs(optimal_cognition_rule$cost[f_rule_now[k]]-optimal_cognition_gen$cost[f_gen_now[k]])<flip){
		if(abs(mean(optimal_cognition_rule$cost[f_rule_now[k]],optimal_cognition_gen$cost[f_gen_now[k]])-optimal_cognition_cat$cost[f_cat_now[[k]]])<flip){
			best_df$diff[k] = 5
		}else if(mean(optimal_cognition_rule$cost[f_rule_now[k]],optimal_cognition_gen$cost[f_gen_now[k]])<optimal_cognition_cat$cost[f_cat_now[[k]]]){
				best_df$diff[k] = 4
			}else{
				best_df$diff[k] = 1
			}
		}else{
			best_df$diff[k] = which.min(c(optimal_cognition_cat$cost[f_cat_now[[k]]],optimal_cognition_rule$cost[f_rule_now[[k]]],optimal_cognition_gen$cost[f_gen_now[[k]]]))
		}
	}	
best_df$height = rep(c(15,35,70),times=7)
contour_diff = ggplot(best_df,aes(x=w,y=N,z=diff,height=height)) +
	geom_tile(aes(fill=diff))+
	scale_fill_gradientn(colours=(red_blue_purp),limits=M_diff,guide='colorbar')+		
	theme_bw()+ theme(text=element_text(family="Helvetica", size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none')+
	scale_x_continuous(expand = c(0,0),breaks=seq(000,3000,by=500),labels=c(250,seq(500,2500,by=500),'All')) + scale_y_continuous(expand = c(0,0),breaks=c(25,50,100),labels=c(25,50,100))+
	xlab('Memory window, w') + ylab('Group size, N')	


diff_fake = data.frame(z =cut(c(M_diff[1]:M_diff[2]),breaks=0:5), x = 1:5,y = 1:5)

diff_legend = ggplot(diff_fake,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_manual(values=(red_blue_purp),labels =c('Recog','Learned','Innate','Either','Any'))+labs(fill='')+theme(text=element_text(family="Helvetica", size=smallfontsize), plot.title=element_text(size=smallfontsize), plot.margin=unit(marg,"cm"), legend.key =element_blank(),legend.text=element_text(size=smallfontsize))

diff_legend = get_legend(diff_legend)

graphics.off()

width = 3.5
height= 2.7

blank <- grid.rect(gp=gpar(col="white"))

pdf('/Users/eleanorbrush/Desktop/best_type_of_learning.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize)

grid.arrange(contour_diff,diff_legend,ncol=2,widths=c(1,0.3))

dev.off()
graphics.off()

