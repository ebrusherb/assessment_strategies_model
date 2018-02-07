library(lattice)
library(pracma)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(RColorBrewer)
library(reshape)
setwd('/Users/eleanorbrush/Dropbox/evo_badgesVSrecognition/code_extended_model')
source('get_legend.R')
load('noise=0.1/all_parameters.Rdata')
source('sub2ind.R')
source('ind2sub.R')
source('model.R')

#### palettes
set1cols = brewer.pal(9,'Set1')
set1cols_xN = brewer.pal(3,'Set1')

#####  how categorization works
N=100
qual_vals = rnorm(N, mean = 0, sd = 1) 
sig_vals = rnorm(N, mean = 0, 1)
sig_qual_corr = 0.5
sig_vals = fixCorr(qual_vals,sig_vals,sig_qual_corr)
confus_prob_cat = Inf

perc_wind = 1

sig_cats_temp = array(NA, dim = N)
left = 1:N
cat_now = 1
cat_med = array(NA,dim=0)
while(length(left)>0){
	if(length(left)>1){
	first = sample(left,size=1)} else{ first = left}
	to_categorize = which(abs(sig_vals[left]-sig_vals[first])<=perc_wind/2)
	if(length(which(exp(-confus_prob_cat*abs(sig_vals[left[to_categorize]]-median(sig_vals[left[to_categorize]])))==0))>0 && confus_prob_cat!=Inf ){
				to_categorize = to_categorize[-which(exp(-confus_prob_cat*abs(sig_vals[left[to_categorize]]-median(sig_vals[left[to_categorize]])))==0)]
			}
	sig_cats_temp[left[to_categorize]] = cat_now
	cat_now = cat_now+1
	cat_med = c(cat_med,mean(sig_vals[left[to_categorize]]))
	left = left[-to_categorize]
}


o = order(cat_med)
sig_cats = array(NA, dim =N)
for(i in 1:length(o)){
	sig_cats[sig_cats_temp==o[i]] = i
}
cat_med = sort(cat_med)
catnum = max(sig_cats)

categories = data.frame(sig_vals=c(sig_vals,sig_vals),sig_cats=c(sig_vals,cat_med[sig_cats]),delta=c(rep(-1,N),rep(0,N)))

perc_wind = 2

sig_cats_temp = array(NA, dim = N)
left = 1:N
cat_now = 1
cat_med = array(NA,dim=0)
while(length(left)>0){
	if(length(left)>1){
	first = sample(left,size=1)} else{ first = left}
	to_categorize = which(abs(sig_vals[left]-sig_vals[first])<=perc_wind/2)
	if(length(which(exp(-confus_prob_cat*abs(sig_vals[left[to_categorize]]-median(sig_vals[left[to_categorize]])))==0))>0 && confus_prob_cat!=Inf ){
				to_categorize = to_categorize[-which(exp(-confus_prob_cat*abs(sig_vals[left[to_categorize]]-median(sig_vals[left[to_categorize]])))==0)]
			}
	sig_cats_temp[left[to_categorize]] = cat_now
	cat_now = cat_now+1
	cat_med = c(cat_med,mean(sig_vals[left[to_categorize]]))
	left = left[-to_categorize]
}


o = order(cat_med)
sig_cats = array(NA, dim =N)
for(i in 1:length(o)){
	sig_cats[sig_cats_temp==o[i]] = i
}
cat_med = sort(cat_med)
catnum = max(sig_cats)

categories = rbind(categories,data.frame(sig_vals=sig_vals,sig_cats=cat_med[sig_cats],delta=c(rep(1,N))))

categories$colors =(round(1/2*(categories$sig_cats+1),2)*100)
num_cols=diff(range(categories$colors))+1
offset = 5
warm_pal = rev(colorRampPalette(brewer.pal(9,'Spectral'))(num_cols+offset))
warm_pal = warm_pal[(1:num_cols)+offset]
warm_pal = warm_pal[sort(unique(categories$colors-min(categories$colors)+1))]
categories$colors=as.factor(categories$colors)

m=min(-1,min(sig_vals))
M=max(1,max(sig_vals))

plot_cats = ggplot(categories, aes(x = delta, y = sig_vals, colour = colors)) + 
		geom_point(aes(size=sig_cats))  +
		theme_bw() +
		theme(
        	 axis.line.y=element_blank(),text=element_text(family="Helvetica", size=10), plot.title=element_text(size=10) ,plot.title=element_text(size=10),legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin=unit(c(1,1,1,1),"cm"),panel.border=element_blank()) +		scale_color_manual(values=warm_pal)+	scale_y_continuous(limit=c(m,M))+scale_x_continuous(limit=c(-1.1,1.5),breaks=c(-1,0,1),labels=c(expression(paste(delta,' = 0',sep='')),expression(paste(delta,' = 1',sep='')),expression(paste(delta,' = 2',sep=''))))+
		xlab("Category width")+ylab('Signal')
				
# pdf(file=paste(wd,'/categories.pdf',sep=''),width=6,height=4)
# print(plot_cats)
# dev.off()

#### how many categories?

N_vals = c(25,50,100)
xN = length(N_vals)
perc_vals = seq(2,0,by=-0.25)
xp = length(perc_vals)
corr_vals = c(0.9)
xc = length(corr_vals)

ybreaks = c(1,2,3,6,12,25,50,100)

cat_mean = array(NA,dim=c(xN,xp,xc))
cat_sd = array(NA,dim=c(xN,xp,xc))

runs = 5000

# # for(n in 1:xN){
	# for(p in 1:xp){
		# for(c2 in 1:xc){
			# N = N_vals[n]
			# perc_wind = perc_vals[p]
			# sig_qual_corr = corr_vals[c2]
			# hold_catnum = array(NA,dim=runs)
			# for(r in 1:runs){
				# qual_vals = rnorm(N, mean = 0, sd = 1) 
				# sig_vals = rnorm(N, mean = 0, 1)
				# sig_vals = fixCorr(qual_vals,sig_vals,sig_qual_corr)
				# sig_cats_temp = array(NA, dim = N)
				# left = 1:N
				# cat_now = 1
				# cat_med = array(NA,dim=0)
				# while(length(left)>0){
					# if(length(left)>1){
					# first = sample(left,size=1)} else{ first = left}
					# to_categorize = which(abs(sig_vals[left]-sig_vals[first])<=perc_wind/2)
					# if(length(which(exp(-confus_prob_cat*abs(sig_vals[left[to_categorize]]-median(sig_vals[left[to_categorize]])))==0))>0 && confus_prob_cat!=Inf ){
								# to_categorize = to_categorize[-which(exp(-confus_prob_cat*abs(sig_vals[left[to_categorize]]-median(sig_vals[left[to_categorize]])))==0)]
							# }
					# sig_cats_temp[left[to_categorize]] = cat_now
					# cat_now = cat_now+1
					# cat_med = c(cat_med,mean(sig_vals[left[to_categorize]]))
					# left = left[-to_categorize]
				# }
				# hold_catnum[r] = max(sig_cats_temp)
			# }
		# cat_mean[n,p,c2] = mean(hold_catnum)
		# cat_sd[n,p,c2] = sd(hold_catnum)	
		# }
	# }
# }

# cat_summary = melt(cat_mean,varnames=c('groupsize','percwind','corr'))
# names(cat_summary)[4]='catnum'
# cat_summary$groupsize = as.factor(rep(N_vals,xp*xc))
# cat_summary$percwind = rep(perc_vals,each=xN,times=xc)
# cat_summary$corr = as.factor(rep(corr_vals,each=xN*xp))
# cat_summary$sd = melt(cat_sd)$value
# cat_summary$min = cat_summary$catnum-cat_summary$sd
# cat_summary$max = cat_summary$catnum+cat_summary$sd
# cat_summary$catnum = log(cat_summary$catnum)
# cat_summary$min = log(cat_summary$min)
# cat_summary$max = log(cat_summary$max)

plot_catnum = ggplot(cat_summary, aes(x=percwind,y=catnum, colour = groupsize)) + 	
	geom_line() + geom_point() + 
	theme_bw() + theme(panel.border=element_blank())+scale_y_continuous(limits=log(c(1.5,101)),breaks=log(ybreaks),labels=ybreaks) +
		theme(text=element_text(family="Helvetica", size=10), plot.title=element_text(size=10) ,legend.key = element_blank(),legend.position = c(.9, .95),
  legend.justification = c("right", "top")) + 
		scale_color_manual(values=set1cols_xN[1:xN]) + 
		labs( colour="Group size")  + xlab(expression(paste('Category width, ',delta,sep='')))+ylab("Number of categories")

		
group_legend = get_legend(plot_catnum)

plot_catnum = 	plot_catnum + geom_ribbon(aes(ymin=min,ymax=max,x=percwind,fill=groupsize,colour=NA),alpha=0.2) + theme(legend.position='none')

# pdf(file=paste(wd,'/number_of_categories.pdf',sep=''),width=5,height=3.4)
# print(plot_catnum)
# dev.off()

pdf(file='/Users/eleanorbrush/Desktop/category_diagram.pdf',width=6.85,height=3.4)
grid.arrange(plot_cats,plot_catnum,group_legend,ncol=3,widths=c(1.2,1,0.1))
dev.off()

