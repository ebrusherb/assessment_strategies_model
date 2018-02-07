direc = '/Users/eleanorbrush/Dropbox/evo_badgesVSrecognition/code_extended_model/noise=0.1'
setwd(direc)

files = list.files(path=direc)

paramsweep = intersect(grep('summary_stats',files,value=FALSE),grep('Rdata',files,value=FALSE))

ord=order(as.numeric(gsub('.Rdata','',gsub('summary_stats_','',files[paramsweep]))))
paramsweep = paramsweep[ord]

error_cat = c()
error_rule = c()
error_genetic = c()
parameters_cat = c()
parameters_rule = c()
parameters_genetic = c()
hold_env = new.env()

for(i in 1:length(paramsweep)){
	to_load = files[paramsweep[i]]
	load(to_load,hold_env)
	error_cat = cbind(error_cat,hold_env$error_cat_mean)
	parameters_cat=rbind(parameters_cat,hold_env$parameters)
	if(!is.na(hold_env$error_rule_mean[1]) || is.nan(hold_env$error_rule_mean[1])){
		error_rule = cbind(error_rule,hold_env$error_rule_mean)
		parameters_rule=rbind(parameters_rule,hold_env$parameters[1:length(unique(hold_env$parameters$N)),])
	}		
	if(!is.na(hold_env$error_genetic_mean[1])){
		error_genetic = cbind(error_genetic,hold_env$error_genetic_mean)
		parameters_genetic=rbind(parameters_genetic,hold_env$parameters[1:length(unique(hold_env$parameters$N)),])
	}	
}

fixed_parameters = hold_env$fixed_parameters

save(error_cat,error_rule,error_genetic,parameters_cat,parameters_rule,parameters_genetic,fixed_parameters,file=paste(direc,'/all_parameters.Rdata',sep=''))


