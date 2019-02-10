library(SQUAREM)



compute_log10_BF<-function(bhat, sdhat, phi_vec){
	return(0.5*log10(sdhat^2/(sdhat^2+phi_vec^2))+0.5*((bhat/sdhat)^2*phi_vec^2/(sdhat^2+phi_vec^2))/log(10))
}

weighted_log10_BF<-function(log10_BFv, wv){

	max = max(log10_BFv)
	vec = log10_BFv - max
	rst = sum(10^vec*wv)
        return(log10(rst)+max)
}

compute_grid_wts<-function(log10_BF_vec, pi_vec, pi0,log10_NC){
	pi_vec[pi_vec<1e-10] = 1e-10
	pi_vec = pi_vec/sum(pi_vec)
	prob_vec = c(pi0, (1-pi0)*pi_vec)
	vec = c(0, log10_BF_vec)
	return(10^(log10(prob_vec)+vec-log10_NC))
}


estimate_wts<-function(wts_matrix, annot_vec, cat){
	sub_wts_matrix = wts_matrix[which(annot_vec==as.factor(cat)),]
	#sub_wts_matrix = wts_matrix
	wts = apply(sub_wts_matrix,2,mean)
	pi_vec = wts[2:length(wts)]/(1-wts[1])
	return(pi_vec)
}

  

torus.loglik<-function(p, BF_matrix, annot_vec){

	size = dim(BF_matrix)[2]
	sp = dim(BF_matrix)[1]
	cat_vec = levels(as.factor(annot_vec))
        enrich_param_size = length(cat_vec)
        index_vec = sapply(1:sp, function(x) which(cat_vec==as.factor(annot_vec[x])))
	pi_matrix = t(sapply(1:length(cat_vec),function(x) p[(1+enrich_param_size+(x-1)*size) : (enrich_param_size+x*size)]))
        
	log10_BF_vec = sapply(1:sp, function(x) weighted_log10_BF(BF_matrix[x,], as.vector(pi_matrix[index_vec[x],])))
	prior_vec  = sapply(1:sp, function(x) p[index_vec[x]])
	log10_NC_vec = sapply(1:sp, function(x) weighted_log10_BF(c(0,log10_BF_vec[x]), c(1-prior_vec[x], prior_vec[x])))
        loglik = sum(log10_NC_vec)/log10(exp(1))
	return(loglik)
}

torus_pool.loglik<-function(p, BF_matrix, annot_vec){

	size = dim(BF_matrix)[2]
	sp = dim(BF_matrix)[1]
	cat_vec = levels(as.factor(annot_vec))
        enrich_param_size = length(cat_vec)
        index_vec = sapply(1:sp, function(x) which(cat_vec==as.factor(annot_vec[x])))
	pi_vec = p[(1+enrich_param_size) : (enrich_param_size+size)]
        
	log10_BF_vec = sapply(1:sp, function(x) weighted_log10_BF(BF_matrix[x,], pi_vec))
	prior_vec  = sapply(1:sp, function(x) p[index_vec[x]])
	log10_NC_vec = sapply(1:sp, function(x) weighted_log10_BF(c(0,log10_BF_vec[x]), c(1-prior_vec[x], prior_vec[x])))
        loglik = sum(log10_NC_vec)/log10(exp(1))
	return(loglik)
}



torus.em<-function(p, BF_matrix, annot_vec){

	pnew = rep(NA, length(p))
	size = dim(BF_matrix)[2]
	sp = dim(BF_matrix)[1]


	cat_vec = levels(as.factor(annot_vec))
	enrich_param_size = length(cat_vec)
	index_vec = sapply(1:sp, function(x) which(cat_vec==as.factor(annot_vec[x])))
        pi_matrix = t(sapply(1:length(cat_vec),function(x) p[(1+enrich_param_size+(x-1)*size) : (enrich_param_size+x*size)]))


        log10_BF_vec = sapply(1:sp, function(x) weighted_log10_BF(BF_matrix[x,], as.vector(pi_matrix[index_vec[x],])))      
	prior_vec  = sapply(1:sp, function(x) p[index_vec[x]])
	log10_NC_vec = sapply(1:sp, function(x) weighted_log10_BF(c(0,log10_BF_vec[x]), c(1-prior_vec[x], prior_vec[x])))
	pip_vec = 1 - 10^(log10(1-prior_vec) - log10_NC_vec)
	wts_matrix = t(sapply(1:sp, function(x) compute_grid_wts(BF_matrix[x,],as.vector(pi_matrix[index_vec[x],]),1-prior_vec[x], log10_NC_vec[x])))
	wts = sapply(cat_vec, function(x) estimate_wts(wts_matrix, annot_vec, x) )
	pi_vec = sapply(cat_vec, function(x) sum(pip_vec[annot_vec == as.factor(x)])/length(which(annot_vec == as.factor(x))))
	 
      
	pnew = c(as.vector(pi_vec), as.vector(wts))
	return(pnew)

}
		      
torus_pool.em<-function(p, BF_matrix, annot_vec){

	pnew = rep(NA, length(p))
	size = dim(BF_matrix)[2]
	sp = dim(BF_matrix)[1]


	cat_vec = levels(as.factor(annot_vec))
	enrich_param_size = length(cat_vec)
	index_vec = sapply(1:sp, function(x) which(cat_vec==as.factor(annot_vec[x])))
        pi_vec = p[(1+enrich_param_size) : (enrich_param_size+size)]


        log10_BF_vec = sapply(1:sp, function(x) weighted_log10_BF(BF_matrix[x,], pi_vec))      
	prior_vec  = sapply(1:sp, function(x) p[index_vec[x]])
	log10_NC_vec = sapply(1:sp, function(x) weighted_log10_BF(c(0,log10_BF_vec[x]), c(1-prior_vec[x], prior_vec[x])))
	pip_vec = 1 - 10^(log10(1-prior_vec) - log10_NC_vec)
	wts_matrix = t(sapply(1:sp, function(x) compute_grid_wts(BF_matrix[x,],pi_vec,1-prior_vec[x], log10_NC_vec[x])))

	wts = as.vector(apply(wts_matrix,2,mean))
	nwts = wts[2:length(wts)]/(1-wts[1])
       	vec = sapply(cat_vec, function(x) sum(pip_vec[annot_vec == as.factor(x)])/length(which(annot_vec == as.factor(x))))
	          


	pnew = c(as.vector(vec), as.vector(nwts) )
	return(pnew)
}
	

torus<-function(betahat, sebetahat, annotation, tol = 1e-1, pool_est = TRUE){

	bhat_vec = betahat
	sd_vec = sebetahat
        annot_vec = annotation


	phi_min = min(sd_vec)/10
	phi_max = 2*sqrt(max(bhat_vec^2-sd_vec^2))

	if(phi_max<phi_min||is.na(phi_max)){
   		phi_max = 8*phi_min
	}	

	phi_vec = sort(exp(seq(log(phi_max), log(phi_min), by = -0.5*log(2))))
	phi_vec = c(phi_vec[1]/sqrt(2), phi_vec)
	pi_vec = rep(1.0/length(phi_vec), length(phi_vec))
	
	cat_vec = levels(as.factor(annot_vec))
	enrich_param_size = length(cat_vec)
	
	if(pool_est == FALSE){
	   p0 = c(rep(1e-3, enrich_param_size), rep(pi_vec, enrich_param_size))
        }else{
	   p0 = c(rep(1e-3, enrich_param_size), pi_vec)
        }	
	BF_matrix = t(sapply(1:length(bhat_vec), function(x) compute_log10_BF(bhat_vec[x], sd_vec[x], phi_vec)))
	
        if(pool_est == FALSE){
	      pf = squarem(p=p0, BF_matrix = BF_matrix, annot_vec = annot_vec, fixptfn=torus.em, objfn=torus.loglik, control=list(tol=tol))
         }else{
              pf = squarem(p=p0, BF_matrix = BF_matrix, annot_vec = annot_vec, fixptfn=torus_pool.em, objfn=torus_pool.loglik, control=list(tol=tol))
	}
        
	p = pf$par
	enrich_param_size = length(cat_vec)
        alpha_vec = p[1:enrich_param_size]
        alpha_vec = log(alpha_vec/(1-alpha_vec))
        a0 = alpha_vec[1]
	alpha_vec = alpha_vec - a0
        alpha_vec[1] = a0
        
	size = dim(BF_matrix)[2]
	sp = dim(BF_matrix)[1]
	index_vec = sapply(1:sp, function(x) which(cat_vec==as.factor(annot_vec[x])))

	if(pool_est == FALSE){
        	pi_matrix = t(sapply(1:length(cat_vec),function(x) p[(1+enrich_param_size+(x-1)*size) : (enrich_param_size+x*size)]))
        	log10_BF_vec = sapply(1:sp, function(x) weighted_log10_BF(BF_matrix[x,], as.vector(pi_matrix[index_vec[x],])))      

		effect_est =  t(pi_matrix)
	}else{
		pi_vec = p[(1+enrich_param_size) : (enrich_param_size+size)]
		log10_BF_vec = sapply(1:sp, function(x) weighted_log10_BF(BF_matrix[x,], pi_vec))  
	        effect_est = as.matrix(pi_vec,nrow=1)	
	}
        prior_vec  = sapply(1:sp, function(x) p[index_vec[x]])
	log10_NC_vec = sapply(1:sp, function(x) weighted_log10_BF(c(0,log10_BF_vec[x]), c(1-prior_vec[x], prior_vec[x])))
	lfdr = 10^(log10(1-prior_vec) - log10_NC_vec)
	

        enrichment_est = alpha_vec
	rst = list(enrichment_est = enrichment_est, effect_est = effect_est, em=pf,lfdr=lfdr)
	
	return(rst)

}




