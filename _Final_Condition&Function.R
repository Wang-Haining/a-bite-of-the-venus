
################################################################
########## Model Fit Evaluation for PCM & GPCM on LD  ##########
################################################################
## PCM & GPCM, thresholds are item-specific
## Simulations factors - General: 
## model = [PCM, GPCM], sample = [500,1000], item = [10,20], category = [3,4,5]
## Simulations factors - LD-specific: 
## nstrait = [0,1,2,3], corr = [.05, .5], pctLD = [.2,.3,.4]
## total conditions = 2*2*2*3 * (1 + 3*2*3) = 456
## data generated from multidim GPCM (a=1 for MPCM) then fitted to Unidim models
################################################################
################################################################

####################### Create simulation conditions #######################
ModelType = c('PCM', 'GPCM'); SampleSize = c(500, 1000); TestLength = c(10, 20)
ItemCategory = c(3,4,5); NStrait = c(0,1,2,3); Corr = c(.05,.5); pctLD = c(.2,.3,.4)
BASE = LD = CONDITION = NULL # Baseline conditions when nstrait=0
for (nt in NStrait) {# Baseline conditions come first
  for (mt in ModelType) {# PCM vs. GPCM
    for (ss in SampleSize) {# Sample size: 500 vs. 1000
      for (tl in TestLength) {# Test length: 10 vs. 20
        for (ic in ItemCategory) {# Item response categories: 3, 4, 5
          for (cr in Corr) {# Correlation(primary, nuisance): .05 vs. .5
            for (pld in pctLD) {# %.LD items on each nuisance trait: .2, .3, .4
              if(nt == 0)
                BASE = c('ModelType'=mt, 'SampleSize'=ss, 'TestLength'=tl, 
                         'ItemCategory'=ic, 'NStrait'=nt, 'Correlation'=0, 'pctLD'=0)
              else
                LD = c('ModelType'=mt, 'SampleSize'=ss, 'TestLength'=tl, 
                       'ItemCategory'=ic, 'NStrait'=nt, 'Correlation'=cr, 'pctLD'=pld)
              ## 456 by 7
              CONDITION = unique(data.frame(rbind(CONDITION, BASE, LD), row.names = NULL))
            }
          }
        }
      }
    }
  }
}
CONDITION = CONDITION %>% mutate_at(vars(names(CONDITION)[-1]), as.numeric)
CONDITION$Condition = c(replicate(24, 'Baseline'), replicate(nrow(CONDITION)-24, 'LD'))
#CONDITION$ModelType = factor(CONDITION$ModelType, labels = c('PCM', 'GPCM'))
### Save the condition file to project directory
write.csv(CONDITION, paste(project.dir, '/_ALL_Conditions_',nrow(CONDITION), ".csv",sep=""))

##############################################################
###### Read conditions
#CONDITION = read.csv(paste(project.dir, list.files(project.dir)[1], sep = ''), row.names = 1)

##############################################################
##############################################################
############# Functions to Obtain Overall Results #############
##############################################################
Results = function(cd, iteration){# specify condition index [cd] and #.iterations [it]
  ##########################
  ## function to generate discrimination
  GenA = function(){
    ## the primary theta for two model conditions
    if(md == 'PCM')
      a.p = rep(1, test)# discrimination=1 on primary theta for the PCM
    if (md == 'GPCM')
      a.p = runif(test, .8, 1.5)# discrimination~unif(.8, 1.5) on primary theta for the GPCM
    
    ## the nuisance theta's
    a.nt = matrix(0, test, nstrait)# a.nt=1 on nuisance theta for LD items only
    
    if(nstrait == 0){# baseline without LD
      discrimination = data.frame('a.p'=a.p)
    }
    if(nstrait > 0)# LD conditions
      for(nt in 1:nstrait){
        # selected LD items: select #.items given pctLD for each nuisance trait
        ld.item = sample(test, round(pctLD * test), replace = F)
        a.nt[ld.item, nt] = 1 # assign a=1 on nuisance traits for LD items
        discrimination = data.frame('a.p'=a.p, 'a.nt'=a.nt)
      }
    return(discrimination)
  }# end GenA()
  
  ##########################
  ## Function to generate item location/steps
  GenB = function(){
    t = matrix(NA, test, (category-1))
    for (i in 1:test) {
      t[i,] = runif((category-1),.8, 2)
    }
    t1 = t(apply(t, 1, cumsum))
    t2 = t1 - rowMeans(t1)
    t3 = t2 + rnorm(test,0,.5)
    step = data.frame('step'=t3)
    return(step)
  }# end GenB()
  
  ##########################
  ## Function to generate thetas
  GenTheta = function(){
    # the primary theta
    theta1 = matrix(rnorm(ss), ss, 1)
    
    if(nstrait == 0)# baseline without LD
      return(theta = data.frame('theta.p' = theta1))
    
    if(nstrait > 0){# LD conditions
      gencorr = function(){
        y = rnorm(ss)
        res = resid(lm(y ~ theta1))
        z = corr * scale(theta1) + sqrt(1-corr^2) * scale(res)
        return(z)
      }
      theta.nt = matrix(NA, ss, nstrait)
      for (nt in 1:nstrait) {
        for (r in corr) {
          theta.nt[,nt] = gencorr()
        }
      }
      theta = data.frame('theta.p'=theta1, 'theta.nt'=theta.nt)
      return(theta)
    }
  }# end GenTheta()
  
  ################# ITERATION STARTS HERE ##############
  ## Function to generate item responses
  GenResp = function(){
    response = matrix(NA, ss, test)# response matrix
    for (n in 1:ss) {
      for (i in 1:test) {
        u = runif(1, 0, 1)# a random number for discretization

        e = NA
        e[1] = exp(0)
        for (k in 2:category) {
          e[k] = exp(sum(discrimination[i,] * theta[n,]) * (k-1) - sum(step[i, c(1:k-1)]))
        }
        e.sum = sum(e)
        crf = e / e.sum # crf of the item
        cum.crf = cumsum(crf) # cumulative crf
        
        for (z in 1:category) {
          if(u <= cum.crf[z]){response[n,i] = (z-1)
          break}
        }
      }
    }
    output = data.frame(response)
    return(output)
  }# end GenResp()
  
  #############################################################################
  #### Functions to fit the response to IRT models and obtain results ####
  FITmodels = function(){
    # fit IRT model and obtain parameter estimates, RMSEs & M2-based GFIs
    if(md == 'PCM')
      irtfit = mirt(Response, 1, itemtype = 'Rasch', method = 'EM', verbose=FALSE)
    if(md == 'GPCM')
      irtfit = mirt(Response, 1, itemtype = 'gpcm', method = 'EM', verbose=FALSE)

    thetahat = fscores(irtfit, method = "EAP", full.scores = TRUE) # estimated theta
    discriminationhat = matrix(coef(irtfit, IRTpars = T, simplify = T)$item[,1]) # category-1 columns
    stephat = coef(irtfit, IRTpars = T, simplify = T)$item[,-1] # category-1 columns
    
    #######################
    ## Function to compute RMSEs of parameter estimates
    Recovery = function(true, estimate){
      difference = na.omit(estimate - true)
      bias = colMeans(difference)
      sqrdiff = difference^2
      rmse = colMeans(sqrdiff)
      return(list('bias'= bias, 'rmse' = rmse))
    }# end Recovery()
  
    
    ######################## IRT-Based results ########################
    r.theta = cor(theta[,'theta.p'], thetahat, method = 'pearson')[1]# cor. of theta recovery
    theta.out = Recovery(theta[,'theta.p'], thetahat)
    discrimination.out = Recovery(discrimination[,'a.p'], discriminationhat)
    step.out = Recovery(step, stephat)
    
    rmse.theta = theta.out$rmse[[1]]; bias.theta = theta.out$bias[[1]]
    rmse.discrimination = discrimination.out$rmse; bias.discrimination = discrimination.out$bias
    rmse.step = step.out$rmse; bias.step = Recovery(step, stephat)$bias
    
    rmse.irt = c('corr.theta'=r.theta, 'rmse.theta'=rmse.theta, 'rmse.a'=rmse.discrimination,
                 'rmse'=rmse.step, 'rmse.step_avg'=mean(rmse.step),
                 'bias.theta'=bias.theta, 'bias.a'=bias.discrimination,
                 'bias'=rmse.step, 'bias.step_avg'=mean(bias.step))
    
    ## IRT model fit indices: M2* for polytomous data
    fitind.irt = M2(irtfit, type = 'C2')[c(1:3, 9, 8, 4, 7)] # C2, CFI, TLI, RMSEA, & SRMSR
    fitind.irt.m2 = M2(irtfit, type = 'M2*')[c(1:3, 9, 8, 4, 7)]# M2*
    # the C2-based GFIs has same df as CFA GFIs, while M2*-based has less df
    # record M2*-based GFIs for comparison purpose

    ## Collect IRT results from this condition ###
    irtind.m2 = t(c('Iter'=it, fitind.irt.m2))
    fit.irt = t(c('Iter'=it, fitind.irt, rmse.irt))
    
    #############################################################################
    #### Functions to fit the response to CFA models and obtain results ####
    ## Function to specify CFA models: depend on test length
    mdCFA = function(data){
      x = substring(colnames(data), 1)# extract variable names
      # make it into linear combinations
      if(md == 'PCM')
        x1 = paste('a*', x, sep='')# constrain all factor loadings = 1
      
      if(md == 'GPCM')
        x1 = paste('a', seq(1:length(x)), '*', x, sep='') 
          
      x2 = paste('+', x1[-1], sep='')
      x3 = paste(c('F=~', x1[1], x2), collapse ='')# constrain loadings=1
      x4 = paste(x3, 'F~~F', sep = '\n')# free factor variance est.
      # the model should be: F =~ 1*X1...+1*Xn & F~~F for PCM
      # the model should be: F =~ a1*X1...+an*Xn & F~~F for GPCM
      md.cfa=x4
      return(md.cfa)
    }# end mdCFA()

    ## Fit single-factor CFA model and obtain parameter estimates, RMSEs & GFIs:call mdCFA()
    cfafit = cfa(mdCFA(Response), Response, ordered = T,
                 estimator = 'WLSMV', parameterization = 'theta')
    est = parameterestimates(cfafit)# col.5 are the estimates
    fhat = lavPredict(cfafit)# factor score estimates
    
    y1 = test + 2; y2 = y1 + (category-1)*test - 1
    loading = est[1:test, 5]# factor loading estimates
    est_values = est[y1:y2, 5]# tau estimates
    stephat = matrix(est_values, test, category-1, byrow = T)
    ## CFA-Based RMSEs
    r.theta = cor(theta[,'theta.p'], fhat, method = 'pearson')[1]
    theta.out = Recovery(theta[,'theta.p'], fhat)
    discrimination.out = Recovery(discrimination[,'a.p'], matrix(loading))
    step.out = Recovery(step, stephat)
    
    rmse.theta = theta.out$rmse[[1]]; bias.theta = theta.out$bias[[1]]
    rmse.discrimination = discrimination.out$rmse; bias.discrimination = discrimination.out$bias
    rmse.step = step.out$rmse; bias.step = Recovery(step, stephat)$bias
    
    rmse.cfa = c('corr.theta'=r.theta, 'rmse.theta'=rmse.theta, 'rmse.a'=rmse.discrimination,
                 'rmse'=rmse.step, 'rmse.step_avg'=mean(rmse.step),
                 'bias.theta'=bias.theta, 'bias.a'=bias.discrimination,
                 'bias'=rmse.step, 'bias.step_avg'=mean(bias.step))
    
    ## CFA-based GFIs: use scaled results
    id = c('chisq.scaled', 'df.scaled', 'pvalue.scaled', 
           'cfi.scaled', 'tli.scaled', 'rmsea.scaled', 'srmr')
    fitind.cfa = fitmeasures(cfafit, id)
    ##### Collect CFA results from this condition ###
    fit.cfa = t(c('Iter'=it, fitind.cfa, rmse.cfa))

    #############################################################################
    #### Return IRT and CFA results
    return(list('irt.M2*'=irtind.m2, 'fit.irt'=fit.irt, 'fit.cfa'=fit.cfa))
  }# end FITmodels()
  ####################################################################
  ################# ITERATION ENDS HERE ##############################
  ####################################################################
  
  ## Function to average results across iterations
  AVG = function(reps){
    file.out = data.frame(matrix(unlist(reps), nrow = iteration))
    names(file.out) = colnames(reps)
    file.avg = c('Iter.'=nrow(reps), colMeans(file.out[,-1]))
    return(file.avg)
  }# end AVG()

  
  ## Function to compute the associations between RMSEs and GFIs
  Association = function(){
    dimnm = c(1:4)
    doc = list(IRT_M2, IRT, CFA)
    

    r.irt = p.irt = NULL
    for (d in 1:2) {
      rmseparam = doc[[2]][,c(10:(12+category-1))]
      fitind = doc[[d]][,c(5:8)]
      asso.r.irt = rcorr(fitind, rmseparam)$r[dimnm, -dimnm] # Pearson's r
      asso.p.irt = rcorr(fitind,rmseparam)$P[dimnm, -dimnm] # P-values for r
      r.irt = rbind(r.irt, asso.r.irt)
      p.irt = rbind(p.irt, asso.p.irt)
    }
    asso.irt = rbind('IRT_M2.correlation' = ' ', r.irt[dimnm,], ' '=' ', 
                     'IRT_M2.p-value'='  ', p.irt[dimnm,], 
                     ' ' = ' ', '****'='****', '  ' = '  ',
                     'IRT.correlation' = ' ', r.irt[-dimnm,], ' '=' ', 
                     'IRT.p-value'='  ', p.irt[-dimnm,])

    for (d in 3) {
      rmseparam = doc[[d]][,c(10:(12+category-1))]
      fitind = doc[[d]][,c(5:8)]
      asso.r.cfa = rcorr(fitind, rmseparam)$r[dimnm, -dimnm] # Pearson's r
      asso.p.cfa = rcorr(fitind,rmseparam)$P[dimnm, -dimnm] # P-values for r
      asso.cfa = rbind('CFA.correlation' = ' ', asso.r.cfa, ' '=' ', 
                       'CFA.p-value'='  ', asso.p.cfa)
    }
    ASSO = rbind(asso.irt, ' ' = ' ', '****'='****', '  ' = '  ', asso.cfa)
    return(ASSO)
  }# end Association()
  
  ## Function to compute PGFIs across iterations
  ## PGFIs: proportion of good-fit model indication under specific condition
  ## Use the data frames: PGFI_abs & PGFI_inc
  PGFIs = function(){
    ###### create a data frame storing PGFIs
    ## cutoffs: 7
    cutoff_abs = c(.01, .02, .03, .04, .05, .07, .09); cutoff_inc = 1 - cutoff_abs
    columns_abs = paste('<',cutoff_abs, sep = '')
    columns_inc = paste('>',cutoff_inc, sep = '')
    rows_abs = 'Absol.Ind'; rows_inc = 'Incre.Ind'
    
    PGFI_abs = data.frame(matrix(NA, ncol=length(cutoff_abs)+1),
                          row.names = paste(rows_abs, sep = ''))
    PGFI_inc = data.frame(matrix(NA, ncol=length(cutoff_inc)+1),
                          row.names = paste(rows_inc, sep = ''))
    PGFI_abs[1,] = c(columns_abs, 'df');  PGFI_inc[1,] = c(columns_inc, 'df')
    
    PGFI_abs = list('IRT_M2'=PGFI_abs, 'IRT'=PGFI_abs,'CFA'=PGFI_abs)
    PGFI_inc = list('IRT_M2'=PGFI_inc, 'IRT'=PGFI_inc,'CFA'=PGFI_inc)
    
    ###### IRT records ######
    doc = list(IRT_M2, IRT, CFA)
    for (d in 1:2) {
      for (ind in c('RMSEA', 'SRMSR')) {
        for (c in 1:length(cutoff_abs)) {
          gfi = doc[[d]][(doc[[d]][,ind] < cutoff_abs[c]), ind]
          PGFI_abs[[d]][ind,c] = length(gfi) / nrow(doc[[d]])
        }
      }
      for (ind in c('CFI', 'TLI')) {
        for (c in 1:length(cutoff_inc)) {
          gfi = doc[[d]][(doc[[d]][,ind] > cutoff_inc[c]), ind]
          PGFI_inc[[d]][ind,c] = length(gfi) / nrow(doc[[d]])
        }
      }
      PGFI_abs[[d]][-1,ncol(PGFI_abs[[d]])] =
        PGFI_inc[[d]][-1,ncol(PGFI_inc[[d]])] = doc[[d]][,'df'][1]
    }
    
    ###### CFA records ######
    for (d in 3) {
      for (ind in c('rmsea.scaled', 'srmr')) {
        for (c in 1:length(cutoff_abs)) {
          gfi = doc[[d]][(doc[[d]][,ind] < cutoff_abs[c]), ind]
          PGFI_abs[[d]][ind,c] = length(gfi) / nrow(doc[[d]])
        }
      }
      for (ind in c('cfi.scaled', 'tli.scaled')) {
        for (c in 1:length(cutoff_inc)) {
          gfi = doc[[d]][(doc[[d]][,ind] > cutoff_inc[c]), ind]
          PGFI_inc[[d]][ind,c] = length(gfi) / nrow(doc[[d]])
        }
      }
      PGFI_abs[[d]][-1,ncol(PGFI_abs[[d]])] = 
        PGFI_inc[[d]][-1,ncol(PGFI_inc[[d]])] = doc[[d]][,'df.scaled'][1]
    }
    #### Combine results
    ABS_Ind = rbind('Absolute' = 'Absolute', 
                    'IRT.m2'=' ', PGFI_abs[[1]], ' '=' ', 'IRT'=' ', 
                    PGFI_abs[[2]], '  '=' ', 'CFA'=' ', PGFI_abs[[3]])
    
    INC_Ind = rbind('Incremental' = 'Incremental',
                    'IRT.m2 '=' ', PGFI_inc[[1]], '   '=' ', 'IRT '=' ', 
                    PGFI_inc[[2]], '    '=' ', 'CFA '=' ', PGFI_inc[[3]])
    
    PGFI = rbind(ABS_Ind, '     '=' ', '****'='****', '      '=' ', INC_Ind)
    
    names(PGFI) = c(paste('cutoff', c(1:7), sep='_'), 'DF')
    return(PGFI)
  }# end PGFIs()
  ################ End of the architecture ##################
  
  ################################################################
  ################################################################
  #### Incorporate simulation conditions ####
  # cd is the condition index
  md = CONDITION[cd, 'ModelType']# model type
  ss = CONDITION[cd,'SampleSize']# sample size
  test = CONDITION[cd,'TestLength']# test length
  category = CONDITION[cd,'ItemCategory']# number of response categories
  nstrait = CONDITION[cd,'NStrait']# number of nuisance trait
  corr = CONDITION[cd,'Correlation']# residual correlation
  pctLD = CONDITION[cd,'pctLD']# percentage of LD items

  set.seed(cd)
  
  ## Generate model parameters for each condition
  discrimination = GenA()
  step = GenB()
  theta = GenTheta()
  
  ## Empty files to collect replication records and averaged results
  IRT_M2 = IRT = CFA = OUT = c()
  #pbi = progress_bar$new(total = iteration, clear = F, width = 100,
   #                      format = "(:spin) [:bar] :percent 
    #                       [Elapsed: :elapsedfull || Estimated remaining: :eta]")
  
  ##### Start iterations for response data and modeling results
  for (it in 1:iteration) {
    set.seed(it)# set.seed for each iteration
    # Adding a progress bar
    #pbi$tick()
    
    # Generate response data
    Response = GenResp()
    ## collect replication results
    ResComb = FITmodels()
    irt_M2 = ResComb$irt.M2
    IRTfit = ResComb$fit.irt
    CFAfit = ResComb$fit.cfa
    #### save iteration history in the records folder #####
    IRT_M2 = rbind(IRT_M2, irt_M2)
    IRT = rbind(IRT, IRTfit)
    CFA = rbind(CFA, CFAfit)
    
    cdname = paste(cd,paste(names(CONDITION),CONDITION[cd,], sep='_', collapse = '_'), sep='_')
    write.csv(IRT_M2, paste(records.dir,'/IRT_M2_', cdname, '.csv', sep=''))
    write.csv(IRT, paste(records.dir,'/IRT_', cdname, '.csv', sep=''))
    write.csv(CFA, paste(records.dir,'/CFA_', cdname, '.csv', sep=''))
  }# end iteration
  
  ## Compute averaged results across iterations
  OUT = data.frame(rbind('IRT'=AVG(IRT), 'CFA'=AVG(CFA)))
  OUT['IRT_M2*', c(1:length(AVG(IRT_M2)))] = AVG(IRT_M2)
  ASSO = Association()
  PGFI = PGFIs()
  
  write.csv(OUT, paste(output.dir,'/AVG_', cdname, '.csv', sep=''))
  write.csv(PGFI, paste(output.dir,'/PGFI_', cdname, '.csv', sep=''))
  write.csv(ASSO, paste(output.dir,'/ASSO_', cdname, '.csv', sep=''))
  
  #Sys.sleep(0.1)
  #setTxtProgressBar(pbi, it)
}# end Results(cd, iteration)
############### End of All Functions #############


