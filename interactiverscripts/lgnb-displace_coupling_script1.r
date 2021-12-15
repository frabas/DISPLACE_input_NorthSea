#######################################################################################
#                                                                                     #
#                                       NORDFO                                        #
#           Coupling NS-LGNB model to DISPLACE: re-evaluate objective function        #
#                               Marie-Christine Rufener                               #
#                               Revised by F. Bastardie                               #
#                                                                                     #
#######################################################################################


# This is the very first script of the LGNB-DISPLACE coupling.
# In order to assure the coupling, we needto retrieve some information from the
# spatio-temporal correlation parameters. 
# To do so, we have first to re-evaluate the obj (objective function) from the LGNB result output,
# and then retrieve the spatio-temporal parameters and the hessian matrix.


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

 
memory.limit(size=56000) #To avoid memory problems 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Load libraries & functions, set working directory and load the results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.1) Load basic libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(TMB)
# first install https://cran.r-project.org/bin/windows/Rtools/rtools40.html
library(gridConstruct)    # devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct") # To install the gridConstruct package CAUTION: devtools is under rtools40 for R>4.0 
library(Matrix)
library(fields)
library(raster)
library(tidyr)


#  load the grid
load("F:/NORDFO_03122020/Results/Cod/No_predictors/results_Gadus_morhua_m1_L2_survey_No_alpha_YearQuarter_highRES_.RData") #Choose an arbitraty result object, as the grid (gr) is the same across species, years, and length groups
rm(list=setdiff(ls(), c("gr")))

  
 ##!!!!!!!!!!##
 ##!!!!!!!!!!##
 ##!!!!!!!!!!##
 sp <- "COD"
 ##!!!!!!!!!!##
 ##!!!!!!!!!!##
 ##!!!!!!!!!!##


# 1.2) Load the NS-LGNB model 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(file.path("F:","NORDFO_03122020")) ##set here your WD where the LGNB cpp file is stored

if(FALSE){
 #Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin") #Run only when on windows
 Sys.setenv(PATH="%PATH%;C:/Rtools_old/mingw64/bin;c:/Rtools_old/bin") #Run only when on windows
 compile(file.path("F:","NORDFO_03122020", "gNORDFO", "model.cpp")) # Compile model only when running script for the first time
 }
 dyn.load(dynlib(file.path("gNORDFO", "model"))) # Load the same C++ model used to generate the loaded results!


if(sp=="COD"){
model <- "m1"
datatype     <- c(L1="survey_No_alpha",L2="survey_No_alpha",L3="survey_No_alpha",L4="survey_No_alpha",L5="survey_No_alpha",L6="survey_No_alpha",L7="both_No_alpha", L8="both_No_alpha",L9="both_No_alpha",L10="both_No_alpha",L11="both_No_alpha",L12="both_No_alpha",L13="both_No_alpha") # caution: LG specific
species      <- "Cod"
spname       <- "Gadus_morhua"
LG_span      <- paste0("L", 1:13)
## Juveniles # see L50 to know
juveniles_span <- paste0("L", 1:6)
## Adults
adults_span <- paste0("L", 7:13)
a_comment <- ""
}

if(model %in% c("m1")) {
   predictor <-  "No_predictors" 
} else{
   predictor <-  "With_predictors" 
}

# 1.3) Load NS-LGNB model results for each size group 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("E:/NORDFO/Results/Cod/No_predictors/") #Note: default result that are loaded are those where no environmental predictors were used.
setwd(file.path("F:","NORDFO_03122020","Results",species, predictor)) #Results folder contains a folder for each species; default is Cod folder
dir.create(path=file.path(getwd(), paste0("LGNB-DISPLACE_coupling_",model)))
## Change the folder name according to species being evaluated - default is cod.
## You have to switch to the /With_predictors folder if you want to load the outputs with the predictions based on env. predicros.

# NOTE: Remember that we never have results beloging to L0 (first size group within DISPLACE, 0-5cm).
# NOTE2: For cod, we have results from L1-L13 - this changes for the other species!!
# NOTE3: It takes a long time to load all the results below!


for(LG in LG_span){
   print(LG)
   load(paste0("results_",spname,"_",model,"_",LG,"_",datatype[LG],"_YearQuarter_highRES_",a_comment,".RData"))
   obj <- env1$obj # Retrieve obj from model; remeber to switch env1 to env2/env3/env4... if results with env. predictors should be used (refer to very last part of the NS_LGNB_cod.R script for a reminder)
   obj$fn(env1$fit$par)
   lpb <- obj$env$last.par.best #Get spatio-temporal parameters
   r <- obj$env$random # Get latent (random) variables
   h <- obj$env$spHess(lpb, random=TRUE) #Hessian matrix based on the spatio-temporal models;This is where we need the loaded model.cpp file (done in section 1.2)
   #image(h)
   env1$lpb <- lpb
   env1$r   <- r
   env1$h   <- h
   
  
  # This is the second step of the LGNB-DISPLACE coupling, where we
  # extract important information from the LGNB model that will be used in the 
  # last LGNB-DISPLACE coupling script (script 2).
  # Specifically, we retrieve information from the: 
  # * spatio-temporal correlation parameters (phi, delta, scale)
  # * Precision matrix (Q)
  # * time periods (to identify the last time-period)
  # * estimated abundances for each time period (transformed to natural scale)
  time_corr <- env1$lpb["time_corr"] #untransformed time-correlation parameter
  delta <- exp(env1$lpb["logdelta"]) # delta (spatial correlation param.)

  lL <- list() # Put everything into a list
  lL$abundance <- as.data.frame(as.list(env1$sdr, "Estimate")[1]);
  lL$abundance <- apply(lL$abundance,2,exp) #Put on natural scale
  lL$gr <- gr
  lL$time_period <- env1$data$time

  lL$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
  lL$delta <- exp(env1$lpb["logdelta"]) # delta (spatial correlation param.)
  lL$scale <- exp(env1$lpb["logscale"]) # scale (spatial correlation param)
  lL$Q <- env1$obj$env$data$Q0 + delta * env1$obj$env$data$I # Precision matrix

  assign(LG, lL)
  save(list=c(LG), file=file.path(getwd(),"LGNB-DISPLACE_coupling", paste0(spname,"_",LG,"_obj_",model,".rds")))
  rm(env1, obj)
  gc()
  
  }



  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Make a list of all length-specific lists
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# re-load
for(LG in LG_span){
   print(LG)
   load(file=file.path(getwd(),"LGNB-DISPLACE_coupling", paste0(spname,"_",LG,"_obj_",model,".rds")))
}


fullres <- list(L1=L1, L2=L2, L3=L3, L4=L4, L5=L5, L6=L6, L7=L7, L8=L8,
                L9=L9, L10=L10, L11=L11, L12=L12, L13=L13)


# Set colnames based on timesteps
for(i in seq_along(fullres)){
  colnames(fullres[[i]]$abundance) <- levels(fullres[[i]]$time_period) 
}



# Save full list (to be used in script 2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save(fullres, file=file.path(getwd(),"LGNB-DISPLACE_coupling",paste0("Parameters_LengthGroups_",species,".RData"))) 



