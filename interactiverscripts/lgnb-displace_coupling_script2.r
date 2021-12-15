#######################################################################################
#                                                                                     #
#                                     NORDFO                                          #
#            Coupling NS-LGNB model to DISPLACE: extract model parameters             #
#                             Marie-Christine Rufener                                 #
#                             Revised by F. Bastardie                                 #
#######################################################################################


# This is the last script of the LGNB-DISPLACE coupling.
# In this script, we used the information retrieved from the previous coupling script
# and use this information to reconstruct the spatio-temporal abundance field for one time-step ahead.
# This script is divided into 5 sections, namely:

# SECTION 1: Load model and model results 
# SECTION 2: Keep original abundance fields from model 
# SECTION 3: Predict abundance field one time-step ahead (core section for DISPLACE coupling)
# SECTION 4: Exract abundances for DISPLACE grid
# SECTION 5: Output for DISPLACE


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

print(R.version)


if(FALSE){
  # reload some utils from gridConstruct for making image() on signature PolygonGrid work
  image.polygons <- function(x,y,add=FALSE,col=heatpal(16), ##col=heat.colors(16),
                           codes=cut(y,length(col)),border=FALSE,
                           strip=colorstrip(col),
                           map=polygon(makeMap(x,map="world")),
                           breaks=NULL,...){
   if(!is.null(breaks)){
     stopifnot(length(breaks)==length(col)+1)
     codes <- cut(y,breaks)
   }
   if(!add)plot(as.data.frame(lapply(x,range)),type="n")
   polygon(splitPolygons(x),col=as.character(col[codes]),border=border,...)
   eval(map)
   eval(strip)
   } 
  # Model matrix.
  ## Assign the average values at the corners to each polygon. model.matrix.polygonGrid <- function(object,...){
  model.matrix.polygonGrid <- function(object,...){
   nsides <- attr(object,"nsides")
   fac <- attr(object,"polygons")
   n <- length(fac)/nsides
   i <- rep(1:n,each=nsides)
   j <- as.numeric(fac)
   A <- spMatrix(n,nrow(object),i=i,j=j,x=0*j+1/nsides)
   A
  }
  image.polygonGrid <- function(x,y,...){
   A <- model.matrix(x)
   polygons <- as.polygons(x)
   image(polygons,as.vector(A%*%y),...)
  }
  
 concTransform <-
function (x) 
{
    i <- order(x)
    ys <- sort(exp(x))  
    p <- ys/sum(ys)
    x[i] <- cumsum(p)
    x
}

  concTransform2 <-
function (x) 
{
    i <- order(x)
    #ys <- sort(exp(x))  # uncomment this, because exp() already done
    ys <- sort(x)
    p <- ys/sum(ys)
    x[i] <- cumsum(p)
    x
}
  
concTransform3 <- function(x) # Based on Bartolino et al. (2011)
{
  i <- order(x)
  ys <- sort(x)
  p <- ys/max(ys)
  x[i] <- p
  x
  #p
}
  

} # end FALSE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load standard inputs for DISPLACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FIXME: Caution, the lines below have been adapted to the North Sea DISPLACE application!


# CALLED FROM DISPLACE when dyn_pop_sce.option(Options::nbcpCoupling)
args <- commandArgs(trailingOnly = TRUE)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied to nbcpCoupling.r: Take pop 0 and tstep 0 and sce scelgnbcoupling and simu simu1")
  pop     <- 0
  tstep   <- 745
  sce     <- "scelgnbcoupling"
  sim     <- "simu1"
  igraph  <- "a_graph1"
}else{
  pop      <- args[1]
  tstep    <- args[2]
  sce      <- args[3]
  sim      <- args[4]
  igraph   <- args[5]
}



igraph <- as.numeric(gsub("a_graph", "", igraph))

application <- "NorthSea"

FRANCOIS <- TRUE ; MARIECHRISTINE <- FALSE



if(.Platform$OS.type == "unix") {
  if(FRANCOIS){
    path <- file.path("~","ibm_vessels", paste0("DISPLACE_input_", application))
  }
  if(MARIECHRISTINE){
    path <- file.path("~:","PHD_projects","Proj_3",paste0("DISPLACE_input_", application))  
  }
}
if(.Platform$OS.type == "windows") {
  if(FRANCOIS){
    path <- file.path("D:", "FBA", paste0("DISPLACE_input_", application))
  }
  if(MARIECHRISTINE){
    path <- file.path("C:","Users","mruf","Documents","PHD_projects","Proj_3",paste0("DISPLACE_input_", application))
  }
}

cat(paste("looking at the R script in ", path, "\n"))




#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  
if(TRUE){
  
  
  # Load basic libraries
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #library(TMB)
  #library(gridConstruct)
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(fields)))
  suppressWarnings(suppressMessages(library(raster)))
  suppressWarnings(suppressMessages(library(tidyr)))

 
  # SECTION 1: Predict abundance field one time-step ahead (core section for DISPLACE coupling)
  # SECTION 2: Extract abundances for DISPLACE grid
  # SECTION 3: Output for DISPLACE
  
  

  # Multivariate normal distribution simulation function (based on precision rather variance)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # @ mu = mean of the abundance field
  # @ prec = precision of the spatio-temporal covariance matrix
  # @ n.sims = number of simulations
  rmvnorm_prec <- function(mu, prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- Cholesky(prec, super=TRUE)
    z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  
  
  
  
  
  #><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
  ########################                SECTION 1                ########################
  #><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
  
  
  # Load results
  #~~~~~~~~~~~~~~~~
  if(pop==4) a_species_name <- "Cod"  # see pop_names_NorthSea.txt
  if(pop==15) a_species_name <- "Herring"
  if(pop==26) a_species_name <- "Plaice"
  if(pop==35) a_species_name <- "Sole"
  if(pop==36) a_species_name <- "Sprat"
  if(pop==43) a_species_name <- "Whiting"
  #...
 
  load(file.path(path, "interactiverscripts", paste0("Parameters_LengthGroups_",a_species_name,".RData"))) #fullres is a list of lists, where each individual list stores the results of a particular size-group
  # names(fullres) #names of the lists(related to the DISPLACE size-groups)
  # names(fullres[[1]]) #names of the objects present in each list    
  if(FALSE){
    plot(gr,type="n",xlab="Longitude",ylab="Latitude",las=1,xlim=c(-3,11), ylim=c(48,60)) 
    fullres[[1]]$abundance[,1][is.infinite(fullres[[1]]$abundance[,1])] <- 0  # remove Inf
    a_sizegr <- "L10"
    a_tstep  <-  "2019 3"
    dd <- cbind(gr, concTransform2(fullres[[a_sizegr]]$abundance[,a_tstep]  ))
    excluded_cells <- which(dd[,3]<0.0)
    fullres[[a_sizegr]]$abundance[excluded_cells,a_tstep]  <- 0
    image.polygonGrid(gr, concTransform3(fullres[[a_sizegr]]$abundance[,a_tstep]  ),
           add=TRUE)
  # =>To see an example of the estimated abundance fields
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1.1) Retrieve spatio-temporal parameters
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Get abundances for the last time-step (2019, with all 4 quarters)
  # FIXME: fix line below if model was run on a temporal resolution different than Year-Quarter basis
  abundances <- vector("list", length(fullres))
  Q1 <- ncol(fullres[[1]]$abundance)-4+1
  Q4 <- ncol(fullres[[1]]$abundance)
  for(i in seq_along(fullres)){
    abundances[[i]] <- fullres[[i]]$abundance[,Q1:Q4]
  }
  names(abundances) <- names(fullres)
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1.2) Predict one-time step ahead
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # To predict the abundance fields in t+1, we need to "draw" a new abundance field that is based
  # on the spatio-temporal correlation parameters + precision matrix that were estimated for each stock-sizeGroup
  # (retrieved in the coupling script part 1 (LGNB_DISPLACE_coupling_script1.R)).
  
  # As such, we need to apply the rmvnorm_prec function, which basically simulates
  # a multivariate Gaussian random distribution based on the input spatio-temporal corr. parameters and precision matrix.
  # The function requires that the user specifies the number of simulations; usually,the more precise the abundance
  # estimates, the more precise will be the forward prediction (thus, only a few simulations, or even only 1, would be required).
  # However, when abundance estimates have high uncertainties (such as those of very small size-groups), the forward-predicted
  # abundance fields tend to be very different from one simulation to another, and therefore a higher number of simulation
  # is required to stabilize the predicted abundance field. Provided that it is a "half-simulation" framework,
  # we can compute forward-predictions for all four quarters of the last time-period. 
  
  
  Allspred <- vector("list", length(fullres))
  
  
  NROW <- dim(fullres[[1]]$abundance)[1]
  NCOL <- ncol(abundances[[1]])
  DIM <- length(fullres)
  tmp <- array(1:(NCOL*NROW), c(NROW,NCOL,DIM)) 
  
  #NSIM <- 1000 #Choose the number of simulations and take the mean across the simulated fields
  NSIM <- 2 # if coupled with DISPLACE we only need one simu, as a new simu will be called at each time step of the coupling. (we put NSIM at 2 to avoid the rowMeans to fail...)
  
  
  
  for(i in seq_along(abundances)){
    for(j in 1:NCOL){ # per year-quarter
      Abundance <- abundances[[i]][,j]
      #tmp[,j,i] <- rmvnorm_prec(mu = Abundance, prec = fullres[[i]]$scale*sqrt(1-fullres[[i]]$phi^2)*fullres[[i]]$Q, n.sims = 1) #Simulate and store results in the array
      tmpres <- as.data.frame(rmvnorm_prec(mu = Abundance, prec = fullres[[i]]$scale*sqrt(1-fullres[[i]]$phi^2)*fullres[[i]]$Q, n.sims = NSIM)) #Simulate and store results in the array
      tmpres <- rowMeans(tmpres[,1:NSIM])
      tmp[,j,i]  <- tmpres
      Allspred <- lapply(seq(dim(tmp)[3]), function(x) tmp[ , , x]) # Convert array back to list (easier to manipulate)
    }
  }
  
  
  # Set abundances back to natural scale
  # for(i in seq_along(Allspred)){
  #   Allspred[[i]] <- apply(Allspred[[i]],2,exp)
  # }
  # ...but too biased toward very high value if tramsned in natural scale...so, keep the exp scale and transform by replacing by:
    for(i in seq_along(Allspred)){
     Allspred[[i]] <- apply(Allspred[[i]],2,function(x) x+(-1*min(x)))
   }  # to remove the negative values...
 
  
  
  names(Allspred) <- names(fullres)
  
  ## Plot to see the progress
  # image(cbind(fullres[[10]]$gr, z=concTransform(Allspred[[10]][,1])), col=tim.colors(99)) #plotting example 
  # replaced by: 
  #plot(gr,type="n",xlab="Longitude",ylab="Latitude",las=1,xlim=c(-3,11), ylim=c(48,60)) 
  #image.polygonGrid(gr, concTransform3(Allspred[[10]][,1]),add=TRUE)

  
  cat("simulate a nbcp pop ",pop," abundance field for t+1...done\n")
  
  
  
 
  #><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
  ########################                SECTION 2                ########################
  #><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
  
  
  
  # First we need to create individual raster files, as the abundance output given in 
  # SECTION 1 isn't a raster file. After that we use the lon/lat of the DISPLACE graphs 
  # and extract the abundances according to each graph node
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.1) Create an empty rasterfile 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gr <- fullres[[1]]$gr # Get the grid from an arbitrary size-class (gr is the same throughout all size-classes)
  e <- extent(as.matrix(gr))
  r <- raster(e,ncol=55,nrow=55) #Adapt raster size according to grid size; These numbers are good for 5x5km grid (as model output)
  
  
  # Include grid positions to dataframes
  coords <- data.frame(lon=gr$lon,lat=gr$lat) #Retrieve grid positions from model output
  Allspred <- lapply(Allspred, cbind, coords) #Bind grid positions
  
  
  # Convert df from wide to long format
  Allspredw <- list(NULL)
  for(i in seq_along(Allspred)){
    Allspredw[[i]] <- gather(Allspred[[i]], key=Quarter,value=Abundance,-c(lon,lat))
  }
  names(Allspredw) <- names(Allspred) #name new list accordingly
  
  
  # Include size-specific column to each df
  for( i in seq_along(Allspredw)){
    Allspredw[[i]]$Size <- rep(names(Allspredw)[i],nrow(Allspredw[[i]]))
  }
  
  Allspredw <- do.call("rbind", Allspredw) #Unlist
  rownames(Allspredw) <- NULL
  
  
  # Identify properly all quarters
  Allspredw$Quarter <- ifelse(Allspredw$Quarter=="1","Q1",
                              ifelse(Allspredw$Quarter=="2","Q2",
                                     ifelse(Allspredw$Quarter=="3","Q3","Q4")))
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.2) Create raster based on abundances
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dflist <- split(Allspredw,list(Allspredw$Quarter,Allspredw$Size))
  
  abufields <- list() 
  for(i in 1:length(dflist)){
    abufields[[i]] <- rasterize(dflist[[i]][,c("lon","lat")], r, dflist[[i]][,"Abundance"], fun=mean)
    abufields[[i]] <- disaggregate(abufields[[i]],2,method="bilinear")
  }
  names(abufields) <- names(dflist) #name new list accordingly
  
  #image(abufields[[47]],col=tim.colors(99)) 
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2.3) Extract abundances based on DISPLACE graph coords
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #loc <- read.table(file="H:/FB_MR/Coupling/coord40.dat",sep="") # DISPLACE coords
  #loc <- read.table(file="C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Manuscript_03/Data/DISPLACE/coord40.dat",sep="") # DISPLACE coords
  loc <- read.table(file=file.path(path, "graphsspe", paste("coord",igraph,".dat",sep=""))) # DISPLACE coords
  loc <- as.matrix(as.vector(loc))
  loc <- matrix(loc, ncol=3)
  loc <- cbind(loc, 1:nrow(loc))
  colnames(loc) <- c('lon', 'lat', 'harb', 'pt_graph')
  
  loc <- as.data.frame(loc)
  loc <- loc[loc$harb==0,] # remove points on harbours (all harb>0)
  loc2 <- loc[,c("lon","lat")] # keep only lon/lat columns
  
  
  
 
  

  # Francois
  #  if(pop==2){
  #   # clip to the real stock area definition
  #   clip <- raster(file.path(path, "interactiverscripts", "sd222324.tif"))
  #   for(i in seq_along(abufields)){
  #     abufields[[i]] <- crop(abufields[[i]], clip)
  #   }
  # }
  # 
  
  dfa <- list()
  for(i in seq_along(abufields)){
    dfa[[i]]<- raster::extract(abufields[[i]],loc2)
    dfa[[i]]<- data.frame(abundance=dfa[[i]], lon=loc2$lon, lat=loc2$lat,pt_graph=loc$pt_graph)
  }
  
  names(dfa) <- names(abufields)
  
 
  
  #><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
  ########################                SECTION 3                ########################
  #><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
  
  # Organize predicted abundances into suitable format for DISPLACE
  # The output will be a dataframe containing the following columns:
  # DISPLACE lon $ lat, Quarter
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3.1) Include Quarter and Size columns
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(i in seq_along(dfa)){
    dfa[[i]]$Quarter <- factor(substr(names(dfa)[i], start = 1, stop = 2)) # Double check if it's correct!!
    dfa[[i]]$Size <- paste0("SG", factor(substr(names(dfa)[i], start=5, stop=7))) # Double check if it's correct!!
  }
  
  dfa <- do.call("rbind",dfa) 
  rownames(dfa) <- NULL
  
  
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3.2) Transform to wide format
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  displace_dat <- spread(dfa, Size, abundance)
  
  #head(displace_dat)
  
  # NOTE: No data for SG0 and SG1 were present in the original dataset; The LGNB model could thus
  # not be run for these size groups and we have therfore to assume that the abundances of SG2 are
  # the same for SG0 and SG1. DISPLACE needs to have these columns informed, otherwise it will cause problems.
  
  
  # Fix:
  displace_dat$SG0 <- displace_dat$SG1
  
  # Put size columns in increasing order
  displace_dat <- displace_dat[,c("lon","lat","pt_graph","Quarter",
                                  "SG0","SG1","SG2","SG3","SG4","SG5",
                                  "SG6","SG7","SG8","SG9","SG10","SG11",
                                  "SG12","SG13")]
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3.3) Normalize abundances
  #~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Normalize function
  #~~~~~~~~~~~~~~~~~~~~
  #NormAbu <- function(x) {
  #  return ((x - min(x)) / (max(x) - min(x)))
  #}
  NormAbu <- function(x){
    return(round(x/sum(x),8))  # CAUTION: values are rounded here to save memory space....
  }
  
  
  idxna <- which(is.na(displace_dat$SG2)) #Can choose any arbitrary size-group column as the result will be the same
  displace_dat <- displace_dat[-idxna,]
  

  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3.4) Convert to DISPLACE POPULATIONS\static_avai format
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # please respect the file naming, also including "_updated" in the name.
  options(scipen = 100) # remove the scientific notation before exporting
   
  the_selected_szgroups <- c(0, 2, 3, 5, 7)
  
  annual_tstep <- as.numeric(tstep) %% 8761
  quarter <- "Q1"
  if(annual_tstep>=2160 && annual_tstep<4344) quarter <- "Q2"
  if(annual_tstep>=4344 && annual_tstep<6552) quarter <- "Q3"
  if(annual_tstep>=6552 && annual_tstep<8761) quarter <- "Q4"
  
  if(quarter=="Q1" || quarter=="Q2"){
    if(quarter=="Q1"){
      dd      <- displace_dat[displace_dat$Quarter=="Q1",]
      dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
      print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
      ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
      library(doBy)
      ddd      <- orderBy(~pt_graph, ddd)
      dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
      dat_full <- ddd[, c("pt_graph", "avai")]
    } else{
      dd      <- displace_dat[displace_dat$Quarter=="Q2",]
      dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
      print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
      ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
      library(doBy)
      ddd      <- orderBy(~pt_graph, ddd)
      dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
      dat_full <- ddd[, c("pt_graph", "avai")]
    } 
    write.table(dat, file=file.path(path, paste0("popsspe_", application), "static_avai",
                                    paste(pop, "spe_avai_szgroup_nodes_semester1_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(dat_full, file=file.path(path, paste0("popsspe_", application), "static_avai",
                                         paste(pop, "spe_full_avai_szgroup_nodes_semester1_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
  } else{
    if(quarter=="Q3"){
      dd      <- displace_dat[displace_dat$Quarter=="Q3",]
      dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
      print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
      ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
      library(doBy)
      ddd      <- orderBy(~pt_graph, ddd)
      dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
      dat_full <- ddd[, c("pt_graph", "avai")]
    } else{
      dd      <- displace_dat[displace_dat$Quarter=="Q4",]
      dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
      print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
      ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
      library(doBy)
      ddd      <- orderBy(~pt_graph, ddd)
      dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
      dat_full <- ddd[, c("pt_graph", "avai")]
    } 
    write.table(dat, file=file.path(path, paste0("popsspe_", application), "static_avai",
                                    paste(pop, "spe_avai_szgroup_nodes_semester2_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(dat_full, file=file.path(path, paste0("popsspe_", application), "static_avai",
                                         paste(pop, "spe_full_avai_szgroup_nodes_semester2_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
  
  
cat("coupling to nbcp...done\n")



} # end FALSE







  