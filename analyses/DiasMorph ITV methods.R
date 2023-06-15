library(ggplot2)
library(dynRB)
library(FD)
library(cati)
library(tidyr)
set.seed(665)

fruits <- read.delim("C:/RD/DiasMorph xFXP.txt", sep = "\t", header = TRUE)
length(unique(fruits$scientificName)) # 1437 spp

# Remove species with less than 5 observations
nxsp <- table(fruits$scientificName)
sp_20obs <- names(nxsp[nxsp>20])

sp_20obs <- fruits[fruits$scientificName %in% sp_20obs, ]
sp <- unique(sp_20obs$scientificName)

# Select 30 species randomly
rn_sp <- sample(sp, size = 30, replace = FALSE)
random_sp <- sp_20obs[sp_20obs$scientificName %in% rn_sp, ]
table(random_sp$scientificName) 

# Dynamic range boxes
DynR_list <- list()
for(i in 3:10){
  DynRB <- dynRB_VPa(A = random_sp[, c(1, 3:i)], pca.corr = T)
  DynR_list[[i]] <- DynRB$result[, c("V1", "V2", "port_mean")]
}

DynRB_df <- do.call("rbind", DynR_list)
DynRB_df$ntraits <- as.factor(sort(rep(1:8, times = nrow(DynR_list[[4]]))))
DynRB_df <- DynRB_df[DynRB_df$port_mean!=1, ] # remove trivial comparisons
DynRB_df$V1V2 <- paste(DynRB_df$V1, DynRB_df$V2, sep = "_")

ggplot(data = DynRB_df, aes(x = ntraits, y = port_mean, col = V1V2, group = V1V2)) + 
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE, alpha = 0.4) 

ggplot(data = DynRB_df, aes(x = ntraits, y = port_mean, group = ntraits)) + 
  geom_boxplot() + xlab("Number of traits") + ylab("Niche overlap")

ggplot(data = DynRB_df, aes(x = ntraits, y = port_mean, group = ntraits)) + 
  geom_jitter(position = position_jitter(0.2), alpha = 0.1) +
  geom_boxplot(alpha = 0.1, col = "darkorchid3", outlier.shape = NA) + xlab("Number of traits") + ylab("Niche overlap") +
  theme_bw()

 write.csv(DynRB_df, "dynRB_corrT_DiasMorph.csv")

 # Kernel density overlap (Geange et al. 2015)
 source("C:/MEE3_070_sm_NicheFunctions.txt")
 no.overall.mat <- list()
 for(i in 3:10){
   # Ensure the first two column names are "id" and "species".
   A.df <- data.frame(id = 1:nrow(random_sp), species = random_sp$scientificName, random_sp[, 3:i])
   # Store some vectors of names:
   spnames <- sort(unique(as.character(A.df$species)))
   no.spp <- length(spnames)
   # Make a vector of variable types to match the variable names (varnames):
   # Choose variables of interest
   varnames <- colnames(A.df)[3:i]  
   no.vars <- length(varnames) 
   vartypes <- rep("cts", no.vars) # Variables continuas
   # Set up a list of objects which are NULL if this is not a resource selection variable, 
   # and with the availability vector if it is resource selection.
   avail.list <- vector("list", no.vars)
   names(avail.list) <- varnames
   # Set up R objects to store results
   alpha.list <- vector("list", no.vars)
   names(alpha.list)<-varnames
   for (vv in 1:no.vars) if (vartypes[vv]=="rsel"){
     choices <- unique(A.df[,vv+2])
     no.ch   <- length(choices)
     alpha.list[[vv]] <- matrix(NA,no.spp,no.ch)
     dimnames(alpha.list[[vv]]) <- list(spnames,choices)
   }
   # Set up an array of niche overlaps (niche overlap function)
   no.array <- array(1, c(no.spp, no.spp, no.vars))
   dimnames(no.array) <- list(spnames,spnames,varnames)
   # Run through each variable in turn, identify its type, calculate the appropriate NO matrix 
   # and store it in the right layer of the no.array
   for (vv in 1:no.vars){
     y <- A.df[,colnames(A.df)==varnames[vv]]
     if (vartypes[vv] == "bin")
       no.array[,,vv] <- no.bin.fn(A.df$species,y)
     if (vartypes[vv] == "cat")
       no.array[,,vv] <- no.cat.fn(A.df$species,y)
     if (vartypes[vv] == "count")
       no.array[,,vv] <- no.count.fn(A.df$species,y)
     if (vartypes[vv] == "cts")
       no.array[,,vv] <- no.cts.fn(A.df$species,y)
     if (vartypes[vv] == "meas")
       no.array[,,vv] <- no.cts.fn(A.df$species,log(y))
     if (vartypes[vv] == "pcent")
       no.array[,,vv] <- no.cts.fn(A.df$species,
                                   log(y/(100-y)))
     if (vartypes[vv] == "propn")
       no.array[,,vv] <- no.cts.fn(A.df$species,
                                   log(y/(1-y)))
     if (vartypes[vv] == "rsel")
     {
       
       # Do Manly's alpha calculations, store.
       no.choices <- length(avail.list[[vv]])
       choicenames <- names(avail.list[[vv]])
       avail.vect <- avail.list[[vv]]
       alpha.mat <- alpha.fn(A.df$species,y,avail.vect)
       alpha.list[[vv]] <- alpha.mat         
       
       # Do niche overlaps, as proportions in categories:
       no.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
     }
   }
   # Also calculate overall NO measures, averaged over the dimensions
   no.overall.mat[[i]] <- apply(no.array, c(1, 2), mean)
 }
 
 # Remove first empty objects
 no.overall.mat <- Filter(length, no.overall.mat)
 
 # Extract upper elements and convert them to a vector
 to.upper.mat <- lapply(no.overall.mat, function(x) x[lower.tri(x, diag = F)])
 NO_vector <- as.vector(unlist(to.upper.mat))
 ntraits <- as.factor(sort(rep(1:8, times = length(to.upper.mat[[1]]))))
 NO_df <- data.frame(NO = NO_vector, ntraits)
 
 # 5x3.5 in
 ggplot(data = NO_df, aes(x = ntraits, y = NO, group = ntraits)) + 
   geom_jitter(position = position_jitter(0.2), alpha = 0.1) +
   geom_boxplot(alpha = 0.1, col = "darkorchid3", outlier.shape = NA) + xlab("Number of traits") + ylab("Niche overlap") +
   theme_bw()

# Euclidean distance
euc_dist_list <- list()
for(i in 3:10){
  euc_dist_list[[i]] <- dist(random_sp[, 3:i], method = "euclidean")
}

# Remove first empty objects
euc_dist_list <- Filter(length, euc_dist_list)

# Extract upper elements and convert them to a vector
to.upper.mat <- lapply(euc_dist_list, function(x) x[lower.tri(x, diag = F)])
dist_vector <- as.vector(unlist(to.upper.mat))
ntraits <- as.factor(sort(rep(1:8, times = length(to.upper.mat[[1]]))))
dist_df <- data.frame(distance = dist_vector, ntraits)

ggplot(data = dist_df, aes(x = ntraits, y = distance, group = ntraits)) + 
  geom_boxplot() + xlab("Number of traits") + ylab("Euclidean distance")

# Gower distance (8x5 in)
gow_dist_list <- list()
for(i in 3:10){
  gow_dist_list[[i]] <- gowdis(as.data.frame(random_sp[, 3:i]))
}

# Remove first empty objects
gow_dist_list <- Filter(length, gow_dist_list)

# Extract upper elements and convert them to a vector
to.upper.mat <- lapply(gow_dist_list, function(x) x[lower.tri(x, diag = F)])
dist_vector <- as.vector(unlist(to.upper.mat))
ntraits <- as.factor(sort(rep(1:8, times = length(to.upper.mat[[1]]))))
dist_df <- data.frame(distance = dist_vector, ntraits)

ggplot(data = dist_df, aes(x = ntraits, y = distance, group = ntraits)) + 
  geom_boxplot() + xlab("Number of traits") + ylab("Gower dissimilarity")
