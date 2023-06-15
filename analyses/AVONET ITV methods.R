library(ggplot2)
library(dynRB)
library(FD)
library(cati)
library(hypervolume)
library(tidyr)
set.seed(662)

avo <- read.delim("C:/RD/AVONET_birdtree_xFXP.txt", sep = "\t", header = TRUE)
avo_na <- na.omit(avo)

# Remove species with less than 5 observations
nxsp <- table(avo_na$Species)
sp_5obs <- names(nxsp[nxsp>4])

avo_5obs <- avo_na[avo_na$Species %in% sp_5obs, ]
sp <- unique(avo_5obs$Species)

# Select species from the same genus (e.g. Turdus)
genus <- avo_5obs[avo_5obs$Genus == "Turdus" | avo_5obs$Genus == "Patagioenas", ]
genus_sp <- unique(genus$Species)
length(genus_sp)
table(genus$Species)
nrow(genus)

##### Dynamic range boxes ######
DynR_list <- list()
for(i in 4:12){
DynRB <- dynRB_VPa(A = genus[, c(2, 4:i)], pca.corr = T) # 2: species column
DynR_list[[i]] <- DynRB$result[, c("V1", "V2", "port_mean")]
}

DynRB_df <- do.call("rbind", DynR_list)
DynRB_df$ntraits <- as.factor(sort(rep(1:9, times = nrow(DynR_list[[4]]))))
DynRB_df <- DynRB_df[DynRB_df$port_mean!=1, ] # remove trivial comparisons
DynRB_df$V1V2 <- paste(DynRB_df$V1, DynRB_df$V2, sep = "_")

# DynRB_df <- read.csv("C:/RD/dynRB_corrT AVONET.csv")

#ggplot(data = DynRB_df, aes(x = ntraits, y = port_mean, col = V1V2, group = V1V2)) + 
#  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE, alpha = 0.4) 

# 5x3.5 in
ggplot(data = DynRB_df, aes(x = ntraits, y = port_mean, group = ntraits)) + 
  geom_jitter(position = position_jitter(0.2), alpha = 0.1) +
  geom_boxplot(alpha = 0.1, col = "darkorchid3") + xlab("Number of traits") + ylab("Niche overlap") +
  theme_bw()

# write.csv(DynRB_df, "dynRB_corrT AVONET.csv")

######  Kernel density overlap (Geange et al. 2015) #########
source("C:/MEE3_070_sm_NicheFunctions.txt")
no.overall.mat <- list()
for(i in 4:12){
# Ensure the first two column names are "id" and "species".
A.df <- data.frame(id = 1:nrow(genus), species = genus$Species, genus[, 4:i])
# Store some vectors of names:
spnames <- sort(unique(as.character(A.df$species)))
no.spp <- length(spnames)
# Make a vector of variable types to match the variable names (varnames):
# Choose variables of interest
varnames <- colnames(A.df)[3:(i-1)]  
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
ntraits <- as.factor(sort(rep(1:9, times = length(to.upper.mat[[1]]))))
NO_df <- data.frame(NO = NO_vector, ntraits)

# 5x3.5 in
ggplot(data = NO_df, aes(x = ntraits, y = NO, group = ntraits)) + 
  geom_jitter(position = position_jitter(0.2), alpha = 0.1) +
  geom_boxplot(alpha = 0.1, col = "darkorchid3") + xlab("Number of traits") + ylab("Niche overlap") +
  theme_bw()

#### Hypervolumes #########
hv_list <- list()
for(i in 4:12){
  hv_list[[i]] <- list()
  for(j in 1:length(genus_sp)){
  sp <- subset(genus, Species == genus_sp[j])
  hv_list[[i]][j] <- hypervolume_gaussian(data = sp[, 4:i])
  }
}
# Remove first empty objects
hv_list <- Filter(length, hv_list)

nsp <- length(genus_sp)
sp_comb <- t(combn(nsp, 2))
n_comparisons <- nrow(sp_comb)
hv_overlap <- array(NA, c(n_comparisons, 4, 9))

for(k in 6:9){ # traits
  hv_join <- hypervolume_join(hv_list[[k]])
  for(i in 1:n_comparisons){ # number of species pairs
  sp1 <- sp_comb[i,1]
  sp2 <- sp_comb[i,2]
  hv_set <- hypervolume_set(hv_join[[sp1]], hv_join[[sp2]], check.memory = F)
  hv_overlap[i,,k] <- hypervolume_overlap_statistics(hv_set)
  }
}

hv_overlap_df <- apply(hv_overlap, 2, c)
hv_overlap_df$ntraits <- sort(rep(1:9, n_comparisons))
write.csv(hv_overlap_df, "overlap.csv") # hasta el trait 5 incluido

############ Distance metrics between individual observations
# Euclidean distance
euc_dist_list <- list()
for(i in 4:12){
  euc_dist_list[[i]] <- dist(genus[, 4:i], method = "euclidean")
}

# Remove first empty objects
euc_dist_list <- Filter(length, euc_dist_list)

# Extract upper elements and convert them to a vector
to.upper.mat <- lapply(euc_dist_list, function(x) x[lower.tri(x, diag = F)])
dist_vector <- as.vector(unlist(to.upper.mat))
ntraits <- as.factor(sort(rep(1:9, times = length(to.upper.mat[[1]]))))
dist_df <- data.frame(distance = dist_vector, ntraits)

ggplot(data = dist_df, aes(x = ntraits, y = distance, group = ntraits)) + 
  geom_jitter(position = position_jitter(0.2), alpha = 0.1) +
  geom_boxplot(alpha = 0.1, col = "darkorchid3", outlier.shape = NA) + xlab("Number of traits") + ylab("Euclidean distance") +
  theme_bw()

# Gower distance (8x5 in)
gow_dist_list <- list()
for(i in 4:12){
  gow_dist_list[[i]] <- gowdis(as.data.frame(genus[, 4:i]))
}

# Remove first empty objects
gow_dist_list <- Filter(length, gow_dist_list)

# Extract upper elements and convert them to a vector
to.upper.mat <- lapply(gow_dist_list, function(x) x[lower.tri(x, diag = F)])
dist_vector <- as.vector(unlist(to.upper.mat))
ntraits <- as.factor(sort(rep(1:9, times = length(to.upper.mat[[1]]))))
dist_df <- data.frame(distance = dist_vector, ntraits)

ggplot(data = dist_df, aes(x = ntraits, y = distance, group = ntraits)) + 
  geom_jitter(position = position_jitter(0.2), alpha = 0.1) +
  geom_boxplot(alpha = 0.1, col = "darkorchid3", outlier.shape = NA) + xlab("Number of traits") + ylab("Gower dissimilarity") +
  theme_bw()




# Rao entropy
n <- nrow(genus)
var_decomp <- decompCTRE(traits = genus[, 4:12], sp = genus$Species, 
                         ind.plot = c(rep("a", 200), rep("b", 16)), 
                         print = FALSE)
var_decomp
barplot(var_decomp)

nind <- 200 # Number of individuals per community
nsp <- length(genus_sp)
sites <- 10

C <- matrix(0, nrow = sites, ncol = nsp)
for(i in 1:10){
  C[i, ] <- as.numeric(sim_sad(s_pool = nsp, n_sim = nind, 
                               sad_type = "lnorm", sad_coef = list("meanlog" = 5, "sdlog" = 0.5), 
                               fix_s_sim = TRUE))
}
comm <- t(C)
rownames(comm) <- unique(genus$Species)
witRao <- c()
betRao <- c()
totRao <- c()
mean.traits <- apply(apply(genus[, 4:12], 2, scale), 2, 
                     function(x) tapply(x, genus$Species, mean, na.rm = TRUE))
for(i in 1:9){
  D <- (as.matrix(dist(mean.traits[, 1:i]))^2)/2
  witRao[i] <- mean(RaoRel(sample = comm, dfunc = D, dphyl = NULL, weight = TRUE, Jost = FALSE, structure = NULL)$FD$Alpha)
  betRao[i] <- RaoRel(sample = comm, dfunc = D, dphyl = NULL, weight = TRUE, Jost = FALSE, structure = NULL)$FD$Beta_add
  totRao[i] <- witRao[i] + betRao[i]
}

Rao_df <- data.frame(witRao, betRao)
Rao_df_perc <- data.frame(perc_witRao = 100*(witRao/totRao), 
                          perc_betRao = 100*(betRao/totRao))
Rao_long <- gather(Rao_df, key = var_source, value = var_explained, witRao:betRao, factor_key = TRUE)
Rao_long$ntraits <- as.factor(rep(1:9, 2))
Rao_long_perc <- gather(Rao_df_perc, key = var_source, value = var_explained, perc_witRao:perc_betRao, factor_key = TRUE)
Rao_long_perc$ntraits <- as.factor(rep(1:9, 2))

ggplot(data = data.frame(ntraits = as.factor(1:9), ITV = Rao_df_perc$perc_witRao), 
       aes(x = ntraits, y = ITV)) + 
  geom_bar(stat = "identity") + 
  xlab("Number of traits") + ylab("Intraspecific trait variability (%)")

ggplot(data = Rao_long, aes(x = ntraits, y = var_explained, fill = var_source)) + 
  geom_bar(position = "stack", stat = "identity") + 
  xlab("Number of traits") + ylab("Variance explained (%)")

ggplot(data = Rao_long_perc, aes(x = ntraits, y = var_explained, fill = var_source)) + 
  geom_bar(position = "stack", stat = "identity") + 
  xlab("Number of traits") + ylab("Variance explained (%)")

