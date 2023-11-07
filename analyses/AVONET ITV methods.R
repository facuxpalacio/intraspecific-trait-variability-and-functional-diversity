library(ggplot2)
library(dynRB)
library(FD)
library(cati)
library(tidyr)
set.seed(662)

avo <- read.delim("C:/RD/AVONET_birdtree_xFXP.txt")
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

######  Effects of scale on ITV #########
# Species widely distributed
nxcountry <- table(avo_na$Species, avo_na$Country)
nxcountry_pa <- 1*(nxcountry>0)
ncountries <- rowSums(nxcountry_pa)
sp <- names(ncountries[ncountries>5])
avo_5countries <- avo_na %>% filter(., Species %in% sp)
nsp <- length(sp)

CV_beak <- list()

# Loop through each species
for (j in 1:nsp) {
  spj <- subset(avo, Species == sp[j])
  country <- unique(spj$Country)
  ncountry <- length(country)
  
  # Initialize a vector to store CV values for this species
  CV_species <- numeric(ncountry)
  
  # Loop through each country
  for (i in 1:ncountry) {
    sp_country <- spj[spj$Country %in% country[1:i], ]
    
    # Calculate coefficient of variation (CV)
    cv_value <- 100 * sd(sp_country$Beak.Length_Culmen, na.rm = TRUE) /
      mean(sp_country$Beak.Length_Culmen, na.rm = TRUE)
    
    # Store the CV value for this country
    CV_species[i] <- cv_value
  }
  
  # Store the CV values for this species in the list
  CV_beak[[j]] <- CV_species
}

# Initialize vectors to store concatenated CV values and corresponding x-values
all_CV_values <- numeric(0)
all_x_values <- numeric(0)

# Create a vector to store species labels
species_labels <- character(0)

# Loop through each species
for (j in 1:length(CV_beak)) {
  CV_values <- CV_beak[[j]]
  num_elements <- length(CV_values)
  
  # Create x-values for the current species
  x_values <- 1:num_elements
  
  # Concatenate CV values and x-values to the overall vectors
  all_CV_values <- c(all_CV_values, CV_values)
  all_x_values <- c(all_x_values, x_values)
  
  # Create a vector of species labels for each data point
  species_labels <- c(species_labels, rep(j, num_elements))
}

# Create a data frame for ggplot2
data <- data.frame(Index = all_x_values, CV = all_CV_values, Species = species_labels)

# Plot CV vs geographic range
ggplot(data, aes(x = Index, y = CV, color = factor(Species))) +
  geom_point(size = 3, alpha = 0.4) + 
  geom_smooth(size = 1.2, method = "loess", se = FALSE, span = 0.8) +  
  xlab("Number of countries a species occurs in") +
  ylab("Coefficient of variation") +
  scale_color_discrete(name = "Species", labels = sp) +  
  theme_minimal()  


