#### Load packages
library(ggplot2) # Plotting
library(GGally) # Plotting
library(fundiv) # Dendrogram-based index
library(geometry) # for TOP index (convex hulls)
library(geozoo) # for TED index (spehere random points)
library(flexmix) # for TED index (KL divergence)
library(TPD) # Trait-probability density
library(alphahull) # Convex hulls
library(BAT) # Hypervolumes
library(lme4) # Mixed models

#### TOP and TED index functions
#### Function to compute TOP index (from Fontana et al. 2015)
TOP.index <- function(traitdat){
  # TOP
  dim1 <- ncol(traitdat)
  #definitions: index i, area as empty vector
  i=0
  area<-matrix(ncol=2,nrow=nrow(traitdat))
  while(nrow(traitdat)>dim1){
    i=i+1
    # use of convhulln function
    # area
    area[i,2] <- convhulln(traitdat,"FA")$area
    # identity of vertices
    vert0<-convhulln(traitdat,"Fx TO 'vert.txt'")
    vert1<-scan("vert.txt",quiet=T)
    vert2<-vert1+1
    
    vertices <- vert2[-1]
    
    traitdat <- traitdat[-vertices,]
    
    area[i,1] <- length(vertices)
    
  }
  area<-na.omit(area)
  # Output (2 numbers): Number of points touched by areas; Sum of the areas (TOP index)
  colSums(area)
}

#### Function to compute TED index (from Fontana et al. 2015)
# List of possible REFERENCES
# Define maximum number of points (max1) and number of traits (dim1)
max1 <- 200
dim1 <- 2 

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)

# TED index calculation
TED.index <- function(traitdat){
# Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  n.sample<-nrow(traitdat)
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
# Transform REFERENCE in data frame
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
# Ev. delete individuals in order to have the same number as in the sample
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  # Compare with sample
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

#### Set simulation parameters
set.seed(665)
# Number of simulations
nsim <- 10

# Number of communities
ncomms <- 10

# Species richness
nsp <- 10 

# Number of individuals sampled per species
nindxsp <- 10

# Community means and variances (the greater the overall variance and the lower the community variance, the greater the external filtering)
min_trait1_comm <- seq(20, 10, length = nsim)
max_trait1_comm <- seq(25, 50, length = nsim)

mean1_comm <- as.data.frame(matrix(NA, nrow = nsim, ncol = ncomms))
for(i in 1:nsim){
mean1_comm[i, ] <- seq(min_trait1_comm[i], max_trait1_comm[i], length = ncomms)
}

mean2_comm <- as.data.frame(matrix(NA, nrow = nsim, ncol = ncomms))
for(i in 1:nsim){
sd_comms <- seq(5, 20, length = nsim)
mean2_comm[i, ] <- rnorm(n = ncomms, mean = 50, sd = sd_comms)
}
colnames(mean1_comm) <- colnames(mean2_comm) <- paste("comm", 1:ncomms, sep = "")

# Overall trait means (fixed for every iteration)
mean_trait1 <- rowMeans(mean1_comm)
mean_trait2 <- rowMeans(mean2_comm)

# Create combinations of each community mean with each community variance (the lower the variance, the greater the external filtering)
mean1_comm_sim <- mean1_comm[rep(1:nsim, nsim), ]
mean1_comm_sim$cv_comm <- sort(rep(seq(0.3, 0.6, length = nsim), nsim)) ##### Parameter to change
SD1comm <- mean1_comm_sim$cv_comm*mean1_comm_sim[, 1:ncomms]

mean2_comm_sim <- mean2_comm[rep(1:nsim, nsim), ]
mean2_comm_sim$cv_comm <- sort(rep(seq(0.3, 0.6, length = nsim), nsim)) ##### Parameter to change
SD2comm <- mean2_comm_sim$cv_comm*mean2_comm_sim[, 1:ncomms]

# Intraspecific variance (the lower the variance, the greater the internal filtering)
cv_sp <- seq(0.05, 0.4, length = nsim) ##### Parameter to change

#### Trait matrix at the individual level
trait_matrix <- list()

# Matrix of species mean traits (community by trait matrix)
N <- nrow(mean1_comm_sim)
mean_trait1_sp <- mean_trait2_sp <- matrix(NA, nrow = N, ncol = nsp*ncomms)

for(i in 1:N){
# Species mean traits per community (community by species matrix)
mean_trait1_sp[i, ] <- rnorm(n = nsp*ncomms, mean = sort(rep(as.numeric(mean1_comm_sim[i, 1:ncomms]), nindxsp)), 
                        sd = as.numeric(SD1comm[i, ])) + 40
mean_trait2_sp[i, ] <- rnorm(n = nsp*ncomms, mean = sort(rep(as.numeric(mean2_comm_sim[i, 1:ncomms]), nindxsp)), 
                        sd = as.numeric(SD2comm[i, ])) + 75
}

# Individual trait values per species 
N <- nsim*nsim
list_names <- 1:N
list_trait1_matrix <- list_trait2_matrix <- list()

for(i in 1:N){
list_trait1_matrix[[i]] <- list_trait2_matrix[[i]] <- list()
  for(j in 1:nsim){
ITV1 <- mean_trait1_sp[i, ]*cv_sp[j] # vector of intraspecific trait variances 
list_trait1_matrix[[i]][[j]] <- rnorm(n = nindxsp*nsp*ncomms, mean = sort(rep(mean_trait1_sp[i, ], nindxsp)), 
                                                       sd = ITV1)
ITV2 <- mean_trait2_sp[i, ]*cv_sp[j] # vector of intraspecific trait variances 
list_trait2_matrix[[i]][[j]] <- rnorm(n = nindxsp*nsp*ncomms, mean = sort(rep(mean_trait2_sp[i, ], nindxsp)), 
                                                       sd = ITV2)

# simulation parameters as names
names(list_trait1_matrix)[i] <- names(list_trait2_matrix)[i] <- paste("mincomm", round(mean1_comm_sim[i, 1], 1), "maxcomm",  round(mean1_comm_sim[i, ncol(mean1_comm_sim)-1], 1), "CVcomm", round(mean1_comm_sim$cv_comm[i], 1), sep = "_")
names(list_trait1_matrix[[i]])[j] <- names(list_trait2_matrix[[i]])[j] <- paste("CVintrasp", round(cv_sp[j], 1), sep = "_")
}
}

trait1_matrix <- as.data.frame(do.call("cbind", do.call("cbind", list_trait1_matrix)))
trait2_matrix <- as.data.frame(do.call("cbind", do.call("cbind", list_trait2_matrix)))
colnames(trait1_matrix) <- names(rapply(list_trait1_matrix, function(x) head(x, 1)))
colnames(trait2_matrix) <- names(rapply(list_trait2_matrix, function(x) head(x, 1)))

n_sim_comms <- nrow(trait1_matrix)*ncol(trait1_matrix) # Total number of simulated communities

comm <- sort(rep(1:ncomms, nindxsp*nsp))
sp <- rep(sort(rep(1:nsp, nindxsp)), ncomms)

# Get simulation scenarios names
sim_names_list <- strsplit(names(trait1_matrix), "_")
range_trait1 <- c()
CVcomm <- c()
CVintrasp <- c()
for(i in 1:length(sim_names_list)){
  range_trait1[i] <- as.numeric(sim_names_list[[i]][4]) - as.numeric(sim_names_list[[i]][2])
  CVcomm[i] <- as.numeric(substring(sim_names_list[[i]][6], 1, 3))
  CVintrasp[i] <-  as.numeric(sim_names_list[[i]][7])
}


#### Plot community distribution function
# plot_comm_trait function for plotting intraspecific trait variability
# trait_matrix = dataframe with a trait, species and community labels
# trait = vector of trait values
# comm_label = vector of community labels
# sp_label = vector of species labels
# species = whether to plot species trait distributions

plot_comm_trait <- function(trait_matrix, trait, comm_label, sp_label, species = FALSE){
plot(x = trait, y = trait, xlab = "Trait", ylab  = "Density",
     xlim = c(20, 100),
     ylim = c(0, 1), type = "n")

# Plot community distributions
col_comm <- sort(rainbow(length(unique(comm_label))))
ncomms <- length(unique(comm_label)) # number of communities
for(i in 1:ncomms){
comm_i <- subset(trait_matrix, comm_label == i)
x <- seq(min(comm_i$trait, na.rm = T) - 3, max(comm_i$trait, na.rm = T) + 3, length = 200)
y <- dnorm(x, mean = mean(comm_i$trait, na.rm = T), sd = sd(comm_i$trait, na.rm = T))
lines(x, y, lwd = 3, col = col_comm[i])
}

if (species == TRUE){
# Plot species distributions
col_sp <- sort(rep(rainbow(length(unique(comm_label))), nsp)) # same color for each spp belonging to the same community
nsp <- length(unique(sp_label)) # number of species
N <- nsp*ncomms
sp_comm <- paste(sp_label, comm_label, sep = ".")
nindxsp <- as.numeric(table(sp_comm)) # number of individuals per species and community
trait_matrix$sp_dist <- sort(rep(1:N, nindxsp))
for(i in 1:N){
sp_i <- subset(trait_matrix, sp_dist == i)
x <- seq(min(sp_i$trait, na.rm = T) - 3, max(sp_i$trait, na.rm = T) + 3, length = 200)
y <- dnorm(x, mean = mean(sp_i$trait, na.rm = T), sd = sd(sp_i$trait, na.rm = T))
lines(x, y, col = col_sp[i])
}
}
}

# Check ecological consistence
# High external filtering, high internal filtering
subtrait <- data.frame(trait = trait1_matrix$mincomm_10_maxcomm_50_CVcomm_0.3.CVintrasp_0, comm = comm, sp = sp)
plot_comm_trait(trait_matrix = subtrait, trait = subtrait$trait, comm_label = subtrait$comm, sp_label = subtrait$sp, species = F)

# High external filtering, low internal filtering
subtrait <- data.frame(trait = trait1_matrix$mincomm_10_maxcomm_50_CVcomm_0.3.CVintrasp_0.2, comm = comm, sp = sp)
plot_comm_trait(trait_matrix = subtrait, trait = subtrait$trait, comm_label = subtrait$comm, sp_label = subtrait$sp, species = F)

# Low external filtering, high internal filtering
subtrait <- data.frame(trait = trait1_matrix$mincomm_20_maxcomm_25_CVcomm_0.4.CVintrasp_0, comm = comm, sp = sp)
plot_comm_trait(trait_matrix = subtrait, trait = subtrait$trait, comm_label = subtrait$comm, sp_label = subtrait$sp, species = F)

# Low external filtering, low internal filtering
subtrait <- data.frame(trait = trait1_matrix$mincomm_20_maxcomm_25_CVcomm_0.3.CVintrasp_0.2, comm = comm, sp = sp)
plot_comm_trait(trait_matrix = subtrait, trait = subtrait$trait, comm_label = subtrait$comm, sp_label = subtrait$sp, species = F)

#### Metric calculations
# Correlations between traits among communities
cor_comms <- NULL
for(j in 1:ncol(trait1_matrix)){
  cor_comms[j] <- cor(trait1_matrix[, j], trait2_matrix[, j], use = "complete.obs")
}

# Extract 1000 simulated communities (100 columns)
one_every_10_column <- rep(c(1, rep(0, 9)), 10)
trait1a_matrix <- trait1_matrix[, one_every_10_column == 1]
trait2a_matrix <- trait2_matrix[, one_every_10_column == 1]

# Check the minimum number of individuals in each simulation (for TPD it must be >4)
min_obs_sim <- c()
for(j in 1:ncol(trait1a_matrix)){ 
   traits <- na.omit(data.frame(trait1 = trait1a_matrix[, j], trait2 = trait2a_matrix[, j], comm, sp))
   min_obs_sim[j] <- min(table(traits$comm,traits$sp))
}
min_obs_sim[min_obs_sim<5]
which(min_obs_sim %in% min_obs_sim[min_obs_sim<5])

# Compute FD metrics
dendroFD <- TOP_comms <- TED_comms <- MVNH_det_comms <- TPD_FRic <- TPD_FEve <- TPD_FDiv <- hv_richness <- hv_regularity <- hv_divergence <- matrix(NA, nrow = ncomms, ncol = ncol(trait1a_matrix))
Ck <- NULL
comm_names <- unique(comm)

for(j in 1:ncol(trait1a_matrix)){ 
  for(i in 1:ncomms){
   traits <- na.omit(data.frame(trait1 = trait1a_matrix[, j], trait2 = trait2a_matrix[, j], comm, sp))
   subcom <- subset(traits, comm == comm_names[i])
   subtrait_matrix <- as.data.frame(scale(subcom[, c("trait1", "trait2")]))
    TOP_comms[i, j] <- TOP.index(subcom[, c("trait1", "trait2")])[2]
    TED_comms[i, j] <- TED.index(subcom[, c("trait1", "trait2")])
    MVNH_det_comms[i, j] <- det(cov(subcom[, c("trait1", "trait2")]))
  }

  # Community by species matrix
   C <- as.matrix(table(traits$comm, traits$sp))
  
  # Community by individual matrix
   traits$sp_ind <- 1:nrow(traits)
   Cind <- as.data.frame.matrix(table(traits$comm, traits$sp_ind))
   colnames(Cind) <- 1:nrow(traits)
   rownames(traits) <- 1:nrow(traits)

  # Dendrogram-based FD
   dendroFD[, j] <- FD_dendro(S = traits[, c("trait1", "trait2")], A = Cind, w = NA, Distance.method = "gower", ord = "podani", Cluster.method = "average", stand.x = TRUE, Weigthedby = "abundance")$FDpg

  # TPDs and TPDc
   TPDs_spp <- TPDs(species = traits$sp, traits = traits[, c("trait1", "trait2")], samples = traits$comm) 
   TPDc_comm <- TPDc(TPDs = TPDs_spp, sampUnit = C)
   #plotTPD(TPDs_spp, nRowCol = c(5,2))
   #plotTPD(TPDc_comm, nRowCol = c(5,2))

  # TPD_FD
   TPD_FD <- REND(TPDc = TPDc_comm)
   TPD_FRic[, j] <- TPD_FD$communities$FRichness
   TPD_FEve[, j] <- TPD_FD$communities$FEvenness
   TPD_FDiv[, j] <- TPD_FD$communities$FDivergence
 
  # Hypervolumes
   hvlist <- kernel.build(comm = Cind, trait = traits[, c("trait1", "trait2")], axes = 0, distance = "euclidean", method = "gaussian", abund = FALSE, samples.per.point = 5)

  # Compute functional diversity metrics
   hv_richness[, j] <- kernel.alpha(hvlist)
   hv_regularity[, j] <- kernel.evenness(hvlist)
   hv_divergence[, j] <- kernel.dispersion(hvlist, func = "divergence")
 }

# Simulation parameters retained
strings_split <- strsplit(names(trait1_matrix),"_")
range_trait1_sim <- c()
CVcomm_sim <- c()
CVintrasp_sim <- c()
for(i in 1:length(strings_split)){
range_trait1_sim[i] <- as.numeric(strings_split[[i]][4]) - as.numeric(strings_split[[i]][2])
CVcomm_sim[i] <- as.numeric(substr(strings_split[[i]][6], 1, 3))
CVintrasp_sim[i] <- as.numeric(strings_split[[i]][7])
}

FD_itv <- data.frame(dendroFD = as.vector(dendroFD), TOP = as.vector(TOP_comms), TED = as.vector(TED_comms), 
  MVNHdet = as.vector(MVNH_det_comms),
  TPD_FRich = as.vector(TPD_FRic),
  TPD_FEve = as.vector(TPD_FEve), 
  TPD_FDiv = as.vector(TPD_FDiv),
  HV_Rich = as.vector(hv_richness),
  HV_Reg = as.vector(hv_regularity),
  HV_Dev = as.vector(hv_divergence), 
  n_sim = sort(rep(1:100, nsim)), 
  range_trait1 = rep(range_trait1_sim, each = nsim), 
  CVcomm = rep(CVcomm_sim, each = nsim), 
  CVintrasp = rep(CVintrasp_sim, each = nsim))

#write.csv(FD_itv, "FD_itv_sims.csv")
#FD_itv <- read.table("FD_itv_sims.txt", header = TRUE)

# Metrics correlations (15x8)
p <- ggpairs(FD_itv[, 1:10], upper = list(continuous = wrap("cor", method = "spearman"))) +
ggplot2::theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))
p

vars <- c(1:10, 1:10)
vals <- unique(vars)
every_comb <- t(combn(vals, 2)) # every possible combination of variables

FD_itv$n_sim <- as.factor(FD_itv$n_sim)
z.coeff <- c()
ci <- matrix(NA, nrow = length(every_comb), ncol = 2)
predictor <- c()
response <- c()
for(i in 1:nrow(every_comb)){
p <- scale(FD_itv[, every_comb[i, 1]])
r <- scale(FD_itv[, every_comb[i, 2]])
df <- data.frame(r, p, FD_itv$n_sim)
model <- lmer(r ~ p + (1|n_sim), data = FD_itv)
z.coeff[i] <- summary(model)$coeff[2,1]
ci[i,] <- confint(model)[4,]
predictor[i] <- colnames(FD_itv)[every_comb[i, 1]]
response[i] <- colnames(FD_itv)[every_comb[i, 2]]
}

coeff_df <- data.frame(predictor, response, z.coeff, ci)

#write.csv(coeff_df, "FD_itv_sims_mixed_models.csv")
