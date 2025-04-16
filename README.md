# Code for: "Integrating intraspecific trait variability in functional diversity: an overview of methods and a guide for ecologists"
This repository contains data and code to reproduce analyses reported in the manuscript *Palacio FX, G Ottaviani, S Mammola, C Graco-Roza, F de Bello and CP Carmona (2025) Integrating intraspecific trait variability in functional diversity: an overview of methods and a guide for ecologists. Ecological Monographs XX: XX-XX.*

# Overview
To assess how different methods provide information on ITV applied to the same dataset, we developed a workflow to simulate intraspecific trait variability in biological assemblages (folder $tutorial$). For each species $i$ of a given assemblage $k$, we first simulated trait mean vectors **μ**<sub>i</sub> for three continuous traits following a normal distribution with arbitrary parameters $μ$ = 10 and $σ$<sup>2</sup> = 36. The great variance used had the purpose of simulating traits with rather different means, reflecting the common use of traits at different scales. Second, we simulated variance-covariance matrices **Σ**<sub>i</sub>, to account for correlations between traits, following a standard Wishart distribution with 3 degrees of freedom. With this information, we simulated trait values for a set of 10 individuals per species, following a multivariate normal distribution with parameters **μ** = **μ**<sub>i</sub> and **Σ** = **Σ**<sub>i</sub>. As a result, the mean and variance of each trait, as well as the covariances between traits, were fixed for each species in the same assemblage, and were allowed to vary across different assemblages reflecting variability in environmental conditions. We repeated this procedure for a set of 20 species and 15 assemblages, giving a trait matrix of 200 individuals by 3 traits. The same assemblage matrix of 15 sites and 20 species was used throughout, with the same overall abundance per assemblage (200 individuals). For each assemblage, abundance was modelled with a log-normal distribution with arbitrary parameters $μ$ = 5 and $σ$<sup>2</sup> = 0.25.

# Repository structure

**analyses**: contains code to reproduce analyses, simulations and figures.

**tutorial**: this walks through a tutorial for different methods describing intraspecific trait variability (ITV) and functional diversity accounting for ITV

The names of all other files should be self-explanatory. If they are not, please open an issue in this repository.
