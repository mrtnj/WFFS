
## Helper functions 


## Performs carrier test for a population, returning genotypes at the
## lethal locus
##
## Parameters:
## * pop -- population object
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * simparam -- simulation parameter object

carrier_test <- function(pop,
                         lethal_ix,
                         lethal_is,
                         simparam) {
    
    if (lethal_is == "snp") {
        result <- pullSnpGeno(pop,
                              simParam = simparam)[, lethal_ix]
    } else if (lethal_is == "qtl") {
        result <- pullQtlGeno(pop,
                              simParam = simparam)[, lethal_ix]
    }
    
    result
}


## Picks a lethal variant from the founder population
##
## Parameters:
## * founderpop -- founder population object
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * simparam -- simulation parameter object

pick_lethal <- function(founderpop,
                        lethal_is,
                        simparam) {

    if (lethal_is == "snp") {
        founder_geno <- pullSnpGeno(founderpop,
                                    simParam = simparam)
    } else if (lethal_is == "qtl") {
        founder_geno <- pullQtlGeno(founderpop,
                                    simParam = simparam)
        effects <- simparam$traits[[1]]@addEff
    }
    
    f <- colSums(founder_geno)/(2 * nrow(founder_geno))
    candidates <- which(f > 0.04 & f < 0.06)
    
    if (lethal_is == "qtl") {
        candidate_effects <- effects[candidates]
        candidates <- candidates[candidate_effects > 0]
    }
    
    if (length(candidates) < 1) {
        stop("Found no candidate loci with appropriate frequency.")
    } else if (length(candidates) == 1) {
        lethal_ix <- candidates
    
    } else {
        lethal_ix <- sample(candidates, 1)
    }
    
    other_snp_ix <- setdiff(candidates, lethal_ix)
 
    list(lethal_ix = lethal_ix,
         other_snp_ix = other_snp_ix,
         founder_lethal_frequency = f[lethal_ix])
}


## Select sires either with a simple one-trait select or on two breeding
## goals after given proportion.
## Parameters:
## * parent generation
## * number of sires
## * offspring distribution expressed as offspring per sire for each decile

select_sires <- function(parent_generation,
                         n_sires,
                         offspring_distribution,
                         divergence = FALSE,
                         prop_goal2 = NULL,
                         simparam) {
    
    if (!divergence) {
        sires <- selectInd(pop = parent_generation,
                           trait = 1,
                           use = "pheno",
                           sex = "M",
                           nInd = n_sires,
                           simParam = simparam)
        
        ## Assign offspring numbers by breeding values
        offspring_per_sire <- rep(offspring_distribution,
                                  each = n_sires/10)
        offspring_per_sire <- offspring_per_sire[order(sires@pheno[,1],
                                                       decreasing = FALSE)]
    } else {
        sires_goal1 <- selectInd(pop = parent_generation,
                                 trait = 1,
                                 use = "pheno",
                                 sex = "M",
                                 nInd = n_sires * (1 - prop_goal2),
                                 simParam = simparam)
        
        offspring_per_sire_goal1 <- rep(offspring_distribution,
                                        each = n_sires * (1 - prop_goal2) / 10)
        offspring_per_sire_goal1 <- offspring_per_sire_goal1[order(sires_goal1@pheno[,1],
                                                                   decreasing = FALSE)]
        
        sires_goal2 <- selectInd(pop = parent_generation,
                                 trait = 2,
                                 use = "pheno",
                                 sex = "M",
                                 nInd = n_sires * prop_goal2,
                                 simParam = simparam)
        
        offspring_per_sire_goal2 <- rep(offspring_distribution,
                                        each = n_sires * prop_goal2 / 10)
        offspring_per_sire_goal2 <- offspring_per_sire_goal2[order(sires_goal2@pheno[,2],
                                                                   decreasing = FALSE)]
        
        sires <- c(sires_goal1,
                   sires_goal2)
        
        offspring_per_sire <- c(offspring_per_sire_goal1,
                                offspring_per_sire_goal2)
    }

    list(sires = sires,
         offspring_per_sire = offspring_per_sire)
}
    
    


## Create a founder pouplation anda simulation parameters object

make_simulation <- function(n_ind,
                            n_chr,
                            n_qtl,
                            n_snp,
                            h2) {
    
    founderpop <- runMacs(nInd = n_ind,
                          nChr = n_chr,
                          segSites = 500)
    
    simparam <- SimParam$new(founderpop)
    simparam$restrSegSites(minQtlPerChr = 100,
                           minSnpPerChr = 400,
                           overlap = FALSE)
    simparam$setSexes("yes_sys")
    
    simparam$addTraitA(nQtlPerChr = n_qtl,
                       mean = 100,
                       var = 20^2)
    simparam$setVarE(h2 = h2)
    
    simparam$addSnpChip(nSnpPerChr = n_snp)
    
    founders <- newPop(founderpop,
                       simParam = simparam)
 
    list(founderpop = founders,
         simparam = simparam)   
}

## Create a founder pouplation and a simulation parameters object
## with two correlated breeding goal traits

make_simulation_divergence <- function(n_ind,
                                       n_chr,
                                       n_qtl,
                                       n_snp,
                                       h2,
                                       corA,
                                       population_history = "generic") {
    
    if (population_history == "generic") {  
      founderpop <- runMacs(nInd = n_ind,
                            nChr = n_chr,
                            segSites = 500)
    } else if (population_history == "constant") {
      founderpop <- runMacs2(nInd = n_ind,
                            nChr = n_chr,
                            Ne = 1000,
                            histNe = NULL,
                            histGen = NULL,
                            segSites = 500)
    }
    
    simparam <- SimParam$new(founderpop)
    simparam$restrSegSites(minQtlPerChr = 100,
                           minSnpPerChr = 400,
                           overlap = FALSE)
    simparam$setSexes("yes_sys")
    
    simparam$addTraitA(nQtlPerChr = n_qtl,
                       mean = c(100, 100),
                       var = c(20^2, 20^2),
                       corA = matrix(c(1, corA,
                                       corA, 1),
                                     byrow = TRUE,
                                     ncol = 2,
                                     nrow = 2))
    simparam$setVarE(h2 = c(h2, h2))
    
    simparam$addSnpChip(nSnpPerChr = n_snp)
    
    founders <- newPop(founderpop,
                       simParam = simparam)
    
    list(founderpop = founders,
         simparam = simparam)   
}


## Create crosses between selected males and females, based on a vector of
## the number of offspring per sire, which should be based on the offspring
## distribution and the rank of the sire

make_crosses <- function(females,
                         males,
                         offspring_per_sire,
                         simParam = simparam) {
  
  n_offspring <- sum(offspring_per_sire)
  
  cross_plan <- data.frame(female_ix = numeric(n_offspring),
                           male_ix = numeric(n_offspring))
  
  ## Sample females with replacement
  cross_plan$female_ix <- sample.int(n = females@nInd,
                                     size = n_offspring,
                                     replace = TRUE)
  
  ## Repeat each sire index as many times as he needs offspring
  cross_plan$male_ix <- unlist(mapply(function(x, i) rep(x, i),
                                      1:males@nInd,
                                      offspring_per_sire,
                                      SIMPLIFY = FALSE))
  makeCross2(females,
             males,
             cross_plan,
             nProgeny = 1,
             simParam = simparam)
}



## Simulate breeding that avoids carrier--carrier matings (but allows
## carrier--noncarrier matings on both sides)
##
## Parameters:
## * parent_generation -- population object with the parents
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * n_sires -- number of sires used each generation
## * divergence -- is the simulation using two goal traits?
## * prop_goal2 -- proportion sires selected for the second goal trait
## * simparam -- simulation parameter object

breed_avoiding_carrier_x_carrier <- function(parent_generation,
                                             lethal_ix,
                                             lethal_is,
                                             n_sires,
                                             offspring_distribution,
                                             simparam,
                                             divergence = FALSE,
                                             prop_goal2 = NULL) {

    ## Take all females as dams and split them in carrier/noncarrier
    dams <- parent_generation[parent_generation@sex == "F"]

    dam_carrier_status <- carrier_test(dams,
                                       lethal_ix,
                                       lethal_is,
                                       simparam)

    carrier_dams <- dams[dam_carrier_status > 0]
    noncarrier_dams <- dams[dam_carrier_status == 0]

    n_carrier_dams <- sum(dam_carrier_status > 0)

    ## Select the top sires and split out noncarriers
    
    sire_selection <- select_sires(parent_generation,
                                   n_sires,
                                   offspring_distribution = offspring_distribution,
                                   divergence = divergence,
                                   prop_goal2 = prop_goal2,
                                   simparam = simparam)
    
    sires <- sire_selection$sires
    offspring_per_sire <- sire_selection$offspring_per_sire

    sire_carrier_status <- carrier_test(sires,
                                        lethal_ix,
                                        lethal_is,
                                        simparam)
    
    carrier_sires <- sires[sire_carrier_status > 0]
    carrier_offspring_per_sire <- offspring_per_sire[sire_carrier_status > 0]
    
    n_carrier_sires <- sum(sire_carrier_status > 0)
    
    noncarrier_sires <- sires[sire_carrier_status == 0]
    noncarrier_offspring_per_sire <- offspring_per_sire[sire_carrier_status == 0]

    ## Create matings for noncarrier/noncarrier and carrier/noncarrier separately
    ## Splitting matings proportionally between the groups
    
    noncarrier_sire_x_any_dam <- make_crosses(females = dams, 
                                              males = noncarrier_sires,
                                              offspring_per_sire = noncarrier_offspring_per_sire,
                                              simParam = simparam)
    
    if (n_carrier_sires > 0) {
      carrier_sire_x_noncarrier_dam <- make_crosses(females = noncarrier_dams, 
                                                    males = carrier_sires,
                                                    offspring_per_sire = carrier_offspring_per_sire,
                                                    simParam = simparam)
      
      offspring <- c(carrier_sire_x_noncarrier_dam,
                     noncarrier_sire_x_any_dam)
    } else {
        offspring <- noncarrier_sire_x_any_dam
    }
    
    offspring
}


## Simulate breeding with an unknown lethal. Individuals homozygous for the
## lethal are not used in breeding.
##
## Parameters:
## * parent_generation -- population object with the parents
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * n_sires -- number of sires used each generation
## * divergence -- is the simulation using two goal traits?
## * prop_goal2 -- proportion sires selected for the second goal trait
## * simparam -- simulation parameter object

breed_unknown_lethal <- function(parent_generation,
                                 lethal_ix,
                                 lethal_is,
                                 n_sires,
                                 offspring_distribution,
                                 simparam,
                                 divergence = FALSE,
                                 prop_goal2 = NULL) {
    
    ## Exclude dams who are affected
    dams <- parent_generation[parent_generation@sex == "F"]
    
    dam_carrier_status <- carrier_test(dams,
                                       lethal_ix,
                                       lethal_is,
                                       simparam)
    
    nonaffected_dams <- dams[dam_carrier_status < 2]
    
    ## Exclude sires who are affected
    potential_sires <- parent_generation[parent_generation@sex == "M"]
    
    sire_carrier_status <- carrier_test(potential_sires,
                                        lethal_ix,
                                        lethal_is,
                                        simparam)
    
    nonaffected_sires <- potential_sires[sire_carrier_status < 2]

    sire_selection <- select_sires(nonaffected_sires,
                                   n_sires,
                                   divergence = divergence,
                                   prop_goal2 = prop_goal2,
                                   offspring_distribution = offspring_distribution,
                                   simparam = simparam)

    sires <- sire_selection$sires
    offspring_per_sire <- sire_selection$offspring_per_sire
    
    ## Create matings
    
    offspring <- make_crosses(females = nonaffected_dams,
                              males = sires,
                              offspring_per_sire = offspring_per_sire,
                              simParam = simparam)
    
    offspring
}


## Get top individuals

get_top_sires <- function(population,
                          number,
                          divergence,
                          prop_goal2) {
 
    if (!divergence) {
        ranking <- order(population@pheno[, 1],
                      decreasing = TRUE)
        to_include <- ranking[1:number]
        top <- population[to_include]
        
    } else {
        ranking_goal1 <- order(population@pheno[, 1],
                               decreasing = TRUE)
        ranking_goal2 <- order(population@pheno[, 2],
                               decreasing = TRUE)
        
        to_include <- unique(c(ranking_goal1[1:number],
                               ranking_goal2[1:number]))
        top <- population[to_include]
    }
    
    top
}


## Breed with selection against carriers
##
## Parameters:
## * parent_generation -- population object with the parents
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * n_sires -- number of sires used each generation
## * n_top_exempt -- let the top # of the sires be included even
##                   if they are carriers
## * divergence -- is the simulation using two goal traits?
## * prop_goal2 -- proportion sires selected for the second goal trait
## * simparam -- simulation parameter object

breed_against_lethal <- function(parent_generation,
                                 lethal_ix,
                                 lethal_is,
                                 n_sires,
                                 offspring_distribution,
                                 simparam,
                                 n_top_exempt = 0,
                                 divergence = FALSE,
                                 prop_goal2 = NULL) {
    
    ## Exclude dams who are affected
    dams <- parent_generation[parent_generation@sex == "F"]
    
    dam_carrier_status <- carrier_test(dams,
                                       lethal_ix,
                                       lethal_is,
                                       simparam)
    
    nonaffected_dams <- dams[dam_carrier_status < 2]
    
    ## Exclude sires who are affected or carriers
    potential_sires <- parent_generation[parent_generation@sex == "M"]
    
    potential_sire_carrier_status <- carrier_test(potential_sires,
                                                  lethal_ix,
                                                  lethal_is,
                                                  simparam)
    
    noncarrier_potential_sires <- potential_sires[potential_sire_carrier_status == 0]
    
    ## Create matings
    matings_with_exempt_sires <- FALSE
    
    if (n_top_exempt > 0) {
      
      ## Find out how many exempt sires there are
      top_sires <- get_top_sires(potential_sires,
                                 number = n_top_exempt,
                                 divergence = divergence,
                                 prop_goal2 = prop_goal2)
      
      
      top_sire_carrier_status <- carrier_test(top_sires,
                                              lethal_ix,
                                              lethal_is,
                                              simparam)
      
      exempt_carrier_sires <- top_sires[top_sire_carrier_status == 1]
      
      if (exempt_carrier_sires@nInd > 0) {
        ## There are exempt carrier sires, mate them separately to
        ## noncarrier dams
        
        noncarrier_dams <- dams[dam_carrier_status == 0]
        
        sire_selection <- select_sires(c(noncarrier_potential_sires,
                                         exempt_carrier_sires),
                                       n_sires,
                                       divergence = divergence,
                                       offspring_distribution = offspring_distribution,
                                       prop_goal2 = prop_goal2,
                                       simparam = simparam)
        
        sire_carrier_status <- carrier_test(sire_selection$sires,
                                            lethal_ix,
                                            lethal_is,
                                            simparam)
        
        noncarrier_sires <- sire_selection$sires[sire_carrier_status == 0]
        noncarrier_offspring_per_sire <- sire_selection$offspring_per_sire[sire_carrier_status == 0]
        
        exempt_carrier_sires <- sire_selection$sires[sire_carrier_status == 1]
        exempt_offspring_per_sire <- sire_selection$offspring_per_sire[sire_carrier_status == 1]
        
        noncarrier_dam_x_exempt_carrier_sire <-
          make_crosses(females = noncarrier_dams,
                       males = exempt_carrier_sires,
                       offspring_per_sire = exempt_offspring_per_sire,
                       simParam = simparam)
        
        nonaffected_dam_x_noncarrier_sire <-
          make_crosses(females = nonaffected_dams,
                       males = noncarrier_sires,
                       offspring_per_sire = noncarrier_offspring_per_sire,
                       simParam = simparam)
        
        offspring <- c(noncarrier_dam_x_exempt_carrier_sire,
                       nonaffected_dam_x_noncarrier_sire)
        matings_with_exempt_sires <- TRUE
      }
    } 
    
    if (!matings_with_exempt_sires) {
      sire_selection <- select_sires(noncarrier_potential_sires,
                                     n_sires,
                                     divergence = divergence,
                                     offspring_distribution = offspring_distribution,
                                     prop_goal2 = prop_goal2,
                                     simparam = simparam)
      
      noncarrier_sires <- sire_selection$sires
      noncarrier_offspring_per_sire <- sire_selection$offspring_per_sire
      
      offspring <- make_crosses(females = nonaffected_dams,
                                males = noncarrier_sires,
                                offspring_per_sire = noncarrier_offspring_per_sire,
                                simParam = simparam)
    }
    
    offspring
}





###########################

## Functions for analysing simulation results

## Read simulation results

read_results <- function(filenames) {

    n_reps <- length(filenames)
    simulation_results <- lapply(filenames,
                                 readRDS)
    names(simulation_results) <- 1:n_reps

    simulation_results
}


## Get genetic trends

get_stats <- function(generations) {
    
    mean_g <- lapply(generations, meanG)
    var_g <- lapply(generations, varG)
    
    if (length(mean_g[[1]]) == 1) {
        ## One breeding goal trait
        stats <- data.frame(generation = 1:length(generations),
                             mean_g = unlist(mean_g),
                             var_g = unlist(var_g),
                             stringsAsFactors = FALSE)
        
    } else if (length(mean_g[[1]]) == 2) {
        mean_g1 <- unlist(lapply(mean_g, "[", 1))
        mean_g2 <- unlist(lapply(mean_g, "[", 2))
        
        var_g1 <- unlist(lapply(var_g, "[", 1))
        var_g2 <- unlist(lapply(var_g, "[", 4))
        cov_g <- unlist(lapply(var_g, "[", 2))
        
        cor_g <- cov_g / (sqrt(var_g1) * sqrt(var_g2))
        
        stats <- data.frame(generation = 1:length(generations),
                            mean_g1 = mean_g1,
                            mean_g2 = mean_g2,
                            var_g1 = var_g1,
                            var_g2 = var_g2,
                            cor_g = cor_g,
                            stringsAsFactors = FALSE)
    }
    
    stats
}

## Get combined stats from a list of simulations

combined_stats <- function(simulation_results) {
    n_reps <- length(simulation_results)
    n_generations <- nrow(simulation_results[[1]]$stats)
    
    stats <- map_df(simulation_results, function(x) x$stats)
    stats$replicate <- rep(1:n_reps, each = n_generations)
    
    stats
}


## Improvement per year expressed in genetic sd

get_sd_per_year <- function(rep_stats) {
    rep_stats$mean_g[-1] - rep_stats$mean_g[-length(rep_stats$mean_g)]
}



## Get carrier numbers for lethal

get_carriers <- function(carrier_status) {
    carriers <- unlist(lapply(carrier_status,
                              function (x) sum(x > 0)))
    cases <- unlist(lapply(carrier_status,
                           function (x) sum(x == 2)))
    n <- unlist(lapply(carrier_status,
                       length))
    data.frame(generation = 1:length(carrier_status),
               carriers = carriers,
               cases = cases,
               n = n)
}

## Get combined carrier numbers for list of simulations

combined_carriers <- function(simulation_results) {
    n_reps <- length(simulation_results)
    n_generations <- nrow(simulation_results[[1]]$stats)
    
    # carriers <- map_df(simulation_results, function(x) get_carriers(x$carrier_status))
    carriers <- map_df(simulation_results, function(x) x$carriers)
    carriers$replicate <- rep(1:n_reps, each = n_generations)

    carriers    
}

## Summarise carriers during last generations from a set of simulations

get_end_carriers <- function(carriers,
                             gen_start = 16,
                             gen_end = 20) {
    do(group_by(carriers, replicate),
       data.frame(average_carriers = mean(.$carriers[.$generation %in% gen_start:gen_end]),
                  average_n = mean(.$n[.$generation %in% gen_start:gen_end])))
}


## Summary statistics for end carriers from a set of simulations

end_carrier_stats <- function(end_carriers) {
    
    f <- end_carriers$average_carriers/2/end_carriers$average_n
    
    data.frame(mean_frequency = mean(f),
               median_frequency = median(f),
               lower = quantile(f, 0.05),
               upper = quantile(f, 0.95))
    
}


## Summary of QTL frequency and lethal frequency for balancing selection cases

get_balancing_stats <- function(simulation_results) {

    founder_variance_explained <- unlist(lapply(simulation_results,
                                                "[",
                                                "founder_variance_explained_by_lethal"))

    end_carriers <- get_end_carriers(combined_carriers(simulation_results))
    

    data.frame(founder_variance_explained = founder_variance_explained,
               end_carrier_frequency = end_carriers$average_carriers/2/end_carriers$average_n)
}


get_balancing_stats_divergence <- function(simulation_results) {
    
    data.frame(founder_variance_explained_g1 = unlist(lapply(simulation_results,
                                                             "[",
                                                             "founder_variance_explained_by_lethal_goal1")),
               founder_variance_explained_g2 = unlist(lapply(simulation_results,
                                                             "[",
                                                             "founder_variance_explained_by_lethal_goal2")),
               lethal_qtl_effect_g1 = unlist(lapply(simulation_results,
                                                    "[",
                                                    "lethal_qtl_effect_goal1")),
               lethal_qtl_effect_g2 = unlist(lapply(simulation_results,
                                                    "[",
                                                    "lethal_qtl_effect_goal2")),
               top_carrier_frequency_g1 = unlist(lapply(simulation_results,
                                                        function(x) x$carriers$top_goal1_frequency[20])),
               top_carrier_frequency_g2 = unlist(lapply(simulation_results,
                                                        function(x) x$carriers$top_goal2_frequency[20])))
}



## Genetic trend plot from stats

genetic_trend_plot <- function(stats) {
    qplot(x = generation,
          y = mean_g,
          colour = replicate,
          group = replicate,
          data = stats,
          geom = "line",
          xlab = "Generation",
          ylab = "Mean genetic value",
          main = "Genetic trend")
}


genetic_variance_plot <- function(stats) {
    qplot(x = generation,
          y = var_g,
          data = stats,
          colour = replicate,
          group = replicate,
          geom = "line",
          xlab = "Generation",
          ylab = "Genetic variance",
          main = "Genetic variance trend")
}

lethal_frequency_plot <- function(carriers) {
    qplot(x = generation,
          y = carriers/2/n,
          data = carriers,
          colour = replicate,
          group = replicate,
          xlab = "Generation",
          ylab = "Lethal allele frequency",
          geom = "line")
}


balancing_stat_plot <- function(balancing_stats) {
    qplot(x = founder_variance_explained,
          y = end_carrier_frequency,
          data = balancing_stats) +
        xlab("Variance explained in breeding goal in founder population") +
        ylab("Lethal allele frequency")
}


affected_plot <- function(carriers) {
    qplot(x = generation,
          y = cases,
          data = carriers,
          colour = replicate,
          group = replicate,
          geom = "line")
}
