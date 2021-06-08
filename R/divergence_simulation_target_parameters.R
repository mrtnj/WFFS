
## Set the effects of the lethal to ranges similar to the real data


library(AlphaSimR)

source("R/simulation_functions.R")


## Read command line arguments
## * selection_rule -- tells us what mating rule to use
## * lethal_is -- tells us whether to use balancing selection or not
## * n_top_exempt -- for the case of selection against the lethal, allows
##                   # of top sires to be included even if they are carriers
## * outfile for populations
## * outfile for results


args <- commandArgs(trailingOnly = TRUE)

selection_rule <- args[1]
lethal_is <- args[2]
n_top_exempt <- as.numeric(args[3])
outfile_populations <- args[4]
outfile_results <- args[5]

print(selection_rule)
print(lethal_is)
print(n_top_exempt)
print(outfile_populations)
print(outfile_results)


## Set up simulation

found_lethal <- FALSE
k <- 0

while (!found_lethal) {
    k <- k + 1
    print(k)
    setup <- make_simulation_divergence(n_ind = 6000,
                                        n_chr = 30,
                                        n_qtl = 3,
                                        n_snp = 130,
                                        h2 = 0.3,
                                        corA = 0.3)
    
    founders <- setup$founderpop
    simparam <- setup$simparam
    
    ## Pick a lethal variant; try--catch in case simulation does not have a candidate
    ## with appropriate frequency
    
    chosen_lethal <- tryCatch({
    
        pick_lethal(founders,
                    lethal_is = lethal_is,
                    simparam)    
    },
    error = function(cond) return(NULL))
    
    if (!is.null(chosen_lethal)) {
        ## Successfully found lethal; move on
        found_lethal <- TRUE   
    }
}

lethal_ix <- chosen_lethal$lethal_ix
other_snp_ix <- chosen_lethal$other_snp_ix
founder_lethal_frequency <- chosen_lethal$founder_lethal_frequency


## Set effects after the fact

modified_g1 <- simparam$traits[[1]]
modified_g2 <- simparam$traits[[2]]

modified_g1@addEff[lethal_ix] <- 5
modified_g2@addEff[lethal_ix] <- 0

simparam$switchTrait(traitPos = 1,
                     lociMap = modified_g1,
                     force = TRUE)
simparam$switchTrait(traitPos = 2,
                     lociMap = modified_g2,
                     force = TRUE)
simparam$setVarE(h2 = c(0.2, 0.2))

## Set offspring distribution

offspring_distribution <- c(1, 1, 1, 3, 4,
                            7, 13, 19, 35, 116)


## Breeding simulation

n_gen <- 20

generations <- list(n_gen)

generations[[1]] <- founders

for (gen_ix in 2:n_gen) {
    
    if (selection_rule == "avoid_carrier_x_carrier") {
    
        generations[[gen_ix]] <-
            breed_avoiding_carrier_x_carrier(generations[[gen_ix - 1]],
                                             lethal_ix,
                                             lethal_is = lethal_is,
                                             n_sires = 300,
                                             offspring_distribution = offspring_distribution,
                                             divergence = TRUE,
                                             prop_goal2 = 0.6,
                                             simparam = simparam)
        
    } else if (selection_rule == "unknown_lethal") {
        
        generations[[gen_ix]] <-
            breed_unknown_lethal(generations[[gen_ix - 1]],
                                 lethal_ix,
                                 lethal_is = lethal_is,
                                 n_sires = 300,
                                 offspring_distribution = offspring_distribution,
                                 divergence = TRUE,
                                 prop_goal2 = 0.6,
                                 simparam = simparam)
           
    } else if (selection_rule == "selection_against") {
        
        generations[[gen_ix]] <-
            breed_against_lethal(generations[[gen_ix - 1]],
                                 lethal_ix,
                                 lethal_is = lethal_is,
                                 n_sires = 300,
                                 n_top_exempt = n_top_exempt,
                                 offspring_distribution = offspring_distribution,
                                 divergence = TRUE,
                                 prop_goal2 = 0.6,
                                 simparam = simparam)
    }
}


## Create simulation object

carrier_status <- lapply(generations,
                         carrier_test,
                         lethal_ix = lethal_ix,
                         lethal_is = lethal_is,
                         simparam)

carriers <- get_carriers(carrier_status)

## Carriers at the 10% top of each goal distribution

carriers$top_goal1_frequency <- 0
carriers$top_goal2_frequency <- 0

for (gen_ix in 1:length(generations)) {
    
    rank_goal1_ix <- order(generations[[gen_ix]]@gv[,1],
                            decreasing = TRUE)
    top_goal1_ix <- rank_goal1_ix[1:(length(rank_goal1_ix) * 0.1)]
    
    rank_goal2_ix <- order(generations[[gen_ix]]@gv[,2],
                           decreasing = TRUE)
    top_goal2_ix <- rank_goal2_ix[1:(length(rank_goal2_ix) * 0.1)]
    
    carriers$top_goal1_frequency[gen_ix] <-
        sum(carrier_status[[gen_ix]][top_goal1_ix])/2/length(top_goal1_ix)
    carriers$top_goal2_frequency[gen_ix] <-
        sum(carrier_status[[gen_ix]][top_goal2_ix])/2/length(top_goal2_ix)
    
}

## Carrier sires actually used

carriers$carrier_sires_used <- 0

for (gen_ix in 1:(length(generations) - 1)) {
    
    sires_used_id <- unique(generations[[gen_ix + 1]]@father)
    sires_used_ix <- which(generations[[gen_ix]]@id %in% sires_used_id)
    
    carriers$carrier_sires_used[gen_ix] <- sum(carrier_status[[gen_ix]][sires_used_ix] == 1)

}

stats <- get_stats(generations)

if (lethal_is == "snp") {

    simulation_results <- list(simparam = simparam,
                               carrier_status = carrier_status,
                               carriers = carriers,
                               stats = stats,
                               lethal_ix = lethal_ix,
                               other_snp_ix = other_snp_ix)

} else if (lethal_is == "qtl") {

    simulation_results <- list(simparam = simparam,
                               carrier_status = carrier_status,
                               carriers = carriers,
                               stats = stats,
                               lethal_ix = lethal_ix,
                               other_snp_ix = other_snp_ix,
                               lethal_qtl_effect_goal1 = simparam$traits[[1]]@addEff[lethal_ix],
                               other_qtl_effects_goal1 = simparam$traits[[1]]@addEff[-lethal_ix],
                               lethal_qtl_effect_goal2 = simparam$traits[[2]]@addEff[lethal_ix],
                               other_qtl_effects_goal2 = simparam$traits[[2]]@addEff[-lethal_ix],
                               founder_lethal_frequency = founder_lethal_frequency,
                               founder_h2_goal1 = varG(founders)[1,1]/varP(founders)[1,1],
                               founder_h2_goal2 = varG(founders)[2,2]/varP(founders)[2,2])
    
    simulation_results$founder_variance_explained_by_lethal_goal1 <-
        2 * simulation_results$lethal_qtl_effect_goal1^2 *
        simulation_results$founder_lethal_frequency *
        (1 - simulation_results$founder_lethal_frequency) *
        simulation_results$founder_h2_goal1
    
    simulation_results$founder_variance_explained_by_lethal_goal2 <-
        2 * simulation_results$lethal_qtl_effect_goal2^2 *
        simulation_results$founder_lethal_frequency *
        (1 - simulation_results$founder_lethal_frequency) *
        simulation_results$founder_h2_goal2
}

saveRDS(generations,
        file = outfile_populations)

saveRDS(simulation_results, 
        file = outfile_results)



