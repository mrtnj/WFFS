
## Simple simulation with these features:
## * Discrete generations
## * 3000 females, 200 males
## * Selection only on the male side
## * One fully recessive lethal ~ 5% starting frequency
## * No carrier x carrier matings


library(AlphaSimR)


source("R/simulation_functions.R")



for (rep_ix in 1:10) {

    ## Set up simulation

    setup <- make_simulation(n_ind = 3200,
                             n_chr = 10,
                             h2 = 0.4)

    founders <- setup$founderpop
    simparam <- setup$simparam


    ## Pick a lethal variant

    chosen_lethal <- pick_lethal(founderpop,
                             lethal_is = "snp",
                             simparam)
    lethal_ix <- chosen_lethal$lethal_ix
    other_snp_ix <- chosen_lethal$other_snp_ix


    ## Breeding simulation

    n_gen <- 40

    generations <- list(n_gen)

    generations[[1]] <- founders

    for (gen_ix in 2:n_gen) {

        generations[[gen_ix]] <-
        breed_avoiding_carrier_x_carrier(generations[[gen_ix - 1]],
                                         lethal_ix,
                                         lethal_is = "snp",
                                         simparam)
    }
    

    ## Create simulation object

    simulation_results <- list(generations = generations,
                               carrier_status = lapply(generations,
                                                       carrier_test,
                                                       lethal_ix = lethal_ix,
                                                       lethal_is = "snp",
                                                       simparam),
                               simparam = simparam,
                               lethal_ix = lethal_ix,
                               other_snp_ix = other_snp_ix)
    
    saveRDS(simulation_results, 
            file = paste("simulations/simple_simulations/results_",
                         rep_ix,
                         ".Rds",
                         sep = ""))
}
            
    
    