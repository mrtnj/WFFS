
## Simple simulation with these features:
## * Discrete generations
## * 3000 females, 200 males
## * Selection only on the male side
## * One fully recessive lethal ~ 5% starting frequency
## * No carrier x carrier matings


library(AlphaSimR)


## Set up simulation

founderpop <- runMacs(nInd = 3200,
                      nChr = 10,
                      segSites = 1000)


simparam <- SimParam$new(founderpop)
simparam$restrSegSites(maxQtl = 100,
                       maxSnp = 900,
                       overlap = FALSE)
simparam$setGender("yes_sys")

simparam$addTraitA(nQtlPerChr = 100)
simparam$setVarE(h2 = 0.4)

simparam$addSnpChip(nSnpPerChr = 10)

founders <- newPop(founderpop,
                   simParam = simparam)


## Pick a lethal variant

founder_geno <- pullSnpGeno(founderpop,
                            simParam = simparam)
f <- colSums(founder_geno)/(2 * nrow(founder_geno))
candidates <- which(f > 0.04 & f < 0.06)

if (length(candidates) < 1) {
    stop("Found no candidate loci with appropriate frequency.")
}

lethal_ix <- sample(candidates, 1)


## Carrier test function

carrier_test <- function(pop, lethal_ix) {

    pullSnpGeno(pop,
                simParam = simparam)[, lethal_ix]

}



## Breeding simulation

n_gen <- 40

generations <- list(n_gen)

generations[[1]] <- founders

for (gen_ix in 2:n_gen) {
    ## Take all females as dams and split them in carrier/noncarrier

    dams <- generations[[gen_ix - 1]][generations[[gen_ix - 1]]@gender == "F"]
    
    dam_carrier_status <- carrier_test(dams,
                                       lethal_ix)

    carrier_dams <- dams[dam_carrier_status > 0]
    noncarrier_dams <- dams[dam_carrier_status == 0]

    n_carrier_dams <- sum(dam_carrier_status > 0)

    ## Select the top sires and split out noncarriers
    sires <- selectInd(pop = generations[[gen_ix - 1]],
                       trait = 1,
                       use = "pheno",
                       gender = "M",
                       nInd = 200,
                       simParam = simparam)

    sire_carrier_status <- carrier_test(sires,
                                        lethal_ix)
    
    noncarrier_sires <- sires[sire_carrier_status == 0]

    ## Create matings for noncarrier/noncarrier and carrier/noncarrier separately
    ## Splitting matings proportionally between the groups
    
    noncarrier_dam_x_any_sire <- randCross2(females = noncarrier_dams,
                                            males = sires,
                                            nProgeny = 1,
                                            nCrosses = 2 * noncarrier_dams@nInd,
                                            simParam = simparam)

    if (n_carrier_dams > 0) {
        carrier_dam_x_noncarrier_sire <- randCross2(females = carrier_dams,
                                                    males = noncarrier_sires,
                                                    nProgeny = 1,
                                                    nCrosses = 2 * carrier_dams@nInd,
                                                    simParam = simparam)
    
        generations[[gen_ix]] <- c(carrier_dam_x_noncarrier_sire,
                                   noncarrier_dam_x_any_sire)
    } else {
        generations[[gen_ix]] <- noncarrier_dam_x_any_sire
    }
}

