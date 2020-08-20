#!/bin/bash

if [ ! -d simulations ]; then
    mkdir simulations
fi

if [ ! -d simulations/divergence_simulations ]; then
    mkdir simulations/divergence_simulations
fi

for SELECTION_RULE in avoid_carrier_x_carrier unknown_lethal; do

    for LETHAL_IS in snp qtl; do
    
        if [ ! -d simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS} ]; then
            mkdir simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}
        fi

        for REP in {1..100}; do

            Rscript R/divergence_simulation.R \
                $SELECTION_RULE \
                $LETHAL_IS \
                0 \
                simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}/populations_$REP.Rds \
                simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}/results_$REP.Rds 
                
        done
        
    done    
    
done


SELECTION_RULE=selection_against

for LETHAL_IS in snp qtl; do
    
    if [ ! -d simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}_${N_TOP_EXEMPT}exempt ]; then
        mkdir simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}_${N_TOP_EXEMPT}exempt
    fi

    for REP in {1..100}; do
    
        for N_TOP_EXEMPT in 0 10 100; do

            Rscript R/divergence_simulation.R \
                $SELECTION_RULE \
                $LETHAL_IS \
                $N_TOP_EXEMPT \
                simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}_${N_TOP_EXEMPT}exempt/populations_$REP.Rds \
                simulations/divergence_simulations/${SELECTION_RULE}_${LETHAL_IS}_${N_TOP_EXEMPT}exempt/results_$REP.Rds 
            
    done
    
done