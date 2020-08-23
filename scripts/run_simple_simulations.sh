#!/bin/bash

if [ ! -d simulations ]; then
    mkdir simulations
fi

if [ ! -d simulations/simple_simulations ]; then
    mkdir simulations/simple_simulations
fi

for SELECTION_RULE in avoid_carrier_x_carrier unknown_lethal; do

    for LETHAL_IS in snp qtl; do
    
        if [ ! -d simulations/simple_simulations/${SELECTION_RULE}_${LETHAL_IS} ]; then
            mkdir simulations/simple_simulations/${SELECTION_RULE}_${LETHAL_IS}
        fi

        for REP in {1..100}; do

            Rscript R/simple_simulation.R \
                $SELECTION_RULE \
                $LETHAL_IS \
                0 \
                simulations/simple_simulations/${SELECTION_RULE}_${LETHAL_IS}/populations_$REP.Rds \
                simulations/simple_simulations/${SELECTION_RULE}_${LETHAL_IS}/results_$REP.Rds 
                
        done
        
    done    
    
done

