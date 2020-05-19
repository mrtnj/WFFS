#!/bin/bash

if [ ! -d simulations ]; then
    mkdir simulations
fi

if [ ! -d simulations/simple_simulations_balancing ]; then
    mkdir simulations/simple_simulations_balancing
fi

for REP in 1..100; do

    Rscript R/simple_simulation_balancing.R $REP
    
done