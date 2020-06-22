#!/bin/bash

if [ ! -d simulations ]; then
    mkdir simulations
fi

if [ ! -d simulations/simple_simulations_unknown ]; then
    mkdir simulations/simple_simulations_unknown
fi

for REP in {1..100}; do

    Rscript R/simple_simulation_unknown.R $REP
    
done