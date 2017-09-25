#!/bin/bash
#cohorts=( GER BEL DEN SWE FIN USA AUS ITA FRA NOR UK )
cohorts=( AUSNZ )
idx=( 5921 6119 4605 1719 2526 4848 1718 5208 1406 4111 5396 2546 2653 2547 1407 5372  6120 3457 )
for co in ${cohorts[*]}
do
for id in ${idx[*]}
do
Rscript MS_full_model_sept16_SNP_analysis.R "$id" "$id"  "$co" 9 &
done
sleep 180
done
