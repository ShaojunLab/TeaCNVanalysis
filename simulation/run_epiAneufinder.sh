#!/bin/bash

R_teacnv="./TeaCNV/simulation/epiAneufinder_diffCNVpct_run.R"
clonal_groups=(Monoclonal Biclonal Triclonal Tetraclonal)
cnvProp=(0.3 0.4 0.6)

for p in "${cnvProp[@]}" ;do
	for group in "${clonal_groups[@]}" ;do
		for i in {1..10} ;do
			Rscript ${R_teacnv} --pct $p --ClonalType $group --nth $i > teacnv.log 2>&1 &
		done
	done
done




