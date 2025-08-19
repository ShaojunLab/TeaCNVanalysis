#!/bin/bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

R_teacnv="./TeaCNV/simulation/TeacNV_for_SimATAC_diffCNVpct_run.R"
clonal_groups=(Monoclonal Biclonal Triclonal Tetraclonal)
cnvProp=(0.3 0.4 0.6)

for p in "${cnvProp[@]}" ;do
	for group in "${clonal_groups[@]}" ;do
		for i in {1..10} ;do
			Rscript ${R_teacnv} --pct $p --ClonalType $group --nth $i --delt_lim 0.35 --DEsegCoeff 0.6 --CorrectByDist FALSE &
		done
	done
done



