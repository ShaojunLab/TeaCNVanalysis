#!/bin/bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

R_teacnv="./TeaCNV/simulation/TeacNV_for_SimATAC_rareClonal_remote.R"
clonal_groups=(Biclonal Triclonal Tetraclonal Pentaclonal)
sizeSet=(50 {100..1000..100} 2000 5000 10000)


for s in "${sizeSet[@]}" ;do
	for group in "${clonal_groups[@]}" ;do
		for i in {1..10} ;do
			Rscript ${R_teacnv} --size $s --ClonalType $group --nth $i \
			--inPath /TeaCNV/simulation/data \
			--outPath /TeaCNV/simulation/rareClonal > run.log 2>&1 &
		done
	done
done
