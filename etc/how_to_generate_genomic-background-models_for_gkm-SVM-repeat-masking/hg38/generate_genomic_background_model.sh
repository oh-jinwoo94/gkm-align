#!/bin/sh


N=$1 # e.g., N gkm-SVM background models will be generated, which will be averaged (N=10 used for the manuscript).
n=$2 # number of DNA sequences to sample for each sample batch (n=30,000 used)
threshold=$3 # Approximate fraction of genomic sequences you intend to mask (0.1 used)


lsgkm_path="./lsgkm/"  # assuming you git cloned lsgkm to this directory. 
gkmalign_path="../../../"  


sdir="${gkmalign_path}/scripts"
genome_dir="${gkmalign_path}/data/genomes/" # running setup.sh generates this genome directory
exclude="files/hg38_DHS_ATAC_cluster_repr_2_superset.bed" # A list of all putative human regulatory regions. set to 'none' if such file not available.




python ${lsgkm_path}/scripts/nrkmers.py 11 nr11mers.fa

printf "" > to_average
for i in $(seq 1 $N)  # Generating N gkm-SVM background models using indepdent sequence samples, which will be averaged.
do

    python ${sdir}/sample_random_genomic_outside.py ${genome_dir}/hg38/ ${exclude}  random_genomic_hg38_L300_n${n}_r${i}.fa 300 ${n} ${i}
    python ${sdir}/gen_random_sequence.py random_generated_L300_n${n}_r${i}.fa 300 ${n} ${i}

    ${lsgkm_path}/src/gkmtrain -l 11 -k 7 -t 2 random_genomic_hg38_L300_n${n}_r${i}.fa  random_generated_L300_n${n}_r${i}.fa  genomic_model_hg38_L300_n${n}_r${i}


    ${lsgkm_path}/src/gkmpredict  nr11mers.fa genomic_model_hg38_L300_n${n}_r${i}.model.txt  nr11mers_r${i}.weights.out
    python ${sdir}/add_rho_kmerweights.py  nr11mers_r${i}.weights.out genomic_model_hg38_L300_n${n}_r${i}.model.txt  nr11mers_genomic_model_hg38_L300_n${n}_r${i}_weights_rho.out  # adjusts for the prediction bias. this step is not crucial for running gkm-align.

    echo  nr11mers_genomic_model_hg38_L300_n${n}_r${i}_weights_rho.out >> to_average
done


python ${sdir}/avg_weights.py to_average nr11mers_genomic_model_hg38_L300_n${n}_rAVG_weights_rho.out


${gkmalign_path}/bin/train_masker nr11mers_genomic_model_hg38_L300_n${n}_rAVG_weights_rho.out files/trainset_background_hg38-subset.txt ${genome_dir} ${threshold} human_genomic_background_model_p_${threshold}.out



