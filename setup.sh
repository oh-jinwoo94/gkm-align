#!/bin/sh


# compile gkm-align
cd src
make clean
make
cd ..



# download genome background models 
mkdir data
printf "\nhg38, mm10 genome background models downloaded to data/. \n"
wget -q https://beerlab.org/gkmalign/masker_model_hg38_outside_union_repr_kmer_cluster_threshold_pieces_p_0.1.out -O data/human_genomic_background_model_p_0.1.out 
wget -q https://beerlab.org/gkmalign/masker_model_mm10_outside_union_repr_kmer_cluster_threshold_pieces_p_0.1.out -O data/mouse_genomic_background_model_p_0.1.out



# download genomes 
while true; do
	printf "\nDownload hg38 and mm10 genomes? (y/n) \n"
	read -p "y recommended for the github tutorial.   " yn
	case $yn in
	[Yy]* )
	 cd data
	 mkdir genomes
	 cd genomes
	 mkdir hg38
	 mkdir mm10

	 cd hg38
	 printf "\nDownloading hg38 to data/genomes/.\n"
	 for c in {1..22} X Y M 
	 do
	     printf "chr${c}\n"
	     wget -q  "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${c}.fa.gz"
	 done
	 gunzip *.fa.gz

	 cd ../mm10
	 printf "\nDownloading mm10 to data/genomes/.\n"
	 for c in {1..19} X Y M
	 do
             printf "chr${c}\n"
	     wget -q "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr${c}.fa.gz"
         done
	 gunzip *.fa.gz
	 cd ../../../../
	 break;;
	[Nn]* ) exit;;
	* ) printf "enter y or n; Download hg38 and mm10 genomes?\n";;
	esac
done

