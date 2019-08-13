python ../grns/binding.py  -r data/krt_enhancer.bed \
                            -o results \
                            -a /home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed \
                            -g hg38 \
                            -p ../data/gimme.vertebrate.v5.1.pfm


python ../grns/interaction.py  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                            -o results \
                            -a /home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed \
                            -g hg38 \
                            -b results/binding.predicted.h5 \
                            -c /home/qxu/projects/regulatoryNetwork/history/cell_trans/human_gene_correlation/expressioncorrelation.txt \
                            -p ../data/gimme.vertebrate.v5.1.pfm

python ../grns/builtnetwork.py -f results/full_features.h5 -o results



