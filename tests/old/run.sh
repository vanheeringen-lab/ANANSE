ananse binding  -n 30 \
                -r data/FB_enhancer.bed \
                -o results/FB_binding.txt \
                -g hg38 \
                --unremove-curated

ananse binding  -n 30 \
                -r data/KRT_enhancer.bed \
                -o results/KRT_binding.txt \
                -g hg38 \
                --unremove-curated

ananse network  -n 30 \
                -e data/FB_rep1_TPM.txt data/FB_rep2_TPM.txt \
                -b results/FB_binding.txt \
                -o results/FB_network.txt \
                -g hg38 \
                --exclude-promoter --include-enhancer

ananse network  -n 30 \
                -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                -b results/KRT_binding.txt \
                -o results/KRT_network.txt \
                -g hg38 \
                --exclude-promoter --include-enhancer

ananse influence    -n 30 \
                    -s results/FB_network.txt \
                    -t results/KRT_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -i 100000
