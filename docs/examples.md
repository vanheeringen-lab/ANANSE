## Examples

```
git clone https://github.com/vanheeringen-lab/ANANSE.git
cd ANANSE/test/
ls -lh data     

# total 3.2M
# -rw-rw-r-- 1 qxu qxu 730K Apr  2 14:51 FB2KRT_degenes.csv
# -rw-rw-r-- 1 qxu qxu 1.2M Apr 11 21:22 FB_enhancer.bed
# -rw-rw-r-- 1 qxu qxu 260K Apr  2 14:51 FB_rep1_TPM.txt
# -rw-rw-r-- 1 qxu qxu 261K Apr 11 21:22 FB_rep2_TPM.txt
# -rw-rw-r-- 1 qxu qxu 253K Apr 11 21:22 KRT_enhancer.bed
# -rw-rw-r-- 1 qxu qxu 267K Apr  2 14:51 KRT_rep1_TPM.txt
# -rw-rw-r-- 1 qxu qxu 268K Apr  2 14:51 KRT_rep2_TPM.txt
```
### Build TF binding network

``` bash
$ ananse binding  -r data/krt_enhancer.bed \
                  -o results/binding.txt \
                  -a ../data/hg38_genes.bed \
                  -g hg38 \
                  -p ../data/gimme.vertebrate.v5.1.pfm
```

### Built gene regulatory network

``` bash
$ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                  -b results/binding.txt \
                  -o results/KRT_features.txt \
                  -g hg38 \
                  -a ../data/hg38_genes.bed 
```

### Infer TF influence score

``` bash
$ ananse influence  -b results/FB_network.txt \
                    -a results/KRT_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -p 
```
