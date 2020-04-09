## Examples

```
git clone https://github.com/vanheeringen-lab/ANANSE.git
cd ANANSE/test/
```
### Build TF binding network

```
$ ananse binding  -r data/krt_enhancer.bed \
                    -o results/binding.txt \
                    -a /data/hg38_genes.bed \
                    -g hg38 \
                    -p /data/gimme.vertebrate.v5.1.pfm
```

### Built gene regulatory network

```
$ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                    -r data/krt_enhancer.bed \
                    -b results/binding.txt \
                    -o results/full_features.txt \
                    -a /data/hg38_genes.bed \
                    -g hg38 \
                    -c expressioncorrelation.txt \
                    -p ../data/gimme.vertebrate.v5.1.pfm
```

### Infer TF influence score

```
$ ananse influence  -a results/full_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -p False
```
