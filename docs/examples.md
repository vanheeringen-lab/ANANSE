## Examples

### Prepare code and dataset
* Activate the environment before you use ANANSE  
```
conda activate ananse
```

* Install the `ANANSE` version used in the paper  
```
conda uninstall ananse
pip install git+https://github.com/vanheeringen-lab/ANANSE.git@d8dbdd2e405558334566386b747867c401f45870
```

* Download the origin data used in the paper  
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

mkdir results
```

### Build TF binding network

``` bash
$ ananse binding  -n 20 \
                  -r data/krt_enhancer.bed \
                  -o results/binding.txt \
                  -g hg38 \
                  --unremove-curated
```

### Built gene regulatory network

``` bash
$ ananse network  -n 20 \
                  -e data/FB_rep1_TPM.txt data/FB_rep2_TPM.txt \
                  -b results/binding.txt \
                  -o results/FB_network.txt \
                  -g hg38 \
                  --exclude-promoter --include-enhancer

$ ananse network  -n 20 \
                  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                  -b results/binding.txt \
                  -o results/KRT_network.txt \
                  -g hg38 \
                  --exclude-promoter --include-enhancer
```

### Infer TF influence score

``` bash
$ ananse influence  -n 20 \
                    -s results/FB_network.txt \
                    -t results/KRT_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -i 100000 
```
