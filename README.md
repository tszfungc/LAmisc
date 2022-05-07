# LAmisc
Local ancestry miscellaneous

## Convert RFMIX2 .fb.tsv to plink format

```bash
# Convert to plink2
python rfmix2pgen.py --file example.rfmix.fb.tsv --out example --popidx 1

# Run genome-wide scan on example phenos
plink2 --pfile example --pheno example.pheno --glm allow-no-covars --out example
```
