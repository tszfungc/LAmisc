# LAmisc
Local ancestry miscellaneous

`pgenlib` is needed. Download [here](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)

## Converting .fb.tsv to plink format

```bash
# Output test.pgen test.psam test.pvar
## Convert RFMIX2 .fb.tsv to plink format
python rfmix2pgen.py --file test.fb.tsv --out test --pop AFR

# Run genome-wide scan on example phenos
plink2 --pfile example --pheno example.pheno --glm allow-no-covars --out example
```
