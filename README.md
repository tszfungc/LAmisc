# LAmisc
Local ancestry miscellaneous

`pgenlib` is needed. Download [here](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)


## To Do

- Global ancestries
- make LAD matrix

## Converting .fb.tsv to plink format

```bash
# Output test.pgen test.psam test.pvar
python rfmix2pgen.py --file test.fb.tsv --out test --pop AFR
```
