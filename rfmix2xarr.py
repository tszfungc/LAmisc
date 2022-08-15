"""
Python script to convert rfmix output to plink2 pgen format

The .fb.tsv file is read and dosages are written in plink2 pgen format.

Three files will be output

.pgen .pvar and .psam

In .pgen, the ancestry probability of the two haplotypes are summed. The resulted dosages are written to .pgen.

In .pvar, variant IDs are default to be chromosome:coordinate

.psam default sex is missing.

"""


import numpy as np, pandas as pd
import xarray as xr
import logging
import sys
import argparse

logging.basicConfig(level=logging.INFO)

def parse_args(argv):

    parser = argparse.ArgumentParser(
        description='Convert rfmix .fb.tsv output to plink2 pgen format'
    )

    parser.add_argument('--file', help='Path to rfmix .fb.tsv file')
    parser.add_argument('--out', help='Prefix of the output files')
    parser.add_argument('--pop', help='the ancestry to output, match to the header of rfmix output', type=str)


    return parser.parse_args(argv)


def main(argv):

    args = parse_args(argv)

    fbtsv = args.file
    outprefix = args.out


    # Count M
    f = open(fbtsv, 'r')
    M = sum(1 for line in f) - 2
    f.close()

    # Count n_pops
    f = open(fbtsv, 'r')
    comment = f.readline()
    pops = comment.strip().split('\t')[1:]
    pops_idx_lookup = {p:(i+1) for i,p in enumerate(pops) }
    idx = pops_idx_lookup[args.pop]

    n_pops = len(pops)

    # Count N
    header = f.readline()
    indv = np.array(
        list(map(lambda x: x.split(':::')[0], header.strip().split('\t')[4:]))[::(2*n_pops)]
    )
    N = indv.shape[0]

    logging.info(f'''
    ==============Start==============

    Read file: {fbtsv}
    Found {N} individuals and {M} variants
    in {n_pops} populations: {','.join(pops)}

    Local ancestry for \"{pops[idx-1]}\" will be output

    ''')

    # pgen
    chrom = None
    pos = []

    #M = 100
    genotype_matrix = np.nan * np.empty((M, N))
    for i, line in enumerate(f):
        if i % 100 == 0:
           logging.info(f'processing {i}-th marker')
    #    if i == 100:
    #        break
        line_split = line.strip().split('\t')
        genotype_matrix[i] = \
            np.float_(line_split[(3+idx)::(2*n_pops)]) + np.float_(line_split[(3+n_pops+idx)::(2*n_pops)])
        if chrom is None:
            chrom = int(line_split[0])
        pos.append(int(line_split[1]))

    vp = np.array(pos)
    vi = [ f'{chrom}:{p}' for p in pos ]
    si = indv

    xarr = xr.DataArray(
        genotype_matrix,
        dims= ['variants', 'samples'],
        coords={
            "pos": ('variants', pos),
            "variant_id": ('variants', vi),
            "sample_id": ('samples', si),
        }
    )


    print(xarr)

    xarr.to_netcdf(f'{args.out}.nc')

    logging.info(f'''
    ===============END================
    ''')

if __name__ == '__main__':
    main(sys.argv[1:])
