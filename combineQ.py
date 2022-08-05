import numpy as np
import sys
import pandas as pd
import argparse


def parse_args(argv):

    parser = argparse.ArgumentParser(
        description='Combine .rfmix.Q global ancestry from all chromosomes'
    )

    parser.add_argument('--file', help='Text file, list of global ancestry output .rfmix.Q')
    parser.add_argument('--out', help='Prefix of the output files')
    parser.add_argument('--pop', help='the ancestry to output, match to the header of rfmix output', type=str)

    return parser.parse_args(argv)



def main(argv):

    args = parse_args(argv)

    rfmixQ = args.file
    outprefix = args.out

    with open(rfmixQ, 'r') as f:
        lines = [l.strip() for l in f.readlines()]

    print(lines)

    dfs = [pd.read_csv(l, sep='\t', skiprows=1) for l in lines]

    Q = np.mean([df[args.pop] for df in dfs], axis=0)

    df_ = dfs[0][['#sample']].copy()
    df_[args.pop] = Q

    df_.to_csv(outprefix+".rfmix.Q", index=None, sep='\t')


if __name__ == '__main__':
    main(sys.argv[1:])
