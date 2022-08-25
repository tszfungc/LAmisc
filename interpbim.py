import argparse
import logging
import sys

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)


def parse_args(argv):
    parser = argparse.ArgumentParser(description="fit LAD decay curve")

    parser.add_argument("--cm-map", help="Genetic map")
    parser.add_argument("--bim", help="bim path")
    parser.add_argument("--out", help="output bim")

    return parser.parse_args(argv)


def interpolate_bim(bim: pd.DataFrame, cm_map: pd.DataFrame) -> pd.DataFrame:
    """Fill genetic distance in .bim given genetic map

    Parameters
    ----------
    bim:
        DataFrame read from plink .bim file
        https://www.cog-genomics.org/plink/1.9/formats#bim
    cm_map:
        DataFrame of genetic map
        chrom, position, cM.
        No header
        tab-delimited

    Returns
    -------
    DataFrame in plink .bim file. Genetic distance are interpolated

    """

    chrom = bim.iloc[0, 0]
    cm_map = cm_map[cm_map.iloc[0] == chrom]
    interp_cM = np.interp(bim[3].values, cm_map.iloc[:, 1], cm_map.iloc[:, 2])

    bim_out = bim.copy()
    bim_out.iloc[:, 2] = interp_cM

    return bim_out


def main(argv):
    args = parse_args(argv)

    bim = pd.read_csv(args.bim, sep="\t", header=None)
    cm_map = pd.read_csv(args.cm_map, sep="\t", header=None)

    bim_fill = interpolate_bim(bim, cm_map)
    bim_fill.to_csv(args.out, sep="\t", header=None, index=False)

    sys.exit(0)


if __name__ == "__main__":
    main(sys.argv[1:])
