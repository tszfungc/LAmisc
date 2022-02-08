import sys
import numpy as np
import pandas as pd
import logging
import argparse

logging.basicConfig(level=logging.INFO)

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Compute local ancestry linkage disequilibrium matrix")

    parser.add_argument("--cm-map", help='recombination map')
    parser.add_argument("--extract", help="list of coordinates to include")
    parser.add_argument("--generations", help="number of generations since admixture",
                        type=int)
    parser.add_argument("--out", help="output prefix")

    return parser.parse_args(argv)



def make_LADmat(args):
    cm_map = pd.read_csv(args.cm_map, sep='\t', header=None)


    # generate interpolated map
    extract_coord_sorted = np.sort(np.loadtxt(args.extract))
    extract_cm = np.interp(extract_coord_sorted, cm_map.iloc[:, 1], cm_map.iloc[:, 2])

    interp_map = pd.DataFrame({0:np.repeat(cm_map.iloc[0, 0], extract_coord_sorted.shape[0]),
                       1:extract_coord_sorted,
                       2:extract_cm})




    # compute LADmat
    lambda_mat = np.abs(extract_cm[:, None] - extract_cm[None, :])
    LADmat = np.exp(-0.01*args.generations*lambda_mat)

    return interp_map,LADmat

def main(argv):
    args = parse_args(argv)

    logging.info(f'''
    ==============Start==============

    Read genetic map: {args.cm_map}
    Read coordinates of target markers: {args.extract}
    Assumed number of generations since admixture: {args.generations}
    Interpolated genetic map will be output to : {args.out}.map
    NumPy matrix of LAD matrix be output to : {args.out}.npy

    ''')

    logging.info(f'Interp')
    interp_map, LADmat = make_LADmat(args)

    logging.info(f'Writing Output files')
    np.savetxt(f'{args.out}.map', interp_map, fmt=('%d', '%d', '%f'), delimiter='\t')
    np.save(f'{args.out}.npy', LADmat)

    logging.info(f'''
    ===============END================
    ''')


if __name__ == '__main__':
    main(sys.argv[1:])
