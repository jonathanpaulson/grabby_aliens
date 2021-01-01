import argparse
import pathlib

import subprocess


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Call the high performance c++ simulator code and save the output to a csv file.')
    parser.add_argument('--D', type=int, required=True,
        help="the dimension of space")
    parser.add_argument('--n', type=float, required=True,
        help="the power of the power-law civ origin time distribution t^n")
    parser.add_argument('--N', type=str, required=True,
        help="the number of civs to sample")
    parser.add_argument('--s', type=float, default=1.0,
        help="the speed of alien expansion")
    parser.add_argument('--c', type=float, default=1.0,
        help="the speed of light")
    parser.add_argument('--L', type=float, default=1.0,
        help="the width of the universe")
    parser.add_argument('--dir', type=str, default=".",
        help="directory to save csv to")

    args = parser.parse_args()
    D,n,N,s,c,L = args.D,args.n,int(float(args.N)),args.s,args.c,args.L
    # if s != 1.0:
    #     print('WARNING: s!=1.0 so graph titles are misleading')
    save_dir =  pathlib.Path(args.dir)
    if not save_dir.exists():
        save_dir.mkdir()
    fname = pathlib.Path(save_dir, f"D={D}_n={float(n)}_N={N:.2e}_s={s:.5e}_L={float(L)}_c={float(c)}.csv")

    subprocess.check_output(f'g++ -std=c++17 -O3 -Wall -Werror -Wextra -Wshadow -Wno-sign-compare simulate.cc && ./a.out {D} {n} {N} {s} {c} {L} > {fname}', shell=True)

