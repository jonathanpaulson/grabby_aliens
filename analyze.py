import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import subprocess
import math

parser = argparse.ArgumentParser(description='Run Grabby Aliens model')
parser.add_argument('--D', type=int, help='Dimensions of space')
parser.add_argument('--n', type=str, help='Origin time power-law; can pass multiple comma-separated values')
parser.add_argument('--N', type=str, help='Number of potential civilizations')
parser.add_argument('--c', type=float, default=1.0, help='The speed of light')
parser.add_argument('--L', type=float, default=1.0, help='The size of the universe')
parser.add_argument('--s', type=float, default=1.0, help='The speed of civ expansion')
parser.add_argument('--seed', type=float, default=0, help='A random seed')
parser.add_argument('--empty_samples', type=int, default=0, help='How many points to sample when estimating how full the universe is')

args = parser.parse_args()
D,n,N,s,c,L,seed,empty_samples = args.D,args.n,int(float(args.N)),args.s,args.c,args.L,args.seed,args.empty_samples
ns = [float(n) for n in n.split(',')]

DATA = {}
for n in ns:
    fname = f'D={D}_n={n}_N={N:.2e}_L={L}_c={c}_seed={seed}'
    subprocess.check_output(f'g++ -std=c++17 -O3 -Wall -Werror -Wextra -Wshadow -Wno-sign-compare simulate.cc && ./a.out {D} {n} {N} {s} {c} {L} {fname} {seed} {empty_samples}', shell=True)
    print(f'Generated data in {fname}.csv and {fname}_years.csv')

    # Read CIV data
    with open(f'{fname}.csv') as csvfile:
        CIVS = list(csv.DictReader(csvfile))

    # Read years data
    with open(f'{fname}_years.csv') as yearfile:
        YEARS = list(csv.DictReader(yearfile))

    DATA[n] = (CIVS, YEARS)

C = ','.join([str(len(CIVS)) for (CIVS,YEARS) in DATA.values()])
# Make graphs
fig, p = plt.subplots(4,2,constrained_layout=True,figsize=(12,12))
fig.suptitle(f'D={D} n={n} N={N:.2e} L={L} c={c} s={s} |C|={C}')

def plot(ax, label, log):
    ax.minorticks_on()
    ax.grid(b=True, which='major', axis='both')
    ax.tick_params(axis='both', which='both', bottom=True, left=True)
    ax.set_xlabel('Percentile')
    for n,(CIVS,YEARS) in DATA.items():
        civs_x = [float(i)/len(CIVS) for i in range(len(CIVS))]
        years_x = [float(i)/len(YEARS) for i in range(len(YEARS))]
        # Rescale model times so median(Origin)=1.0
        T50 = np.median([float(row['OriginTime']) for row in CIVS])
        if label == 'Origin':
            x = civs_x
            y = [float(row['OriginTime'])/T50 for row in CIVS]
        elif label == 'Origin (Gyr)':
            x = years_x
            y = sorted([float(row['OriginTime']) for row in YEARS])
        elif label == 'MinArrival':
            x = civs_x
            y = sorted([float(row['MinArrival'])/T50 for row in CIVS])
        elif label == 'MinWait (Gyr)':
            x = years_x
            y = sorted([float(row['MinWait']) for row in YEARS])
        elif label == 'MinSee':
            x = civs_x
            y = sorted([float(row['MinSee'])/T50 for row in CIVS])
        elif label == 'MinSETI (Gyr)':
            x = years_x
            y = sorted([float(row['MinSETI']) for row in YEARS])
        elif label == 'MaxAngle':
            x = civs_x
            y = sorted([float(row['MaxAngle']) for row in CIVS])
        elif label == '% Empty':
            x = civs_x
            y = list(reversed(sorted([float(row['PctEmpty']) for row in CIVS])))
        else:
            assert False, f'Unknown label={label}'
        ax.plot(x, y, label=f'n={n}')
    ax.legend(loc='lower right');
    ax.set_ylabel(label)
    if log:
        ax.set_yscale('log')

plot(p[0,0], 'Origin', log=False)
plot(p[0,1], 'Origin (Gyr)', log=True)
plot(p[1,0], 'MinArrival', log=False)
plot(p[1,1], 'MinWait (Gyr)', log=True)
plot(p[2,0], 'MinSee', log=False)
plot(p[2,1], 'MinSETI (Gyr)', log=True)
plot(p[3,0], 'MaxAngle', log=False)
if empty_samples:
    plot(p[3,1], '% Empty', log=False)
else:
    fig.delaxes(p[3,1])

plt.savefig(f'{fname}.png')
# Open PNG in windows
# Switch which line is commented for Linux
subprocess.check_output(f'cmd.exe /C start {fname}.png', shell=True)
#subprocess.check_output(f'display {fname}.png', shell=True)
