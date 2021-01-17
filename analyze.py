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
parser.add_argument('--n', type=float, help='Origin time power-law')
parser.add_argument('--N', type=str, help='Number of potential civilizations')
parser.add_argument('--c', type=float, default=1.0, help='The speed of light')
parser.add_argument('--L', type=float, default=1.0, help='The size of the universe')
parser.add_argument('--s', type=float, default=1.0, help='The speed of civ expansion')
parser.add_argument('--seed', type=float, default=0, help='A random seed')
parser.add_argument('--empty_samples', type=int, default=0, help='How many points to sample when estimating how full the universe is')

args = parser.parse_args()
D,n,N,s,c,L,seed,empty_samples = args.D,args.n,int(float(args.N)),args.s,args.c,args.L,args.seed,args.empty_samples
if s != 1.0:
    print('WARNING: s!=1.0 so graph titles are misleading')

fname = f'D={D}_n={n}_N={N:.2e}_L={L}_c={c}_seed={seed}'
subprocess.check_output(f'g++ -std=c++17 -O3 -Wall -Werror -Wextra -Wshadow -Wno-sign-compare simulate.cc && ./a.out {D} {n} {N} {s} {c} {L} {fname} {seed} {empty_samples}', shell=True)

print(f'Generated data in {fname}.csv and {fname}_years.csv')

# Read CIV data
CIVS = []
with open(f'{fname}.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        CIVS.append(row)

# Read years data
YEARS = []
with open(f'{fname}_years.csv') as yearfile:
    reader = csv.DictReader(yearfile)
    for row in reader:
        YEARS.append(row)

# Rescale model times so median(Origin)=1.0
T50 = np.median([float(row['OriginTime']) for row in CIVS])
TS = [float(row['OriginTime'])/T50 for row in CIVS]
WS = [float(row['MinArrival'])/T50 for row in CIVS]
SS = [float(row['MinSee'])/T50 for row in CIVS]

# Make graphs
fig, p = plt.subplots(4,2,constrained_layout=True,figsize=(12,12))
fig.suptitle(f'D={D} n={n} N={N:.2e} L={L} |C|={len(TS)} |C|/L^D={len(TS)/L**D}')

def plot(ax, x, y, ylabel, log):
    ax.minorticks_on()
    ax.grid(b=True, which='major', axis='both')
    ax.tick_params(axis='both', which='both', bottom=True, left=True)
    ax.set_xlabel('Percentile')
    ax.plot(x, y)
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale('log')

civs_x = [float(i)/len(CIVS) for i in range(len(CIVS))]
years_x = [float(i)/len(YEARS) for i in range(len(YEARS))]

plot(p[0,0], civs_x, TS, 'Origin', log=False)
origin_years = sorted([row['OriginTime'] for row in YEARS])
plot(p[0,1], years_x, origin_years, 'Origin (BYr)', log=True)

plot(p[1,0], civs_x, sorted(WS), 'MinArrival', log=False)
wait_years = sorted([row['MinWait'] for row in YEARS])
plot(p[1,1], years_x, wait_years, 'MinWait (BYr)', log=True)

plot(p[2,0], civs_x, sorted(SS), 'MinSee', log=False)

seti_years = sorted([row['MinSETI'] for row in YEARS])
plot(p[2,1], years_x, seti_years, 'MinSETI (BYr)', log=True)

angles = sorted([float(row['MaxAngle']) for row in CIVS])
plot(p[3,0], civs_x, angles, 'MaxAngle', log=False)
E = [float(row['PctEmpty']) for row in CIVS]
if any([e for e in E]): # if we generated PctEmpty data
    plot(p[3,1], civs_x, list(reversed(sorted(E))), '% Empty', log=False)
else:
    fig.delaxes(p[3,1])

plt.savefig(f'{fname}.png')
# Open PNG in windows
# Switch which line is commented for Linux
subprocess.check_output(f'cmd.exe /C start {fname}.png', shell=True)
#subprocess.check_output(f'display {fname}.png', shell=True)
