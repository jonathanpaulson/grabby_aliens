import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import subprocess
import math

parser = argparse.ArgumentParser()
parser.add_argument('--D', type=int)
parser.add_argument('--n', type=float)
parser.add_argument('--N', type=str)
parser.add_argument('--c', type=float, default=1.0)
parser.add_argument('--L', type=float, default=1.0)
parser.add_argument('--s', type=float, default=1.0)
parser.add_argument('--seed', type=float, default=0)
parser.add_argument('--empty_samples', type=int, default=0)

args = parser.parse_args()
D,n,N,s,c,L,seed,empty_samples = args.D,args.n,int(float(args.N)),args.s,args.c,args.L,args.seed,args.empty_samples
if s != 1.0:
    print('WARNING: s!=1.0 so graph titles are misleading')

fname = f'D={D}_n={n}_N={N:.2e}_L={L}_c={c}_seed={seed}'
subprocess.check_output(f'g++ -std=c++17 -O3 -Wall -Werror -Wextra -Wshadow -Wno-sign-compare simulate.cc && ./a.out {D} {n} {N} {s} {c} {L} {fname} {seed} {empty_samples}', shell=True)

print(f'{fname}.csv')

XT = []
T = []
W = []
A = []
E = []
S = []
with open(f'{fname}.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        T.append(float(row['OriginTime']))
        W.append(float(row['MinArrival']))
        A.append(float(row['MaxAngle']))
        E.append(float(row['PctEmpty']))
        S.append(float(row['MinSee']))
        XT.append((float(row['X']), float(row['OriginTime'])))

YEARS = []
with open(f'{fname}_years.txt') as yearfile:
    reader = csv.DictReader(yearfile)
    for row in reader:
        YEARS.append(row)

T50 = np.median(T)
TS = [x/T50 for x in T]
WS = [x/T50 for x in W]
SS = [x/T50 for x in S]

print(np.mean(E))

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

civs_x = [float(i)/len(T) for i in range(len(T))]
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

plot(p[3,0], civs_x, sorted(A), 'MaxAngle', log=False)
if any([e for e in E]):
    plot(p[3,1], civs_x, list(reversed(sorted(E))), '% Empty', log=False)
else:
    fig.delaxes(p[3,1])

plt.savefig(f'{fname}.png')
subprocess.check_output(f'cmd.exe /C start {fname}.png', shell=True)
#subprocess.check_output(f'display {fname}.png', shell=True)
