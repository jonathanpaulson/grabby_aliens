import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import argparse
import subprocess
import math

parser = argparse.ArgumentParser()
parser.add_argument('--D', type=int)
parser.add_argument('--n', type=float)
parser.add_argument('--N', type=int)
parser.add_argument('--s', type=float)

args = parser.parse_args()
D,n,N,s = args.D,args.n,args.N,1.0/args.s

fname = f'D={D}_n={n}_N={N:.2e}_s={s}'
subprocess.check_output(f'./simulate {D} {n} {N} {1.0/s} > {fname}.csv', shell=True)

T = []
W = []
A = []
with open(f'{fname}.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        T.append(float(row['OriginTime']))
        W.append(float(row['MinWait']))
        A.append(float(row['MaxAngle']))

T50 = np.median(T)
TS = [x/T50 for x in T]
WS = [x/T50 for x in W]
AS = [x/T50 for x in A]

fig, p = plt.subplots(3)
fig.suptitle(f'D={D} n={n} N={N:.2e} s={s} |C|={len(TS)} |C|/s^D={len(TS)/s**D}')

p[0].plot(TS)
p[0].set_ylabel('OriginTime')
p[0].set_xlabel('Index')

p[1].plot(sorted(WS))
p[1].set_ylabel('MinWait')
p[1].set_xlabel('Index')

p[2].plot(sorted(AS))
p[2].set_ylabel('MaxAngle')
p[2].set_xlabel('Index')

plt.savefig(f'{fname}.png')
subprocess.check_output(f'cmd.exe /C start {fname}.png', shell=True)
