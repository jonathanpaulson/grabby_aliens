There are two files: analyze.py and simulate.cc. You are intended to run analyze.py, which compiles and runs simulate.cc itself.

# Dependencies
You will need python3 and g++. I am using "Python 3.6.9" and "g++ (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0", although any recent version should work.

# Inputs
analyze.py takes the following arguments<br>
--n The power in the origin-time power-law. You may pass a comma-separated list of values<br>
--N The number of potential civilizations.<br>
--sc The ratio s/c - how fast civs expand relative to the speed of light. You may pass a comma-separated list of values.<br>
--s The speed of civilization expansion (default 1.0)<br>
--m The power in the universe expansion scale factor (default 2/3; see section "8 Cosmology")<br>
--D The dimensions of space. Should be 1,2 or 3 (default 3)<br>
--L The size of the universe (default 1.0)<br>
--seed A random seed (default 0)<br>
--empty_samples How precisely to estimate how full the universe is (default 0, meaning no estimate at all)<br>

# Outputs
Each run of simulate.cc will produce two output files, one ending in "civs.csv" and one ending in "years.csv"<br>
analyze.py will run simulate.cc once for each combination of parameters you give it. For instance, if you give "--n 10,12" and "--c 1.0,2.0",
then simulate.cc will run four times: once with n=10,c=1.0, once with n=10,c=2.0, once with n=12,c=1.0 and once with n=12,c=2.0

The "civs" file contains the following columns, with one row per "surviving" civilization:<br>
- D columns (X,Y,Z) for the position in space<br>
- OriginTime: The model/conformal OriginTime<br>
- MinArrival: The model time when another civilization first arrives at our origin point<br>
- MinSee: The model time when we first see signals from another civilization<br>
- NumberSeen: The number of alien civilizations we see signals from at our origin time<br>
- MaxAngle: How much of the sky another civilization takes up in our sky at our origin time (section 10 H; between 0 and pi)<br>
- PctEmpty: An estimate of how much of the universe is empty at our origin time<br>
- Volume: An estimate of what fraction of the universe this civ eventually controls<br>

The "years" file contains the following columns, all in "clock" years:<br>
- OriginTime: Sampled distribution of civilization origin times<br>
- MinWait: Sampled distribution of how long civilizations will have to wait to meet another civilization<br>
- MinSETI: Sampled distribution of how long civilizations will have to wait to see signals from another civilization<br>

analyze.py will combine statistics from these simulations into a graph, with different values of "n" as different datasets and different values of "c" as columns of graphs.

# Examples
To generate figure 13:
`time python3 analyze.py --D 3 --n 1.5,3,6,12 --N 1e8 --sc 1.0,0.75,0.5,0.25 --L 1 --s 1 --seed 0 --empty 100`

To generate table 1:
`time python3 analyze.py --D 3 --n 6,12 --N 1e8 --sc 0.5,0.75 --L 1 --s 1 --seed 0 --empty 100 --table_1`

To generate data for figure 12:
`time python3 analyze.py --D 3 --n 1.5,3,6,12 --N 1e8 --sc 1.0 --L 1 --s 1 --seed 0 --figure_12`

Help:
`python3 analyze.py --help`

You can reach me at jonathanpaulson@gmail.com
