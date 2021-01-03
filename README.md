time python3 analyze.py --D 3 --n 0 --N 1e9 --c 2 --L 1

To only run a sim (creating a csv file) and not create any figures, you can use `simulate.py`. This simply provides a friendlier interface for calling the C++ code:
```
$ python simulate.py --help
usage: simulate.py [-h] --D D --n N --N N [--s S] [--c C] [--L L] [--dir DIR]

Call the high performance c++ simulator code and save the output to a csv
file.

optional arguments:
  -h, --help  show this help message and exit
  --D D       the dimension of space
  --n N       the power of the power-law civ origin time distribution t^n
  --N N       the number of civs to sample
  --s S       the speed of alien expansion
  --c C       the speed of light
  --L L       the width of the universe
  --dir DIR   directory to save csv to
```
