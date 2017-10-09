# TMscore.py
------------
A simple class to wrap the original TMscore Fortran77 code.
The source code for TMscore (required) can be found here: https://zhanglab.ccmb.med.umich.edu/TM-score/

## Usage

``` Python
import TMscore


tmscore = TMscore("/home/h3nnn4n/TMscore.py/src/TMscore")
print(tmscore("/home/h3nnn4n/1crn.pdb", "/home/h3nnn4n/pspp.jl/src/simulated_annealing.py/best.pdb"))
```

## License
See LICENSE. For more info on TMscore see https://zhanglab.ccmb.med.umich.edu/TM-score/
