# TMscore.py
------------
A simple class to wrap the original TMscore Fortran77 code.
The source code for TMscore (required) can be found here: https://zhanglab.ccmb.med.umich.edu/TM-score/
It is also available in this repository.

## Usage

``` Python
from TMscore import TMscore


tmscore = TMscore(path_to_the_binary)  # Point where the binary is
print(tmscore(pdb1, pdb2))  # paths for 2 pdb files for comparison
```

## License
See LICENSE. For more info on TMscore see https://zhanglab.ccmb.med.umich.edu/TM-score/
