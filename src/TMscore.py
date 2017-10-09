import re
import os.path
import subprocess


class TMscore():
    def __init__(self, path):
        if os.path.isfile(path):
            self.path = path
            pass
        else:
            raise Exception("%s was not found" % path)

    def __call__(self, prot_a, prot_b):
        if os.path.isfile(prot_a) and os.path.isfile(prot_b):
            out = subprocess.check_output([self.path, prot_a, prot_b])
            data = str(out).split("\\n")
            for d in data:
                x = re.sub("\s\s+", " ", d).split(' ')
                if x[0] == "TM-score" and x[1] == "=":
                    return float(x[2])
        else:
            raise Exception("Check that %s and %s exists" % (prot_a, prot_b))


if __name__ == "__main__":
    tmscore = TMscore("/home/h3nnn4n/TMscore.py/src/TMscore")
    print(tmscore("/home/h3nnn4n/1crn.pdb", "/home/h3nnn4n/pspp.jl/src/simulated_annealing.py/best.pdb"))
