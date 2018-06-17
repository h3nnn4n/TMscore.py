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

        self.rmsd = None
        self.gdt_ts = None
        self.gdt_ha = None
        self.tm_score = None

    def __call__(self, prot_a, prot_b):
        if os.path.isfile(prot_a) and os.path.isfile(prot_b):
            out = subprocess.check_output([self.path, prot_a, prot_b])
            data = str(out).split("\\n")
            for d in data:
                x = re.sub("\s\s+", " ", d).split(' ')
                if x[0] == "TM-score" and x[1] == "=":
                    self.tm_score = float(x[2])
                elif x[0] == "GDT-TS-score=":
                    self.gdt_ts = float(x[1])
                elif x[0] == "GDT-HA-score=":
                    self.gdt_ha = float(x[1])
                elif x[0] == "RMSD":
                    self.rmsd = float(x[5])
        else:
            raise Exception("Check that %s and %s exists" % (prot_a, prot_b))

    def info(self):
        print(self.rmsd)
        print(self.tm_score)
        print(self.gdt_ts)
        print(self.gdt_ha)


if __name__ == "__main__":
    tmscore = TMscore("./TMscore")
    tmscore("./1crn.pdb", "./best.pdb")
    tmscore.info()
