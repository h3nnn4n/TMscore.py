import re
import os
import os.path
import errno
import subprocess


class TMscore():
    def __init__(self, path=None):
        self.setup_path(path)

        self.rmsd = None
        self.gdt_ts = None
        self.gdt_ts_info = None
        self.gdt_ha = None
        self.gdt_ha_info = None
        self.tm_score = None
        self.maxsub = None

    def setup_path(self, path):
        if path is None:
            path = './TMscore'

        if os.path.isfile(path):
            self.path = os.path.join(os.getcwd(), path)
        else:
            raise FileNotFoundError(
                errno.ENOENT,
                os.strerror(errno.ENOENT),
                path
            )

    def __call__(self, prot_a, prot_b):
        if os.path.isfile(prot_a) and os.path.isfile(prot_b):
            out = subprocess.check_output([self.path, prot_a, prot_b])
            data = str(out).split("\\n")
            for d in data:
                x = re.sub(r"\s\s+", " ", d).split(' ')
                if x[0] == "TM-score" and x[1] == "=":
                    self.tm_score = float(x[2])
                elif x[0] == "GDT-TS-score=":
                    self.gdt_ts = float(x[1])
                    a = float(x[2].split('=')[1])
                    b = float(x[3].split('=')[1])
                    c = float(x[4].split('=')[1])
                    d = float(x[5].split('=')[1])
                    self.gdt_ts_info = (a, b, c, d)
                elif x[0] == "GDT-HA-score=":
                    self.gdt_ha = float(x[1])
                    a = float(x[2].split('=')[1])
                    b = float(x[3].split('=')[1])
                    c = float(x[4].split('=')[1])
                    d = float(x[5].split('=')[1])
                    self.gdt_ha_info = (a, b, c, d)
                elif x[0] == "RMSD":
                    self.rmsd = float(x[5])
                elif x[0] == "MaxSub-score=":
                    self.maxsub = float(x[1])
        else:
            raise Exception("Check that %s and %s exists" % (prot_a, prot_b))

    def print_info(self):
        print(self.rmsd)
        print(self.tm_score)
        print(self.maxsub)
        print(self.get_gdt_ts_info())
        print(self.get_gdt_ha_info())

    def get_rmsd(self):
        return self.rmsd

    def get_gdt_ts(self):
        return self.gdt_ts

    def get_gdt_ts_info(self):
        return (self.gdt_ts, self.gdt_ts_info)

    def get_gdt_ha(self):
        return self.gdt_ha

    def get_gdt_ha_info(self):
        return (self.gdt_ha, self.gdt_ha_info)

    def get_tm_score(self):
        return self.tm_score

    def get_maxsub(self):
        return self.maxsub

    def get_all(self):
        a = {}
        a['rmsd'] = self.rmsd
        a['gdt_ts'] = self.get_gdt_ts_info()
        a['gdt_ha'] = self.get_gdt_ha_info()
        a['maxsub'] = self.maxsub
        a['tm_score'] = self.tm_score

        return a

if __name__ == "__main__":
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")
    tmscore.print_info()
