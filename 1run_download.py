import os 
import pandas as pd
import sys
from multiprocessing import Pool

metafile = sys.argv[1]
download_path = sys.argv[2]

data = pd.read_csv(metafile, sep = ",")
runs = data["Run"]

cmds = []
for i in runs:
    cmd = "prefetch {0} --output-directory {1}".format(i, download_path)
    cmds.append(cmd)

pool = Pool(6)
for i in cmds:
    pool.apply_async(os.system, args=(i,))
pool.close()
pool.join()
