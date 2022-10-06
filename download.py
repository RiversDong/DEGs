import pandas as pd
from multiprocessing import Pool
import os 
def get_data(metafile, download_path):
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
