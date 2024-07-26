import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, t
import os

results = []
for filename in os.listdir("final_ts"):
    filename = "final_ts/" + filename
    df = pd.read_csv(filename, header=None)

    x = np.arange(df.size)
    y = df.values.ravel()

    for i in range(1, df.size):
        yi = y[::i]
        xi = x[::i]
        m = ((yi - yi.mean()) * (xi - xi.mean())).sum() / ((xi - xi.mean()) ** 2).sum()
        q = yi.mean() - m * xi.mean()
        res = yi - (q + m * xi)
        if pearsonr(res[1:], res[:-1])[1] > 0.01:
            sse = (res**2).sum()
            tst = m * np.sqrt((xi.size - 2) * ((xi - xi.mean()) ** 2).sum() / sse)
            p = (1 - t.cdf(abs(tst), (xi.size - 2))) * 2
            results.append([filename, m, p, i])
            break

results = pd.DataFrame(data=results, columns=["filename", "trend", "pval", "delta_t"])
PTHRESHOLD = .01
significant = results[results.pval < PTHRESHOLD / len(results)]
