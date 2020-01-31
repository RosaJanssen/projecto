import numpy as np 
import matplotlib as mpl 
import pandas as pd
import statistics
from scipy.ndimage.filters import gaussian_filter1d


## agg backend is used to create plot as a .png file
#mpl.use('agg')

import matplotlib.pyplot as plt

data_to_plot_la2 = pd.read_csv('resultfile_lookahead_2', sep=',', index_col =0)
data_to_plot_la3 = pd.read_csv('resultfile_lookahead_drie', sep=',', index_col =0)
data_to_plot_la4 = pd.read_csv('resultfile_lookahead_vier', sep=',', index_col =0)
data_to_plot_la5 = pd.read_csv('resultfile_lookahead_vijf', sep=',', index_col =0)
data_to_plot_la6 = pd.read_csv('resultfile_lookahead_zes', sep=',', index_col =0)

averages = [data_to_plot_la2.Stability.mean(), data_to_plot_la3.Stability.mean(), data_to_plot_la4.Stability.mean(), data_to_plot_la5.Stability.mean(), data_to_plot_la6.Stability.mean()]
ysmoothed = gaussian_filter1d(averages, sigma=1)

lookaheads = [2,3,4,5,6]





# Set interval on x-axis 
#plt.xticks(np.arange(min(data["Stability"]), max(data["Stability"])+1, 1.0))

# Name labels
plt.xlabel('Number of steps')
plt.ylabel('Stability')
#plt.xticks(np.arange(0,6, 1.0))
plt.xticks([2, 3, 4, 5, 6], ['2', '3', '4', '5', '6'])

# Name title
plt.title('Lookahead algorithm')
plt.plot([2, 3, 4, 5, 6], ysmoothed)

plt.show()