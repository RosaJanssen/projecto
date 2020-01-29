import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Walk through data created in the algorithm
data = pd.read_csv('resultfile', sep=',', index_col =0)

# Make histogram
plt.hist(data["Stability"], bins=range(int(max(data["Stability"])+ 1)), color='#0504aa', alpha=0.7, rwidth=0.85)

# Set interval on x-axis 
plt.xticks(np.arange(min(data["Stability"]), max(data["Stability"])+1, 1.0))

# Name labels
plt.ylabel('Frequency')
plt.xlabel('Stability')

# Name title
plt.title('Lookahead algorithm')

plt.show()