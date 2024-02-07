import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data = pd.read_csv('TissueExpConcDifferential.csv', sep=',',header=0)

data['Kidney'].plot(kind='hist', bins = 400)

plt.show()
