import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_csv('data.csv',sep='\t').sort_values('complete')
df['wait'] = df['launched'] - df['added']
df['delay'] = df['starting'] - df['launched']
df['runtime'] = df['complete'] - df['starting']



df['runtime'].plot(kind='hist',bins=100)



areafig, ax1 = plt.subplots()
x= np.arange(len(df['added'].values))
y =df.loc[:,['added','wait','delay','runtime']].T.values
ax1.stackplot(x,y,edgecolor='none',labels=['added','wait','loading','runtime'])



cumfig, ax2 = plt.subplots()
df2 = pd.DataFrame({'const': 1,
                   'time': df.starting.values,
                    'task': df.task.values})
df3 = pd.DataFrame({'const': -1,
                   'time': df.complete.values,
                   'task': df.task.values})
df2=df2.append(df3)
df2.sort_values('time',inplace=True)
df2['cumsum']= df2['const'].cumsum()
ax2.plot(df2['time'],df2['cumsum'].values)

plt.legend()
plt.tight_layout()
plt.show()

