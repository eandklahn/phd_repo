import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from process_ac import *

def get_single_line(filename, line_number):

    f = open(filename, 'r')
    for n in range(line_number):
        line = f.readline()
    f.close()
    
    return line


filename1 = '20180209DyII_ac0.dat'
filename2 = '20171111BF4_AC1000.dat'
start_of_data = 20

df = pd.read_csv(filepath_or_buffer=filename1,
                 header=start_of_data
                 )

for h in df.columns:
    if np.all(np.isnan(df[h])):
        df.drop(h, axis=1, inplace=True)

temp_subsets = []
fq_range = set(df['Frequency (Hz)'])
num_meas_freqs = len(fq_range)
num_meas_temps = int(len(df['Temperature (K)'])/num_meas_freqs)

frequencies = df['Frequency (Hz)'][:num_meas_freqs]
temperature_subsets = []
for n in range(num_meas_temps):
    temp_subsets.append(df.iloc[n*num_meas_freqs:n*num_meas_freqs+num_meas_freqs])

temps = [subset['Temperature (K)'].mean() for subset in temp_subsets]

print(temp_subsets[0])
temp_subsets = []
new_column = np.ones(len(df[df.columns[0]]))
df['ones'] = new_column
for n in range(num_meas_temps):
    temp_subsets.append(df.iloc[n*num_meas_freqs:n*num_meas_freqs+num_meas_freqs])
print(temp_subsets[0])

#processes = 1
#
#ccinput = open('ccinput.dat', 'w')
#
#ccinput.write('{} {} {}\n'.format(processes, num_meas_temps, num_meas_freqs))
#
#ccinput.write('parameters guesses\n')
#
#for n in range(num_meas_freqs):
#
#    what_to_write = [frequencies[n]]
#
#    for i, sub in enumerate(temp_subsets):
#        what_to_write += [sub["M' (emu)"][i*num_meas_freqs+n], sub["M'' (emu)"][i*num_meas_freqs+n]]
#    what_to_write.append('\n')
#
#    line = ' '.join(str(i) for i in what_to_write)
#    ccinput.write(line)
#
#ccinput.close()
    
#set = 0
#sub = temp_subsets[set]
#freq = sub['Frequency (Hz)']
#mpp = sub["M'' (emu)"]
#fig, ax = plt.subplots()
#ax.plot(freq, mpp, label="M''")
#ax.plot(freq, Xpp_(freq, 100, 1, 0.01, 3), label=r"$\chi$''")
#ax.legend()
#ax.set_xscale('log')
#plt.show()

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#V, T = np.meshgrid(np.log10(frequencies), temps)
#M = np.array(df["M'' (emu)"]).reshape((num_meas_temps, num_meas_freqs))
#surf = ax.plot_surface(T, V, M, cmap=cm.coolwarm)
#
#
#ax.set_xlabel('Temp')
#ax.set_ylabel('Freq')
#ax.set_zlabel("M''")
#
#old_ticks = ax.get_yticks()
#new_ticks = 10**old_ticks
#ax.set_yticklabels('{:3.1e}'.format(i) for i in new_ticks)
#
#fig.colorbar(surf)
#plt.show()

#fig, ax = plt.subplots()
#for set in range(num_meas_temps):
#    sub = temp_subsets[set]
#    ax.plot(sub["M' (emu)"], sub["M'' (emu)"], label='{:4.2f}'.format(temps[set]))
#ax.legend()
#plt.show()