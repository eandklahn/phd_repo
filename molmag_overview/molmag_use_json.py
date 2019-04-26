import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ion_colors = {'Dy3+': 'c',
              'Fe1+': 'm',
              'Fe2+': 'm',
              'Co2+': 'm',
              'Mn': 'b',
              'None': '#A9A9A9'}
              
def plot_year_vs_ueff(c_l):
    
    # Determining what to plot
    lookup_these = [c_l.index(d) for d in c_l
                    if d.get('year')!=None
                    and d.get('ueff_cm')!=None]
    
    # Collecting data
    years = [c_l[index].get('year') for index in lookup_these]
    ueffs = [c_l[index].get('ueff_cm') for index in lookup_these]
    color = [ion_colors.get(str(c_l[index].get('ion'))) for index in lookup_these]

    # Making the plot
    f, ax = plt.subplots()
    ax.scatter(years, ueffs, color=color)
    
    # Y-axis manipulations
    ax.set_ylim((1, max(ueffs)+0.1*max(ueffs)))
    ax.set_ylabel(r'$U_{eff} [cm^{-1}]$')
    
    # X-axis manipulations
    ax.set_xlim(1993-1, max(years)+1)
    ax.set_xlabel('Year')
    
    # Constructing the legend
    
    legend_colors = list(set(color))
    legend_names = ['']*len(legend_colors)
    legend_handles = ['']*len(legend_colors)
    
    for key, val in ion_colors.items():
        if val in legend_colors:
            index = legend_colors.index(val)
            legend_names[index] += '{} '.format(key)
    
    legend_names = [l_n.strip() for l_n in legend_names]
    
    for i, c in enumerate(legend_colors):
        legend_handles[i] = Line2D([0], [0], color=c, label=legend_names[i], markersize=15)
    
    # Final plotting commands
    ax.legend(handles=legend_handles)
    ax.grid()
    plt.show()
    
    
if __name__ == '__main__':

    f = open('compounds.json', 'r')
    compounds = json.load(f)
    f.close()
    
    plot_year_vs_ueff(compounds)