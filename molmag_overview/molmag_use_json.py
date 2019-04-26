import json
import matplotlib.pyplot as plt

def plot_year_vs_ueff(c_l):

    y = []
    u = []
    for d in c_l:
        y.append(d.get('year'))
        u.append(d.get('ueff_cm'))
    
    
    
    
    
    
if __name__ == '__main__':

    f = open('compounds.json', 'r')
    compounds = json.load(f)
    f.close()
    
    plot_year_vs_ueff(compounds)