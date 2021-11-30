import pandas as pd

def read_ppms_file(filename):
    
    with open(filename, 'r') as f:
        d = f.readlines()
    
    header_start, data_start = 0,0
    for i, line in enumerate(d):
        if '[Header]' in line:
            header_start = i+1
        elif '[Data]' in line:
            data_start = i+1
    
    header = d[header_start:data_start-1]
    header = [h.strip().split(',') for h in header if not h.startswith(';') and h.startswith('INFO')]
    header = {h[2]: h[1] for h in header}
    
    df = pd.read_csv(filename,
                    header=data_start,
                    engine='python')

    return header, df
    
