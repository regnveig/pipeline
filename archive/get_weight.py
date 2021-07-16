import cooler
import cooltools
import sys
import pandas as pd

chr = str(sys.argv[2]).split('chr')[1]
path = sys.argv[3]
table = cooler.Cooler(sys.argv[1]).bins().fetch('chr'+chr)
weights = table['weight'][:]
new_name = 'vector  c-tale_normalization  '+chr+'    5000 BP'
weights.name = new_name
weights.to_csv(path+'/vector.txt', index=False)
