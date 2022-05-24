#!/usr/bin/env python

import sys
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csc_matrix

if sys.argv[1] == '-':
    df = pd.read_csv(sys.stdin, sep='\t', index_col=0)
else:
    df = pd.read_csv(sys.argv[1], sep='\t', index_col=0)

mtx = csc_matrix(df)
mmwrite(sys.argv[2] + '/matrix.mtx', mtx, field='integer')

with open(sys.argv[2] + '/features.tsv', 'w') as f:
    f.write('\n'.join(df.index.values))
    f.write('\n')

with open(sys.argv[2] + '/barcodes.tsv', 'w') as f:
    f.write('\n'.join(df.columns.values))
    f.write('\n')
