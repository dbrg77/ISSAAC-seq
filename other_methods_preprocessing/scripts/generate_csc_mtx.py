#!/usr/bin/env python

import gzip
from scipy.io import mmwrite
from scipy.sparse import lil_matrix

peak_idx = {}
with open('macs2_pk/aggregate_peaks_formatted.bed') as fh:
    for i,line in enumerate(fh):
        _, _, _, pn = line.strip().split('\t')
        peak_idx[pn] = i

num_peaks = len(peak_idx.keys())
num_cells = len(open('outs/raw_mtx/barcodes.tsv').readlines())

mtx = lil_matrix((num_peaks, num_cells), dtype=int)

with gzip.open('outs/peak_read_ov.tsv.gz', 'rt') as fh:
    for i,line in enumerate(fh):
        col_n = i
        count_info = line.strip().split('\t')[1]
        items = count_info.split(',')
        for pn_count in items:
            pn, count = pn_count.split(':')
            row_n = peak_idx[pn]
            mtx[row_n, col_n] = int(count)

mtx = mtx.tocsc()

mmwrite('outs/raw_mtx/matrix.mtx', mtx, field='integer')
