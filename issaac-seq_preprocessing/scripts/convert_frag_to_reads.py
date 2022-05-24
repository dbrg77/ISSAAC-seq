#!/usr/bin/env python

import gzip
from collections import defaultdict

cell_rc = defaultdict(int)

frags = gzip.open('outs/fragments.tsv.gz', 'rt')
cb_summary = open('outs/per_cell_barcode_total_fragment_count.tsv', 'w')

for i,j in enumerate(frags):
    chrom, start, end, cb, _, = j.strip().split('\t')
    cell_rc[cb] += 1
    print('%s\t%s\t%s\t%s\t.\t+' % (chrom, int(start)-100, int(start)+100, cb))
    print('%s\t%s\t%s\t%s\t.\t-' % (chrom, int(end)-100, int(end)+100, cb))

frags.close()

for cb, count in cell_rc.items():
    cb_summary.write('%s\t%s\n' % (cb, count))

cb_summary.close()
