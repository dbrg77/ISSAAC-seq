#!/usr/bin/env python

import sys
import gzip

r2 = gzip.open(sys.argv[1], mode='rt')
i2 = gzip.open(sys.argv[2], mode='rt')

for un, us, _, uq, cbn, cbs, _, cbq in zip(r2, r2, r2, r2, i2, i2, i2, i2):
    n1 = un.split(' ')[0]
    n2 = cbn.split(' ')[0]
    assert n1 == n2
    print(un.strip())
    print(cbs.strip() + us[:10])
    print('+')
    print(cbq.strip() + uq[:10])
