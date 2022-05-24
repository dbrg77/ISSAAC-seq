#!/usr/bin/env python

import gzip
import sys

i7 = sys.argv[1]
i5 = sys.argv[2]

i7f = gzip.open(i7, 'rt')
i5f = gzip.open(i5, 'rt')

for i7n, i7s, _, i7q, i5n, i5s, _, i5q in zip(i7f, i7f, i7f, i7f, i5f, i5f, i5f, i5f):
    n1 = i7n.split(' ')[0]
    n2 = i5n.split(' ')[0]
    assert n1 == n2
    print(i7n.strip())
    print(i7s.strip() + i5s.strip())
    print('+')
    print(i7q.strip() + i5q.strip())

i7f.close()
i5f.close()
