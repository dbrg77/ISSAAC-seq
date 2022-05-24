#!/usr/bin/env python

import gzip
import sys

i7 = sys.argv[1]
i5 = sys.argv[2]
umi = sys.argv[3]

i7f = gzip.open(i7, 'rt')
i5f = gzip.open(i5, 'rt')
umif = gzip.open(umi, 'rt')

for i7n, i7s, _, i7q, i5n, i5s, _, i5q, un, us, _, uq in zip(i7f, i7f, i7f, i7f, i5f, i5f, i5f, i5f, umif, umif, umif, umif):
        n1 = i7n.split(' ')[0]
        n2 = i5n.split(' ')[0]
        n3 = un.split(' ')[0]
        assert n1 == n2 and n1 == n3
        print(i7n.strip())
        print(i7s.strip() + i5s.strip() + us[:10])
        print('+')
        print(i7q.strip() + i5q.strip() + uq[:10])

i7f.close()
i5f.close()
umif.close()
