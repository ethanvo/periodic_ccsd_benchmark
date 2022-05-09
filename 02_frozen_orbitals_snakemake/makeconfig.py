#!/usr/bin/env python
import pandas as pd
import itertools
data = {
'material': ['c', 'si', 'sic', 'bn', 'bp', 'aln', 'alp', 'mgo', 'mgs', 'lih', 'lif', 'licl'],
'basis': ['ccpvdz'],
'kdensity': ['1', '2', '3', '4']
}

rows = ['_'.join(element) for element in itertools.product(*data.values())]
print(rows)
