#!/bin/python3
import sys
import pandas as pd
from Bio.Seq import Seq

# Read TSV file into DataFrame
df = pd.read_table(sys.argv[1])
df['genes'] = df.groupby('Seq')["Align"].transform("nunique")
df = df.sort_values(by=['genes'], ascending=False)
primers = []
while len(df.index) > 0:
  seq = df.iloc[0]["Seq"]
primers.append(seq)
added = df.loc[df['Seq'] == seq]["Align"]
df = df[~df.Align.isin(added)]


# open file in write mode
with open(sys.argv[1], 'w') as fp:
    for seq in primers:
    	seq_tailed = "GCTCTTCCGATCTa" + seq.reverse_complement() + "CACTGCAGACCACTAA"
        fp.write("%s\n" % seq_tailed)


