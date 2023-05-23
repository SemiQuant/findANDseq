#!/bin/python3
import sys
import pandas as pd
from Bio.Seq import Seq
# Read TSV file into DataFrame
df = pd.read_table(sys.argv[1])
df['genes'] = df.groupby('Seq')["Align"].transform("nunique")
df = pd.concat([df.reset_index(drop=True), df['ID'].str.split('_', expand=True)], axis=1)
opt_tm = int(sys.argv[3])

df = df.astype({3: float})
df.loc[df['Start'] >= 300, 'range'] = 2
df.loc[df['Start'] < 300, 'range'] = 3
df.loc[(df['Start'] > 300) & (df['Start'] < 600), 'range'] = 1 

df_best = df[(df["Start"] > 200) & df["Start"] < 300]
df_best = df_best[(df_best[3] >= opt_tm-1) & (df_best[3] >= opt_tm+3)]
df_best = df_best.sort_values(by=['genes'], ascending=False)

df = df.sort_values(by=["range", 'genes', 3], ascending=[True, False, True])
df = pd.concat([df_best, df])

# Calculate the median, mean, SD, IRQ and range of values
median = df[3].median()
mean = df[3].mean()
sd = df[3].std()
iqr1 = df[3].quantile(0.25)
iqr4 = df[3].quantile(0.75)
min = df[3].min()
max = df[3].max()
print(f"Median: {median}")
print(f"Mean: {mean}")
print(f"Standard Deviation: {sd}")
print(f"Interquartile Range: {iqr1}, {iqr4}")
print(f"Range of Values: {min}, {max}")

primers = []
while len(df.index) > 0:
    seq = df.iloc[0]["Seq"]
    primers.append(seq)
    added = df.loc[df['Seq'] == seq]["Align"]
    df = df[~df.Align.isin(added)]

# open file in write mode
with open(sys.argv[2], 'w') as fp:
    for seq in primers:
        seq_tailed = "CGTGTGCTCTTCCGATCTaa" + str(Seq(seq).reverse_complement()) + "CACTGCAGACCACTAA" #seq_tailed = "GCTCTTCCGATCTa" + seq.reverse_complement() + "CACTGCAGACCACTAA"
        fp.write("%s\n" % seq_tailed)



with open(sys.argv[4], 'w') as fp:
    count = 0
    for seq in primers:
        count += 1
        fp.write(">%s\n" % count)
        fp.write("%s\n" % str(Seq(seq).reverse_complement()))