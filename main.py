import pandas as pd

data = pd.read_csv("mutect.snv.res", sep="\t", low_memory=False)

def vaf(sample):
    variant = []
    for i in sample:
        x = i.split('|')
        variant.append(float(x[0]))
    return variant


CRCTumorOA1maf = vaf(data['CRCTumorOA1maf'])
CRCTumorOA2maf = vaf(data['CRCTumorOA2maf'])
CRCTumorOA3maf = vaf(data['CRCTumorOA3maf'])
CRCTumorOA4maf = vaf(data['CRCTumorOA4maf'])
CRCTumorOB1maf = vaf(data['CRCTumorOB1maf'])
CRCTumorUA_Lmaf = vaf(data['CRCTumorUA_Lmaf'])
CRCTumorUB_Lmaf = vaf(data['CRCTumorUB_Lmaf'])
OAmaf = vaf(data['OAmaf'])
OBmaf = vaf(data['OBmaf'])
UAmaf = vaf(data['UAmaf'])
UBmaf = vaf(data['UBmaf'])

samplesVAF = [CRCTumorOA1maf, CRCTumorOA2maf, CRCTumorOA3maf, CRCTumorOA4maf, CRCTumorOB1maf, CRCTumorUA_Lmaf,
              CRCTumorUB_Lmaf, OAmaf, UAmaf, UBmaf]

df = pd.DataFrame(samplesVAF)
df = df.transpose()
df.columns = ['CRCTumorOA1maf', 'CRCTumorOA2maf', 'CRCTumorOA3maf', 'CRCTumorOA4maf', 'CRCTumorOB1maf', 'CRCTumorUA_Lmaf',
              'CRCTumorUB_Lmaf', 'OAmaf', 'UAmaf', 'UBmaf']

df.insert(loc=0, column='#chr', value=data["#chr"])
df.insert(loc=1, column='pos', value=data["pos"])
df.insert(loc=2, column='ref', value=data["ref"])
df.insert(loc=3, column='alt', value=data["alt"])

df.to_csv('VAF_Data.csv', index=False)
