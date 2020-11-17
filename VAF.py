import pandas as pd

data = pd.read_csv("venv/bin/Input/mutect.snv.res", sep="\t", low_memory=False)

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
              CRCTumorUB_Lmaf, OAmaf, OBmaf, UAmaf, UBmaf]

df = pd.DataFrame(samplesVAF)
df = df.transpose()
df.columns = ['CRCTumorOA1maf', 'CRCTumorOA2maf', 'CRCTumorOA3maf', 'CRCTumorOA4maf', 'CRCTumorOB1maf', 'CRCTumorUA_Lmaf',
              'CRCTumorUB_Lmaf', 'OAmaf','OBmaf', 'UAmaf', 'UBmaf']

df.insert(loc=0, column='chr', value=data["#chr"])
df.insert(loc=1, column='pos', value=data["pos"])
df.insert(loc=3, column='CRCTumorOA1d', value=data["CRCTumorOA1d"])
df.insert(loc=4, column='CRCTumorOA2d', value=data["CRCTumorOA2d"])
df.insert(loc=5, column='CRCTumorOA3d', value=data["CRCTumorOA3d"])
df.insert(loc=6, column='CRCTumorOA4d', value=data["CRCTumorOA4d"])
df.insert(loc=7, column='CRCTumorOB1d', value=data["CRCTumorOB1d"])
df.insert(loc=8, column='CRCTumorUA_Ld', value=data["CRCTumorUA_Ld"])
df.insert(loc=9, column='CRCTumorUB_Ld', value=data["CRCTumorUB_Ld"])
df.insert(loc=10, column='OAd', value=data["OAd"])
df.insert(loc=11, column='OBd', value=data["OBd"])
df.insert(loc=12, column='UAd', value=data["UAd"])
df.insert(loc=13, column='UBd', value=data["UBd"])


df.to_csv('VAF_Data.txt', index=False)
