import pandas as pd

data = pd.read_csv("mutect.snv.res.filtered.classified.founds.nopara.somatic", sep="\t", low_memory=False)

def input(sample):
    
    vaf = []
    endsratio = []
    cmedian = []
    strandFisherP = []
    badQualFrac = []
    logLikelihood = []
    indelsMean = []

    for i in sample:
        y = i.replace(',', '|')
        x = y.split('|')
        if len(x) != 1:
            vaf.append(float(x[0]))
            endsratio.append(float(x[1]))
            cmedian.append(float(x[3]))
            strandFisherP.append(float(x[6]))
            badQualFrac.append(float(x[7]))
            logLikelihood.append(float(x[8]))
            indelsMean.append(float(x[13]))
        else:
            vaf.append(0)
            endsratio.append(0)
            cmedian.append(0)
            strandFisherP.append(0)
            badQualFrac.append(0)
            logLikelihood.append(0)
            indelsMean.append(0)

    return [vaf, endsratio, cmedian, strandFisherP, badQualFrac, logLikelihood, indelsMean]


OAmaf = input(data['OAmaf'])
OBmaf = input(data['OBmaf'])
UAmaf = input(data['UAmaf'])
UBmaf = input(data['UBmaf'])


CRCTumorOA1maf = input(data['CRCTumorOA1maf'])
CRCTumorOA2maf = input(data['CRCTumorOA2maf'])
CRCTumorOA3maf = input(data['CRCTumorOA3maf'])
CRCTumorOA4maf = input(data['CRCTumorOA4maf'])
CRCTumorOB1maf = input(data['CRCTumorOB1maf'])
CRCTumorUA_Lmaf = input(data['CRCTumorUA_Lmaf'])
CRCTumorUB_Lmaf = input(data['CRCTumorUB_Lmaf'])


df_1 = pd.DataFrame(CRCTumorOA1maf)
df_1 = df_1.transpose()
df_1.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_1.insert(loc=0, column='chr', value=data["chr"])
df_1.insert(loc=1, column='pos', value=data["pos"])
df_1.to_csv('CRCTumorOA1maf.txt', index=False)

df_2 = pd.DataFrame(CRCTumorOA2maf)
df_2 = df_2.transpose()
df_2.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_2.insert(loc=0, column='chr', value=data["chr"])
df_2.insert(loc=1, column='pos', value=data["pos"])
df_2.to_csv('CRCTumorOA2maf.txt', index=False)

df_3 = pd.DataFrame(CRCTumorOA3maf)
df_3 = df_3.transpose()
df_3.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_3.insert(loc=0, column='chr', value=data["chr"])
df_3.insert(loc=1, column='pos', value=data["pos"])
df_3.to_csv('CRCTumorOA3maf.txt', index=False)

df_4 = pd.DataFrame(CRCTumorOA4maf)
df_4 = df_4.transpose()
df_4.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_4.insert(loc=0, column='chr', value=data["chr"])
df_4.insert(loc=1, column='pos', value=data["pos"])
df_4.to_csv('CRCTumorOA4maf.txt', index=False)

df_5 = pd.DataFrame(OAmaf)
df_5 = df_5.transpose()
df_5.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_5.insert(loc=0, column='chr', value=data["chr"])
df_5.insert(loc=1, column='pos', value=data["pos"])
df_5.to_csv('OAmaf.txt', index=False)

df_6 = pd.DataFrame(CRCTumorOB1maf)
df_6 = df_6.transpose()
df_6.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_6.insert(loc=0, column='chr', value=data["chr"])
df_6.insert(loc=1, column='pos', value=data["pos"])
df_6.to_csv('CRCTumorOB1maf.txt', index=False)

df_7 = pd.DataFrame(OBmaf)
df_7 = df_7.transpose()
df_7.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_7.insert(loc=0, column='chr', value=data["chr"])
df_7.insert(loc=1, column='pos', value=data["pos"])
df_7.to_csv('OBmaf.txt', index=False)

df_8 = pd.DataFrame(UAmaf)
df_8 = df_8.transpose()
df_8.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_8.insert(loc=0, column='chr', value=data["chr"])
df_8.insert(loc=1, column='pos', value=data["pos"])
df_8.to_csv('UAmaf.txt', index=False)

df_9 = pd.DataFrame(UBmaf)
df_9 = df_9.transpose()
df_9.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_9.insert(loc=0, column='chr', value=data["chr"])
df_9.insert(loc=1, column='pos', value=data["pos"])
df_9.to_csv('UBmaf.txt', index=False)

df_10 = pd.DataFrame(CRCTumorUB_Lmaf)
df_10 = df_10.transpose()
df_10.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_10.insert(loc=0, column='chr', value=data["chr"])
df_10.insert(loc=1, column='pos', value=data["pos"])
df_10.to_csv('CRCTumorUB_Lmaf.txt', index=False)

df_11 = pd.DataFrame(CRCTumorUA_Lmaf)
df_11 = df_11.transpose()
df_11.columns = ['vaf','endsratio','cmedian','strandFisherP','badQualFrac','logLikelihood','indelsMean']
df_11.insert(loc=0, column='chr', value=data["chr"])
df_11.insert(loc=1, column='pos', value=data["pos"])
df_11.to_csv('CRCTumorUA_Lmaf.txt', index=False)


