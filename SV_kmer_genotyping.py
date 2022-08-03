import sys
import os
import pandas as pd
import fnmatch
import os
import warnings
import shutil

warnings.filterwarnings("ignore")

folder=sys.argv[1]
inputfile_SV_infor=sys.argv[2]
inputfile_pos=sys.argv[3]
inputfile_isolate=sys.argv[4]
outputfile_genotype=sys.argv[5]

os.chdir(folder)




df_isolate=pd.read_csv(inputfile_isolate,sep="\t",header=None)
df_isolate.columns=['isolate']
isolate_list=df_isolate['isolate'].tolist()


for i in isolate_list:
    outputfile_kmer=folder+"/kmer."+i+".count"
    kmer_dir=folder+"/kmer_counts_in_each_isolate/"+i
    os.chdir(kmer_dir)
    df_append_iso = pd.DataFrame()
    for filename in os.listdir(kmer_dir):
        if fnmatch.fnmatch(filename,"*txt"):
            df_iso_SV=pd.read_csv(filename,delim_whitespace=True,header=None)
            df_iso_SV.columns=['kmer','count']
            df_iso_SV['SV']=filename
            df_append_iso=df_append_iso.append(df_iso_SV)

    df_iso_all= pd.DataFrame(df_append_iso)
    df_iso_all['SV_ID']=df_iso_all['SV'].replace('.txt', '', regex=True)
    df_iso_all = df_iso_all.drop(["SV"], axis=1)
    df_iso_all=df_iso_all[df_iso_all['count']>0]
    df_iso_all_filter_size=df_iso_all.groupby(['SV_ID']).size().reset_index(name='present_kmer_counts')
    df_iso_all_filter_size=df_iso_all_filter_size.rename(columns = {'present_kmer_counts':i})
    df_iso_all_filter_size.to_csv(outputfile_kmer,sep="\t",index=False)











os.chdir(folder)
df_merge = pd.DataFrame(columns=['SV_ID'])
for i in isolate_list:
    df=pd.read_csv("kmer."+i+".count",sep="\t")
    df_merge=df_merge.merge(df, on='SV_ID', how='outer')
df_merge.fillna(0, inplace=True)

df_pos=pd.read_csv(inputfile_pos,sep="\t")
df_pos['SV_ID']=df_pos['ID'].str.split("::").str[0]

df_pos['group']=df_pos['SV_ID'].str.split("_").str[0]
df_pos['group'] = df_pos['group'].str.extract('(\d+)', expand=False).astype(int)

df_pos['pos']=df_pos['SV_ID'].str.split("_").str[1]
df_pos['pos']=df_pos['pos'].str.split("-").str[0].astype(int)

df_pos=df_pos[['SV_ID','group','pos']]
df_pos=df_pos.drop_duplicates(keep='first')

df_pos=df_pos.groupby('SV_ID').head(1).reset_index(drop=True)

df_merge_all=pd.merge(df_pos,df_merge,how="inner",on="SV_ID")
df_merge_all = df_merge_all.sort_values(['group','pos'], ascending=[True,True])







df_merge_all_sort = df_merge_all.sort_values(['group','pos'], ascending=[True,True])


df_merge_all_sort['cat']=df_merge_all_sort['SV_ID'].str.split("_").str[2]
df_merge_all_sort['SV_ID_p1']=df_merge_all_sort['SV_ID'].str.split("_").str[0]
df_merge_all_sort['SV_ID_p2']=df_merge_all_sort['SV_ID'].str.split("_").str[1]
df_merge_all_sort['SV_ID_p3']=df_merge_all_sort['SV_ID'].str.split("_").str[3]
df_merge_all_sort['SV_ID']=df_merge_all_sort['SV_ID_p1']+"_"+df_merge_all_sort['SV_ID_p2']+"_"+df_merge_all_sort['SV_ID_p3']
df_merge_all_sort = df_merge_all_sort.drop(["SV_ID_p1", "SV_ID_p2",'SV_ID_p3'], axis=1)


df_merge_all_sort_presence=df_merge_all_sort[df_merge_all_sort['cat']=='presence']
for i in isolate_list:
    df_merge_all_sort_presence=df_merge_all_sort_presence.rename(columns = {i:i+"_presence"})
df_merge_all_sort_presence = df_merge_all_sort_presence.drop(["cat"], axis=1)

df_merge_all_sort_absence=df_merge_all_sort[df_merge_all_sort['cat']=='absence']
for i in isolate_list:
    df_merge_all_sort_absence=df_merge_all_sort_absence.rename(columns = {i:i+"_absence"})
df_merge_all_sort_absence = df_merge_all_sort_absence.drop(["cat"], axis=1)

df_concat_merge=pd.merge(df_merge_all_sort_presence,df_merge_all_sort_absence,how="inner",on=['SV_ID','group','pos'])
df_concat_merge['SV_ID_p1']=df_concat_merge['SV_ID'].str.split("_").str[0]
df_concat_merge['SV_ID_p2']=df_concat_merge['SV_ID'].str.split("_").str[1]
df_concat_merge['SV_ID']=df_concat_merge['SV_ID_p1']+"_"+df_concat_merge['SV_ID_p2']
df_concat_merge = df_concat_merge.drop(["SV_ID_p1", "SV_ID_p2"], axis=1)


df_SV_infor=pd.read_csv(inputfile_SV_infor,sep="\t")
df_SV_infor=df_SV_infor[['name','SV_type','size']]
df_SV_infor.columns=['name','SV_type','SV_size']
df_SV_infor=df_SV_infor.rename(columns = {'name':'SV_ID'})
df_SV_infor=df_SV_infor.drop_duplicates(keep='first')
df_SV_infor=df_SV_infor.groupby('SV_ID').head(1).reset_index(drop=True)


df_concat_merge_final=pd.merge(df_SV_infor,df_concat_merge,how='inner',on='SV_ID')
df_concat_merge_final = df_concat_merge_final.sort_values(['group','pos'], ascending=[True,True])

for i in isolate_list:
    df_concat_merge_final.loc[((df_concat_merge_final[i+'_presence'] >= 15) & (df_concat_merge_final[i+'_absence'] == 0) &
        (df_concat_merge_final['SV_type'].str.contains('insertion|deletion'))) , [i+'_genotype']] = '1/1'
    
    df_concat_merge_final.loc[((df_concat_merge_final[i+'_presence'] >= 15) & (df_concat_merge_final[i+'_absence'] >= 5) &
        (df_concat_merge_final['SV_type'].str.contains('insertion|deletion'))) , [i+'_genotype']] = '1/0'
    
    df_concat_merge_final.loc[((df_concat_merge_final[i+'_presence'] < 5 ) & (df_concat_merge_final[i+'_absence'] >= 5) &
        (df_concat_merge_final['SV_type'].str.contains('insertion|deletion'))) , [i+'_genotype']] = '0/0'
    
    df_concat_merge_final = df_concat_merge_final.drop([i+'_presence',i+'_absence'], axis=1)
    df_concat_merge_final=df_concat_merge_final.rename(columns = {i+"_genotype":i})

df_concat_merge_final=df_concat_merge_final.fillna('./.')
df_concat_merge_final.to_csv(outputfile_genotype,sep="\t",index=False)


files = os.listdir(folder)
for item in files:
    if item.startswith(("kmer.")):
        os.remove(os.path.join(folder, item))