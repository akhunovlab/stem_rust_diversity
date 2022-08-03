import pandas as pd
import numpy as np
import os
import sys
import pybedtools
import fnmatch
import warnings
import shutil

warnings.filterwarnings("ignore")


folder=sys.argv[1]
inputfile_haplotype=sys.argv[2]
hap1_ID=sys.argv[3]
hap2_ID=sys.argv[4]
hap3_ID=sys.argv[5]
hap4_ID=sys.argv[6]
outputfile_SV_info=sys.argv[7]

os.chdir(folder)


df_haplotype=pd.read_csv(inputfile_haplotype,sep="\t",header=None)
df_haplotype.columns=['hap1','hap2','hap3','hap4']

haplotype_list=df_haplotype[['hap1','hap2','hap3','hap4']].values.tolist()

for i in haplotype_list:
    hap1=i[0]
    hap2=i[1]
    hap3=i[2]
    hap4=i[3]
    
    input_coords_1_2="ref_"+hap1+".que_"+hap2+".filter.coords"
    input_coords_1_3="ref_"+hap1+".que_"+hap3+".filter.coords"
    input_coords_1_4="ref_"+hap1+".que_"+hap4+".filter.coords"
    input_coords_2_3="ref_"+hap2+".que_"+hap3+".filter.coords"
    input_coords_2_4="ref_"+hap2+".que_"+hap4+".filter.coords"

    output_for_genotype_SV_coords_1_2="output_for_genotype_SV_coords."+hap1+"."+hap2_ID+".txt"
    output_for_genotype_SV_coords_1_3="output_for_genotype_SV_coords."+hap1+"."+hap3_ID+".txt"
    output_for_genotype_SV_coords_1_4="output_for_genotype_SV_coords."+hap1+"."+hap4_ID+".txt"


    ###hap_1_2
    df_coords_1_2=pd.read_csv(input_coords_1_2,sep="\t",header=None)
    df_coords_1_2.columns=['que','a','b','c','d','ref','que_start','que_end','ref_start','ref_end','e','f','g','h','i','j','k','strand','m','n','o']
    df_coords_1_2=df_coords_1_2[['que','ref','que_start','que_end','ref_start','ref_end','strand']]
    df_coords_1_2=df_coords_1_2[df_coords_1_2['que_start']<=df_coords_1_2['que_end']]

    df_coords_1_2['que_end_previous_1'] = df_coords_1_2['que_end'].shift(1,fill_value=0)
    df_coords_1_2['que_end_previous_2'] = df_coords_1_2['que_end'].shift(2,fill_value=0)
    df_coords_1_2['que_end_previous_3'] = df_coords_1_2['que_end'].shift(3,fill_value=0)
    df_coords_1_2['que_end_previous_4'] = df_coords_1_2['que_end'].shift(4,fill_value=0)

    df_coords_1_2['que_start_next_1'] = df_coords_1_2['que_start'].shift(-1,fill_value=0)
    df_coords_1_2['que_start_next_2'] = df_coords_1_2['que_start'].shift(-2,fill_value=0)
    df_coords_1_2['que_start_next_3'] = df_coords_1_2['que_start'].shift(-3,fill_value=0)
    df_coords_1_2['que_start_next_4'] = df_coords_1_2['que_start'].shift(-4,fill_value=0)

    df_coords_1_2=df_coords_1_2[
       (df_coords_1_2['que_start']-df_coords_1_2['que_end_previous_1'] > -10  ) &
        (df_coords_1_2['que_start']-df_coords_1_2['que_end_previous_2'] > -10  ) &
        (df_coords_1_2['que_start']-df_coords_1_2['que_end_previous_3'] > -10  ) &
        (df_coords_1_2['que_start']-df_coords_1_2['que_end_previous_4'] > -10  ) &

        (df_coords_1_2['que_start_next_1']-df_coords_1_2['que_end'] > -10  ) &
        (df_coords_1_2['que_start_next_2']-df_coords_1_2['que_end'] > -10  ) &
        (df_coords_1_2['que_start_next_3']-df_coords_1_2['que_end'] > -10  ) &
        (df_coords_1_2['que_start_next_4']-df_coords_1_2['que_end'] > -10  ) &

        (df_coords_1_2['que_end_previous_1'] != 0) &
        (df_coords_1_2['que_end_previous_2'] != 0) &
        (df_coords_1_2['que_end_previous_3'] != 0) &
        (df_coords_1_2['que_end_previous_4'] != 0) &

        (df_coords_1_2['que_start_next_1'] != 0) &
        (df_coords_1_2['que_start_next_2'] != 0) &
        (df_coords_1_2['que_start_next_3'] != 0) &
        (df_coords_1_2['que_start_next_4'] != 0)
        ]
    
    df_coords_1_2_align=df_coords_1_2[['que','que_start','que_end']]
    df_coords_1_2_align.loc[((df_coords_1_2_align['que_start'] <= df_coords_1_2_align['que_end'])), 'que_start_final'] = df_coords_1_2_align['que_start']
    df_coords_1_2_align.loc[((df_coords_1_2_align['que_start'] <= df_coords_1_2_align['que_end'])), 'que_end_final'] = df_coords_1_2_align['que_end']
    df_coords_1_2_align.loc[((df_coords_1_2_align['que_start'] > df_coords_1_2_align['que_end'])), 'que_start_final'] = df_coords_1_2_align['que_end']
    df_coords_1_2_align.loc[((df_coords_1_2_align['que_start'] > df_coords_1_2_align['que_end'])), 'que_end_final'] = df_coords_1_2_align['que_start']
    df_coords_1_2_align=df_coords_1_2_align[['que','que_start_final','que_end_final']]
    df_coords_1_2_align=pybedtools.BedTool.from_dataframe(df_coords_1_2_align)

    df_coords_1_2=df_coords_1_2[['que','ref','que_start','que_end', 'ref_start','ref_end', 'strand']]
    df_coords_1_2['que_end_previous'] = df_coords_1_2['que_end'].shift(1,fill_value=0)
    df_coords_1_2['ref_end_previous'] = df_coords_1_2['ref_end'].shift(1,fill_value=0)

    df_coords_1_2=df_coords_1_2[
        (df_coords_1_2['ref_end_previous'] != 0)
        ]

    df_coords_1_2['gap_que']=df_coords_1_2['que_start']-df_coords_1_2['que_end_previous']
    df_coords_1_2['gap_ref']=df_coords_1_2['ref_start']-df_coords_1_2['ref_end_previous']
    df_coords_1_2=df_coords_1_2[(abs(df_coords_1_2['gap_que'])<=30000) |(abs(df_coords_1_2['gap_ref'])<=30000) ]
    
    df_coords_1_2=df_coords_1_2[(abs(df_coords_1_2['gap_que']-df_coords_1_2['gap_ref'])>10000) & (abs(df_coords_1_2['gap_que']-df_coords_1_2['gap_ref'])<=1500000)]
    df_coords_1_2['gap']=abs(df_coords_1_2['gap_que']-df_coords_1_2['gap_ref'])

    df_coords_1_2.loc[((df_coords_1_2['ref_end_previous'] <= df_coords_1_2['ref_start'])), 'ref_end_previous_final'] = df_coords_1_2['ref_end_previous']
    df_coords_1_2.loc[((df_coords_1_2['ref_end_previous'] <= df_coords_1_2['ref_start'])), 'ref_start_final'] = df_coords_1_2['ref_start']
    df_coords_1_2.loc[((df_coords_1_2['ref_end_previous'] > df_coords_1_2['ref_start'])), 'ref_end_previous_final'] = df_coords_1_2['ref_start']
    df_coords_1_2.loc[((df_coords_1_2['ref_end_previous'] > df_coords_1_2['ref_start'])), 'ref_start_final'] = df_coords_1_2['ref_end_previous']

    df_coords_1_2.loc[((df_coords_1_2['que_end_previous'] <= df_coords_1_2['que_start'])), 'que_end_previous_final'] = df_coords_1_2['que_end_previous']
    df_coords_1_2.loc[((df_coords_1_2['que_end_previous'] <= df_coords_1_2['que_start'])), 'que_start_final'] = df_coords_1_2['que_start']
    df_coords_1_2.loc[((df_coords_1_2['que_end_previous'] > df_coords_1_2['que_start'])), 'que_end_previous_final'] = df_coords_1_2['que_start']
    df_coords_1_2.loc[((df_coords_1_2['que_end_previous'] > df_coords_1_2['que_start'])), 'que_start_final'] = df_coords_1_2['que_end_previous']

    df_coords_1_2['ref_end_previous_final']=df_coords_1_2['ref_end_previous_final'].astype(int)
    df_coords_1_2['ref_start_final']=df_coords_1_2['ref_start_final'].astype(int)
    df_coords_1_2['que_end_previous_final']=df_coords_1_2['que_end_previous_final'].astype(int)
    df_coords_1_2['que_start_final']=df_coords_1_2['que_start_final'].astype(int)
    
    df_coords_1_2_SV=df_coords_1_2[['que','que_end_previous_final','que_start_final']]
    df_coords_1_2_SV=pybedtools.BedTool.from_dataframe(df_coords_1_2_SV)
    df_coords_1_2_SV_final=df_coords_1_2_SV.subtract(df_coords_1_2_align,A=True)
    
    df_coords_1_2_SV_final = df_coords_1_2_SV_final.to_dataframe()
    df_coords_1_2_SV_final.columns=['que','que_end_previous_final','que_start_final']
    df_coords_1_2_SV_final_merge=pd.merge(df_coords_1_2_SV_final,df_coords_1_2,how="inner",on=['que','que_end_previous_final','que_start_final'])
    
    df_coords_1_2_SV_final_merge_output=df_coords_1_2_SV_final_merge[['ref','ref_end_previous_final','ref_start_final','que','que_end_previous_final','que_start_final']]
    df_coords_1_2_SV_final_merge_output.to_csv(output_for_genotype_SV_coords_1_2,sep="\t",index=False)









    ###hap_1_3
    df_coords_1_3=pd.read_csv(input_coords_1_3,sep="\t",header=None)
    df_coords_1_3.columns=['que','a','b','c','d','ref','que_start','que_end','ref_start','ref_end','e','f','g','h','i','j','k','strand','m','n','o']
    df_coords_1_3=df_coords_1_3[['que','ref','que_start','que_end','ref_start','ref_end','strand']]
    df_coords_1_3=df_coords_1_3[df_coords_1_3['que_start']<=df_coords_1_3['que_end']]

    df_coords_1_3['que_end_previous_1'] = df_coords_1_3['que_end'].shift(1,fill_value=0)
    df_coords_1_3['que_end_previous_2'] = df_coords_1_3['que_end'].shift(2,fill_value=0)
    df_coords_1_3['que_end_previous_3'] = df_coords_1_3['que_end'].shift(3,fill_value=0)
    df_coords_1_3['que_end_previous_4'] = df_coords_1_3['que_end'].shift(4,fill_value=0)

    df_coords_1_3['que_start_next_1'] = df_coords_1_3['que_start'].shift(-1,fill_value=0)
    df_coords_1_3['que_start_next_2'] = df_coords_1_3['que_start'].shift(-2,fill_value=0)
    df_coords_1_3['que_start_next_3'] = df_coords_1_3['que_start'].shift(-3,fill_value=0)
    df_coords_1_3['que_start_next_4'] = df_coords_1_3['que_start'].shift(-4,fill_value=0)

    df_coords_1_3=df_coords_1_3[
       (df_coords_1_3['que_start']-df_coords_1_3['que_end_previous_1'] > -10  ) &
        (df_coords_1_3['que_start']-df_coords_1_3['que_end_previous_2'] > -10  ) &
        (df_coords_1_3['que_start']-df_coords_1_3['que_end_previous_3'] > -10  ) &
        (df_coords_1_3['que_start']-df_coords_1_3['que_end_previous_4'] > -10  ) &

        (df_coords_1_3['que_start_next_1']-df_coords_1_3['que_end'] > -10  ) &
        (df_coords_1_3['que_start_next_2']-df_coords_1_3['que_end'] > -10  ) &
        (df_coords_1_3['que_start_next_3']-df_coords_1_3['que_end'] > -10  ) &
        (df_coords_1_3['que_start_next_4']-df_coords_1_3['que_end'] > -10  ) &

        (df_coords_1_3['que_end_previous_1'] != 0) &
        (df_coords_1_3['que_end_previous_2'] != 0) &
        (df_coords_1_3['que_end_previous_3'] != 0) &
        (df_coords_1_3['que_end_previous_4'] != 0) &

        (df_coords_1_3['que_start_next_1'] != 0) &
        (df_coords_1_3['que_start_next_2'] != 0) &
        (df_coords_1_3['que_start_next_3'] != 0) &
        (df_coords_1_3['que_start_next_4'] != 0)
        ]
    
    df_coords_1_3_align=df_coords_1_3[['que','que_start','que_end']]
    df_coords_1_3_align.loc[((df_coords_1_3_align['que_start'] <= df_coords_1_3_align['que_end'])), 'que_start_final'] = df_coords_1_3_align['que_start']
    df_coords_1_3_align.loc[((df_coords_1_3_align['que_start'] <= df_coords_1_3_align['que_end'])), 'que_end_final'] = df_coords_1_3_align['que_end']
    df_coords_1_3_align.loc[((df_coords_1_3_align['que_start'] > df_coords_1_3_align['que_end'])), 'que_start_final'] = df_coords_1_3_align['que_end']
    df_coords_1_3_align.loc[((df_coords_1_3_align['que_start'] > df_coords_1_3_align['que_end'])), 'que_end_final'] = df_coords_1_3_align['que_start']
    df_coords_1_3_align=df_coords_1_3_align[['que','que_start_final','que_end_final']]
    df_coords_1_3_align=pybedtools.BedTool.from_dataframe(df_coords_1_3_align)

    df_coords_1_3=df_coords_1_3[['que','ref','que_start','que_end', 'ref_start','ref_end', 'strand']]
    df_coords_1_3['que_end_previous'] = df_coords_1_3['que_end'].shift(1,fill_value=0)
    df_coords_1_3['ref_end_previous'] = df_coords_1_3['ref_end'].shift(1,fill_value=0)

    df_coords_1_3=df_coords_1_3[
        (df_coords_1_3['ref_end_previous'] != 0)
        ]

    df_coords_1_3['gap_que']=df_coords_1_3['que_start']-df_coords_1_3['que_end_previous']
    df_coords_1_3['gap_ref']=df_coords_1_3['ref_start']-df_coords_1_3['ref_end_previous']
    df_coords_1_3=df_coords_1_3[(abs(df_coords_1_3['gap_que'])<=30000) |(abs(df_coords_1_3['gap_ref'])<=30000) ]
    
    df_coords_1_3=df_coords_1_3[ (abs(df_coords_1_3['gap_que']-df_coords_1_3['gap_ref'])>10000) & (abs(df_coords_1_3['gap_que']-df_coords_1_3['gap_ref'])<=1500000)]
    df_coords_1_3['gap']=abs(df_coords_1_3['gap_que']-df_coords_1_3['gap_ref'])

    df_coords_1_3.loc[((df_coords_1_3['ref_end_previous'] <= df_coords_1_3['ref_start'])), 'ref_end_previous_final'] = df_coords_1_3['ref_end_previous']
    df_coords_1_3.loc[((df_coords_1_3['ref_end_previous'] <= df_coords_1_3['ref_start'])), 'ref_start_final'] = df_coords_1_3['ref_start']
    df_coords_1_3.loc[((df_coords_1_3['ref_end_previous'] > df_coords_1_3['ref_start'])), 'ref_end_previous_final'] = df_coords_1_3['ref_start']
    df_coords_1_3.loc[((df_coords_1_3['ref_end_previous'] > df_coords_1_3['ref_start'])), 'ref_start_final'] = df_coords_1_3['ref_end_previous']

    df_coords_1_3.loc[((df_coords_1_3['que_end_previous'] <= df_coords_1_3['que_start'])), 'que_end_previous_final'] = df_coords_1_3['que_end_previous']
    df_coords_1_3.loc[((df_coords_1_3['que_end_previous'] <= df_coords_1_3['que_start'])), 'que_start_final'] = df_coords_1_3['que_start']
    df_coords_1_3.loc[((df_coords_1_3['que_end_previous'] > df_coords_1_3['que_start'])), 'que_end_previous_final'] = df_coords_1_3['que_start']
    df_coords_1_3.loc[((df_coords_1_3['que_end_previous'] > df_coords_1_3['que_start'])), 'que_start_final'] = df_coords_1_3['que_end_previous']

    df_coords_1_3['ref_end_previous_final']=df_coords_1_3['ref_end_previous_final'].astype(int)
    df_coords_1_3['ref_start_final']=df_coords_1_3['ref_start_final'].astype(int)
    df_coords_1_3['que_end_previous_final']=df_coords_1_3['que_end_previous_final'].astype(int)
    df_coords_1_3['que_start_final']=df_coords_1_3['que_start_final'].astype(int)
    
    df_coords_1_3_SV=df_coords_1_3[['que','que_end_previous_final','que_start_final']]
    df_coords_1_3_SV=pybedtools.BedTool.from_dataframe(df_coords_1_3_SV)
    df_coords_1_3_SV_final=df_coords_1_3_SV.subtract(df_coords_1_3_align,A=True)
    
    df_coords_1_3_SV_final = df_coords_1_3_SV_final.to_dataframe()
    df_coords_1_3_SV_final.columns=['que','que_end_previous_final','que_start_final']
    df_coords_1_3_SV_final_merge=pd.merge(df_coords_1_3_SV_final,df_coords_1_3,how="inner",on=['que','que_end_previous_final','que_start_final'])
    
    df_coords_1_3_SV_final_merge_output=df_coords_1_3_SV_final_merge[['ref','ref_end_previous_final','ref_start_final','que','que_end_previous_final','que_start_final']]
    df_coords_1_3_SV_final_merge_output.to_csv(output_for_genotype_SV_coords_1_3,sep="\t",index=False)









    ###species_1_4
    df_coords_1_4=pd.read_csv(input_coords_1_4,sep="\t",header=None)
    df_coords_1_4.columns=['que','a','b','c','d','ref','que_start','que_end','ref_start','ref_end','e','f','g','h','i','j','k','strand','m','n','o']
    df_coords_1_4=df_coords_1_4[['que','ref','que_start','que_end','ref_start','ref_end','strand']]
    df_coords_1_4=df_coords_1_4[df_coords_1_4['que_start']<=df_coords_1_4['que_end']]

    df_coords_1_4['que_end_previous_1'] = df_coords_1_4['que_end'].shift(1,fill_value=0)
    df_coords_1_4['que_end_previous_2'] = df_coords_1_4['que_end'].shift(2,fill_value=0)
    df_coords_1_4['que_end_previous_3'] = df_coords_1_4['que_end'].shift(3,fill_value=0)
    df_coords_1_4['que_end_previous_4'] = df_coords_1_4['que_end'].shift(4,fill_value=0)

    df_coords_1_4['que_start_next_1'] = df_coords_1_4['que_start'].shift(-1,fill_value=0)
    df_coords_1_4['que_start_next_2'] = df_coords_1_4['que_start'].shift(-2,fill_value=0)
    df_coords_1_4['que_start_next_3'] = df_coords_1_4['que_start'].shift(-3,fill_value=0)
    df_coords_1_4['que_start_next_4'] = df_coords_1_4['que_start'].shift(-4,fill_value=0)

    df_coords_1_4=df_coords_1_4[
       (df_coords_1_4['que_start']-df_coords_1_4['que_end_previous_1'] > -10  ) &
        (df_coords_1_4['que_start']-df_coords_1_4['que_end_previous_2'] > -10  ) &
        (df_coords_1_4['que_start']-df_coords_1_4['que_end_previous_3'] > -10  ) &
        (df_coords_1_4['que_start']-df_coords_1_4['que_end_previous_4'] > -10  ) &

        (df_coords_1_4['que_start_next_1']-df_coords_1_4['que_end'] > -10  ) &
        (df_coords_1_4['que_start_next_2']-df_coords_1_4['que_end'] > -10  ) &
        (df_coords_1_4['que_start_next_3']-df_coords_1_4['que_end'] > -10  ) &
        (df_coords_1_4['que_start_next_4']-df_coords_1_4['que_end'] > -10  ) &

        (df_coords_1_4['que_end_previous_1'] != 0) &
        (df_coords_1_4['que_end_previous_2'] != 0) &
        (df_coords_1_4['que_end_previous_3'] != 0) &
        (df_coords_1_4['que_end_previous_4'] != 0) &

        (df_coords_1_4['que_start_next_1'] != 0) &
        (df_coords_1_4['que_start_next_2'] != 0) &
        (df_coords_1_4['que_start_next_3'] != 0) &
        (df_coords_1_4['que_start_next_4'] != 0)
        ]
    df_coords_1_4_align=df_coords_1_4[['que','que_start','que_end']]
    df_coords_1_4_align.loc[((df_coords_1_4_align['que_start'] <= df_coords_1_4_align['que_end'])), 'que_start_final'] = df_coords_1_4_align['que_start']
    df_coords_1_4_align.loc[((df_coords_1_4_align['que_start'] <= df_coords_1_4_align['que_end'])), 'que_end_final'] = df_coords_1_4_align['que_end']
    df_coords_1_4_align.loc[((df_coords_1_4_align['que_start'] > df_coords_1_4_align['que_end'])), 'que_start_final'] = df_coords_1_4_align['que_end']
    df_coords_1_4_align.loc[((df_coords_1_4_align['que_start'] > df_coords_1_4_align['que_end'])), 'que_end_final'] = df_coords_1_4_align['que_start']
    df_coords_1_4_align=df_coords_1_4_align[['que','que_start_final','que_end_final']]
    df_coords_1_4_align=pybedtools.BedTool.from_dataframe(df_coords_1_4_align)

    df_coords_1_4=df_coords_1_4[['que','ref','que_start','que_end', 'ref_start','ref_end', 'strand']]
    df_coords_1_4['que_end_previous'] = df_coords_1_4['que_end'].shift(1,fill_value=0)
    df_coords_1_4['ref_end_previous'] = df_coords_1_4['ref_end'].shift(1,fill_value=0)

    df_coords_1_4=df_coords_1_4[
        (df_coords_1_4['ref_end_previous'] != 0)
        ]

    df_coords_1_4['gap_que']=df_coords_1_4['que_start']-df_coords_1_4['que_end_previous']
    df_coords_1_4['gap_ref']=df_coords_1_4['ref_start']-df_coords_1_4['ref_end_previous']
    df_coords_1_4=df_coords_1_4[(abs(df_coords_1_4['gap_que'])<=30000) |(abs(df_coords_1_4['gap_ref'])<=30000) ]
    
    df_coords_1_4=df_coords_1_4[ (abs(df_coords_1_4['gap_que']-df_coords_1_4['gap_ref'])>10000) & (abs(df_coords_1_4['gap_que']-df_coords_1_4['gap_ref'])<=1500000) ]
    df_coords_1_4['gap']=abs(df_coords_1_4['gap_que']-df_coords_1_4['gap_ref'])

    df_coords_1_4.loc[((df_coords_1_4['ref_end_previous'] <= df_coords_1_4['ref_start'])), 'ref_end_previous_final'] = df_coords_1_4['ref_end_previous']
    df_coords_1_4.loc[((df_coords_1_4['ref_end_previous'] <= df_coords_1_4['ref_start'])), 'ref_start_final'] = df_coords_1_4['ref_start']
    df_coords_1_4.loc[((df_coords_1_4['ref_end_previous'] > df_coords_1_4['ref_start'])), 'ref_end_previous_final'] = df_coords_1_4['ref_start']
    df_coords_1_4.loc[((df_coords_1_4['ref_end_previous'] > df_coords_1_4['ref_start'])), 'ref_start_final'] = df_coords_1_4['ref_end_previous']

    df_coords_1_4.loc[((df_coords_1_4['que_end_previous'] <= df_coords_1_4['que_start'])), 'que_end_previous_final'] = df_coords_1_4['que_end_previous']
    df_coords_1_4.loc[((df_coords_1_4['que_end_previous'] <= df_coords_1_4['que_start'])), 'que_start_final'] = df_coords_1_4['que_start']
    df_coords_1_4.loc[((df_coords_1_4['que_end_previous'] > df_coords_1_4['que_start'])), 'que_end_previous_final'] = df_coords_1_4['que_start']
    df_coords_1_4.loc[((df_coords_1_4['que_end_previous'] > df_coords_1_4['que_start'])), 'que_start_final'] = df_coords_1_4['que_end_previous']

    df_coords_1_4['ref_end_previous_final']=df_coords_1_4['ref_end_previous_final'].astype(int)
    df_coords_1_4['ref_start_final']=df_coords_1_4['ref_start_final'].astype(int)
    df_coords_1_4['que_end_previous_final']=df_coords_1_4['que_end_previous_final'].astype(int)
    df_coords_1_4['que_start_final']=df_coords_1_4['que_start_final'].astype(int)
    
    df_coords_1_4_SV=df_coords_1_4[['que','que_end_previous_final','que_start_final']]
    df_coords_1_4_SV=pybedtools.BedTool.from_dataframe(df_coords_1_4_SV)
    df_coords_1_4_SV_final=df_coords_1_4_SV.subtract(df_coords_1_4_align,A=True)
    
    df_coords_1_4_SV_final = df_coords_1_4_SV_final.to_dataframe()
    df_coords_1_4_SV_final.columns=['que','que_end_previous_final','que_start_final']
    df_coords_1_4_SV_final_merge=pd.merge(df_coords_1_4_SV_final,df_coords_1_4,how="inner",on=['que','que_end_previous_final','que_start_final'])
    
    df_coords_1_4_SV_final_merge_output=df_coords_1_4_SV_final_merge[['ref','ref_end_previous_final','ref_start_final','que','que_end_previous_final','que_start_final']]
    df_coords_1_4_SV_final_merge_output.to_csv(output_for_genotype_SV_coords_1_4,sep="\t",index=False)








df_hap2_append = pd.DataFrame()
for filename in os.listdir(folder):
    if fnmatch.fnmatch(filename,"output_for_genotype_SV_*."+hap2_ID+".txt"):
        df_hap2=pd.read_csv(filename,sep="\t")
        df_hap2.columns=['ref','ref_start','ref_stop','que','que_start','que_stop']
        
        df_hap2.loc[((df_hap2['ref_start'] <= df_hap2['ref_stop'])), 'ref_stop_final'] = df_hap2['ref_stop']
        df_hap2.loc[((df_hap2['ref_start'] <= df_hap2['ref_stop'])), 'ref_start_final'] = df_hap2['ref_start']
        df_hap2.loc[((df_hap2['ref_start'] > df_hap2['ref_stop'])), 'ref_stop_final'] = df_hap2['ref_start']
        df_hap2.loc[((df_hap2['ref_start'] > df_hap2['ref_stop'])), 'ref_start_final'] = df_hap2['ref_stop']

        df_hap2.loc[((df_hap2['que_start'] <= df_hap2['que_stop'])), 'que_stop_final'] = df_hap2['que_stop']
        df_hap2.loc[((df_hap2['que_start'] <= df_hap2['que_stop'])), 'que_start_final'] = df_hap2['que_start']
        df_hap2.loc[((df_hap2['que_start'] > df_hap2['que_stop'])), 'que_stop_final'] = df_hap2['que_start']
        df_hap2.loc[((df_hap2['que_start'] > df_hap2['que_stop'])), 'que_start_final'] = df_hap2['que_stop']
        
        df_hap2=df_hap2[['ref','ref_start_final','ref_stop_final','que','que_start_final','que_stop_final']]
        
        df_hap2['ref_gap']=df_hap2['ref_stop_final']-df_hap2['ref_start_final']
        df_hap2['que_gap']=df_hap2['que_stop_final']-df_hap2['que_start_final']
        
        df_hap2.loc[((df_hap2['ref_gap'] <30 ) & (df_hap2['que_gap'] >=30)), 'type'] = 'insertion'
        df_hap2.loc[((df_hap2['ref_gap'] >=30 ) & (df_hap2['que_gap'] <30)), 'type'] = 'deletion'
        df_hap2.loc[((df_hap2['ref_gap'] >=30 ) & (df_hap2['que_gap'] >=30) & (df_hap2['ref_gap'] < df_hap2['que_gap'])), 'type'] = 'expansion'
        df_hap2.loc[((df_hap2['ref_gap'] >=30 ) & (df_hap2['que_gap'] >=30) & (df_hap2['ref_gap'] > df_hap2['que_gap'])), 'type'] = 'contraction'
        df_hap2['size']=abs(df_hap2['ref_gap']-df_hap2['que_gap'])
        
        df_hap2_append=df_hap2_append.append(df_hap2)
df_hap2_all= pd.DataFrame(df_hap2_append)
df_hap2_all = df_hap2_all.drop(["ref_gap", "que_gap"], axis=1)
df_hap2_all.columns=['ref','ref_start_final','ref_stop_final','que_RKQQC','que_start_RKQQC','que_stop_RKQQC','type_RKQQC','size_RKQQC']
        

        

df_hap3_append = pd.DataFrame()
for filename in os.listdir(folder):
    if fnmatch.fnmatch(filename,"output_for_genotype_SV_*."+hap3_ID+".txt"):
        df_hap3=pd.read_csv(filename,sep="\t")
        df_hap3.columns=['ref','ref_start','ref_stop','que','que_start','que_stop']
        
        df_hap3['que']=df_hap3['que'].replace('a1_', 'A1_', regex=True)
        
        df_hap3.loc[((df_hap3['ref_start'] <= df_hap3['ref_stop'])), 'ref_stop_final'] = df_hap3['ref_stop']
        df_hap3.loc[((df_hap3['ref_start'] <= df_hap3['ref_stop'])), 'ref_start_final'] = df_hap3['ref_start']
        df_hap3.loc[((df_hap3['ref_start'] > df_hap3['ref_stop'])), 'ref_stop_final'] = df_hap3['ref_start']
        df_hap3.loc[((df_hap3['ref_start'] > df_hap3['ref_stop'])), 'ref_start_final'] = df_hap3['ref_stop']

        df_hap3.loc[((df_hap3['que_start'] <= df_hap3['que_stop'])), 'que_stop_final'] = df_hap3['que_stop']
        df_hap3.loc[((df_hap3['que_start'] <= df_hap3['que_stop'])), 'que_start_final'] = df_hap3['que_start']
        df_hap3.loc[((df_hap3['que_start'] > df_hap3['que_stop'])), 'que_stop_final'] = df_hap3['que_start']
        df_hap3.loc[((df_hap3['que_start'] > df_hap3['que_stop'])), 'que_start_final'] = df_hap3['que_stop']
        
        df_hap3=df_hap3[['ref','ref_start_final','ref_stop_final','que','que_start_final','que_stop_final']]
        
        df_hap3['ref_gap']=df_hap3['ref_stop_final']-df_hap3['ref_start_final']
        df_hap3['que_gap']=df_hap3['que_stop_final']-df_hap3['que_start_final']
        
        df_hap3.loc[((df_hap3['ref_gap'] <30 ) & (df_hap3['que_gap'] >=30)), 'type'] = 'insertion'
        df_hap3.loc[((df_hap3['ref_gap'] >=30 ) & (df_hap3['que_gap'] <30)), 'type'] = 'deletion'
        df_hap3.loc[((df_hap3['ref_gap'] >=30 ) & (df_hap3['que_gap'] >=30) & (df_hap3['ref_gap'] < df_hap3['que_gap'])), 'type'] = 'expansion'
        df_hap3.loc[((df_hap3['ref_gap'] >=30 ) & (df_hap3['que_gap'] >=30) & (df_hap3['ref_gap'] > df_hap3['que_gap'])), 'type'] = 'contraction'
        df_hap3['size']=abs(df_hap3['ref_gap']-df_hap3['que_gap'])
        
        df_hap3_append=df_hap3_append.append(df_hap3)
df_hap3_all= pd.DataFrame(df_hap3_append)
df_hap3_all = df_hap3_all.drop(["ref_gap", "que_gap"], axis=1)
df_hap3_all.columns=['ref','ref_start_final','ref_stop_final','que_pgt_a1','que_start_pgt_a1','que_stop_pgt_a1','type_pgt_a1','size_pgt_a1']
 
    
    
    
    
    
df_hap4_append = pd.DataFrame()
for filename in os.listdir(folder):
    if fnmatch.fnmatch(filename,"output_for_genotype_SV_*."+hap4_ID+".txt"):
        df_hap4=pd.read_csv(filename,sep="\t")
        df_hap4.columns=['ref','ref_start','ref_stop','que','que_start','que_stop']
        
        df_hap4['que']=df_hap4['que'].replace('b1_', 'B1_', regex=True)
        
        df_hap4.loc[((df_hap4['ref_start'] <= df_hap4['ref_stop'])), 'ref_stop_final'] = df_hap4['ref_stop']
        df_hap4.loc[((df_hap4['ref_start'] <= df_hap4['ref_stop'])), 'ref_start_final'] = df_hap4['ref_start']
        df_hap4.loc[((df_hap4['ref_start'] > df_hap4['ref_stop'])), 'ref_stop_final'] = df_hap4['ref_start']
        df_hap4.loc[((df_hap4['ref_start'] > df_hap4['ref_stop'])), 'ref_start_final'] = df_hap4['ref_stop']

        df_hap4.loc[((df_hap4['que_start'] <= df_hap4['que_stop'])), 'que_stop_final'] = df_hap4['que_stop']
        df_hap4.loc[((df_hap4['que_start'] <= df_hap4['que_stop'])), 'que_start_final'] = df_hap4['que_start']
        df_hap4.loc[((df_hap4['que_start'] > df_hap4['que_stop'])), 'que_stop_final'] = df_hap4['que_start']
        df_hap4.loc[((df_hap4['que_start'] > df_hap4['que_stop'])), 'que_start_final'] = df_hap4['que_stop']
        
        df_hap4=df_hap4[['ref','ref_start_final','ref_stop_final','que','que_start_final','que_stop_final']]
        
        df_hap4['ref_gap']=df_hap4['ref_stop_final']-df_hap4['ref_start_final']
        df_hap4['que_gap']=df_hap4['que_stop_final']-df_hap4['que_start_final']
        
        df_hap4.loc[((df_hap4['ref_gap'] <30 ) & (df_hap4['que_gap'] >=30)), 'type'] = 'insertion'
        df_hap4.loc[((df_hap4['ref_gap'] >=30 ) & (df_hap4['que_gap'] <30)), 'type'] = 'deletion'
        df_hap4.loc[((df_hap4['ref_gap'] >=30 ) & (df_hap4['que_gap'] >=30) & (df_hap4['ref_gap'] < df_hap4['que_gap'])), 'type'] = 'expansion'
        df_hap4.loc[((df_hap4['ref_gap'] >=30 ) & (df_hap4['que_gap'] >=30) & (df_hap4['ref_gap'] > df_hap4['que_gap'])), 'type'] = 'contraction'
        df_hap4['size']=abs(df_hap4['ref_gap']-df_hap4['que_gap'])
        
        df_hap4_append=df_hap4_append.append(df_hap4)
df_hap4_all= pd.DataFrame(df_hap4_append)
df_hap4_all = df_hap4_all.drop(["ref_gap", "que_gap"], axis=1)
df_hap4_all.columns=['ref','ref_start_final','ref_stop_final','que_pgt_b1','que_start_pgt_b1','que_stop_pgt_b1','type_pgt_b1','size_pgt_b1']


df_merge=pd.merge(pd.merge(df_hap2_all,df_hap3_all,how="outer",on=['ref','ref_start_final','ref_stop_final']),df_hap4_all,how="outer",on=['ref','ref_start_final','ref_stop_final'])


df_merge['que_RKQQC']=df_merge["que_RKQQC"].fillna(df_merge["que_pgt_a1"])
df_merge['que_start_RKQQC']=df_merge["que_start_RKQQC"].fillna(df_merge["que_start_pgt_a1"])
df_merge['que_stop_RKQQC']=df_merge["que_stop_RKQQC"].fillna(df_merge["que_stop_pgt_a1"])
df_merge['type_RKQQC']=df_merge['type_RKQQC'].fillna(df_merge["type_pgt_a1"])
df_merge['size_RKQQC']=df_merge['size_RKQQC'].fillna(df_merge["size_pgt_a1"])

df_merge['que_RKQQC']=df_merge["que_RKQQC"].fillna(df_merge["que_pgt_b1"])
df_merge['que_start_RKQQC']=df_merge["que_start_RKQQC"].fillna(df_merge["que_start_pgt_b1"])
df_merge['que_stop_RKQQC']=df_merge["que_stop_RKQQC"].fillna(df_merge["que_stop_pgt_b1"])
df_merge['type_RKQQC']=df_merge['type_RKQQC'].fillna(df_merge["type_pgt_b1"])
df_merge['size_RKQQC']=df_merge['size_RKQQC'].fillna(df_merge["size_pgt_b1"])

df_merge=df_merge[['ref','ref_start_final','ref_stop_final','que_RKQQC','que_start_RKQQC','que_stop_RKQQC','type_RKQQC','size_RKQQC']]
df_merge.columns=['ref','ref_start','ref_stop','que','que_start','que_stop','SV_type','size_calculate']
df_merge['name']=df_merge['ref']+"_"+df_merge['ref_start'].astype(str)+"-"+df_merge['ref_stop'].astype(str)
df_merge=df_merge.rename(columns = {'size_calculate':'size'})

df_merge_new=df_merge.drop_duplicates(keep='first')
df_merge_new=df_merge_new[['ref','ref_start','ref_stop','que','que_start','que_stop','SV_type','size','name']]


df_merge_new.to_csv(outputfile_SV_info,sep="\t",index=False)






files = os.listdir(folder)
for item in files:
    if item.startswith(("output_for_genotype_SV_coords.")):
        os.remove(os.path.join(folder, item))