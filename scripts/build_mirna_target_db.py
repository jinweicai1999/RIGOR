import os 
import pandas as pd 
import sys 
from pathlib import Path 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .project_paths import ANNOTATION_DIR ,PROCESSED_DIR 

input_file =str (ANNOTATION_DIR /"Oryza_sativa_targetGene.txt")
output_dir =str (PROCESSED_DIR )
output_file =os .path .join (output_dir ,"clean_mirna_targets.csv")

if not os .path .exists (output_dir ):
    os .makedirs (output_dir )

print (f"🚀 start cleaning raw targetsFile: {input_file }")
data =[]

with open (input_file ,'r',encoding ='utf-8',errors ='ignore')as f :
    for line in f :
    #Skip header or empty lines
        if line .lower ().startswith ("mirna")or not line .strip ():
            continue 

            #Split by force, ignoring spaces and tab characters
        parts =line .strip ().split ()
        if len (parts )>=6 :
        #1. Lowercase miRNA ID
            mirna_clean =parts [0 ].lower ().strip ()

            #2. Core cleaning: remove LOC_ and .1 etc. transcript suffixes
            target_clean =parts [1 ].split ('.')[0 ]

            #3. Extract evidence for subsequent sorting
            deg_raw =parts [4 ]
            psr_raw =parts [2 ]
            lib_raw =parts [5 ]

            #Calculate library count
            lib_count =len (lib_raw .split (','))if lib_raw .upper ()!='NA'else 0 

            #Convert to numeric type for sorting
            try :deg_score =float (deg_raw )
            except :deg_score =999.0 

            try :psr_score =float (psr_raw )
            except :psr_score =999.0 

            data .append ({
            'miRNA':mirna_clean ,
            'Target_ID':target_clean ,#Clean MSU/RAP compatible ID
            'psRNAtarget':psr_raw ,
            'psRNAtarget_score':psr_score ,
            'Degradome':deg_raw ,
            'Degradome_score':deg_score ,
            'Libraries':lib_raw ,
            'Lib_count':lib_count 
            })

df =pd .DataFrame (data )
df .to_csv (output_file ,index =False )
print (f"✅ CleaningCompleted！totalTreatment {len (df )} target records。")
print (f"📁 cleanDatabaseHassave to: {output_file }")
