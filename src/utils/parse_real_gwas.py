import pandas as pd 
import os 
from src .project_paths import PROCESSED_DIR 

def scan_gene_across_traits (gene_info ,traits ,gwas_dir =str (PROCESSED_DIR )):
    """Radar Scan: Extract signal for target gene from GWAS results of the phenotype"""
    print (f"\n📡 CurrentlyStart GWAS radar，panoramicScanGene: {gene_info ['name']} (Chr{gene_info ['chr']}: {gene_info ['start']}-{gene_info ['end']})")

    #Expand 2kb flanking sequence to capture variations in the start region
    upstream =gene_info ['start']-2000 
    downstream =gene_info ['end']+2000 

    all_signals =[]

    for trait in traits :
        file_path =os .path .join (gwas_dir ,f"gwas_{trait }.assoc.linear")
        if not os .path .exists (file_path ):
            continue 

        print (f"  -> CurrentlyScantrait: {trait } ...",end ="")

        #The 1.6GB file is too large; we use an iterator to read it in chunks, saving memory significantly
        chunk_iter =pd .read_csv (file_path ,sep ='\s+',chunksize =100000 ,usecols =['CHR','SNP','BP','BETA','P'])

        best_snp =None 
        min_p =1.0 

        for chunk in chunk_iter :
        #Filter SNPs that fall within the target gene interval
            gene_snps =chunk [(chunk ['CHR']==gene_info ['chr'])&(chunk ['BP']>=upstream )&(chunk ['BP']<=downstream )]

            if not gene_snps .empty :
            #Find the SNP with the smallest P-value in this region
                chunk_best =gene_snps .loc [gene_snps ['P'].idxmin ()]
                if chunk_best ['P']<min_p :
                    min_p =chunk_best ['P']
                    best_snp =chunk_best 

                    #If a SNP is found and is reasonably significant (here relaxed to P<0.05 to see trends; strict threshold can be set to 1e-5)
        if best_snp is not None and min_p <0.05 :
            print (f" 🎯 Hit! strongestSignal P={min_p :.2e}, effect size BETA={best_snp ['BETA']:.3f}")
            all_signals .append ({
            "Gene_ID":gene_info ['name'],
            "Trait":trait ,
            "Association_Type":"GWAS_Lead_SNP",
            "P_Value":min_p ,
            "Effect_Size_or_Correlation":best_snp ['BETA'],
            "Evidence_Level":"High"if min_p <1e-5 else ("Medium"if min_p <1e-3 else "Low")
            })
        else :
            print ("⚪ No significant signal detected")

    return pd .DataFrame (all_signals )

if __name__ =="__main__":
#The 10 traits we are investigating
    traits =[
    "Heading_date","Plant_height","Num_panicles","Num_effective_panicles",
    "Yield","Grain_weight","Spikelet_length","Grain_length","Grain_width","Grain_thickness"
    ]

    #The gene we are testing (please ensure SAPK2's chromosome and position are correct; here assumed to be Chr7: 12000000)
    #In practice, you should consult RAP-DB to replace these with actual start and end coordinates
    test_gene ={"name":"SAPK2","chr":7 ,"start":25715835 ,"end":25721775 }

    # RunningScan
    final_df =scan_gene_across_traits (test_gene ,traits )

    if not final_df .empty :
    #Save it, ready to feed into the large model at any time
        output_path =PROCESSED_DIR /"rice_gene_trait_stats.csv"
        final_df .to_csv (output_path ,index =False )
        print (f"\n✅ radarScandone！Hasextract real hardData，generate feeding table：{output_path }")
        try :
            from tabulate import tabulate 
            print (final_df .to_markdown (index =False ))
        except :
            print (final_df )
    else :
        print ("\n⚠️ No significant signal (P<0.05) was scanned in this gene interval. You may want to check if the coordinates are accurate or test other genes.")
