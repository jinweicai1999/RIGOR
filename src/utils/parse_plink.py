import pandas as pd 
import os 
from src .project_paths import PROCESSED_DIR 

def extract_gene_stats (plink_file ,trait_name ,genes_info ):
    """Extract statistical information for target gene from Plink file"""
    print (f"🔍 Currentlyfrom {plink_file } extract fromGeneassociationSignal...")
    df =pd .read_csv (plink_file ,sep ='\s+')
    parsed_data =[]
    for gene in genes_info :
    #Define gene territory (extend upstream and downstream by 2000 bp to capture regulation regions)
        upstream =gene ['start']-2000 
        downstream =gene ['end']+2000 

        #Identify SNPs falling within this interval
        gene_snps =df [(df ['CHR']==gene ['chr'])&(df ['BP']>=upstream )&(df ['BP']<=downstream )]

        if not gene_snps .empty :
        #Find the most significant SNP (with the smallest P-value)
            lead_snp =gene_snps .loc [gene_snps ['P'].idxmin ()]
            pval =lead_snp ['P']

            #Filter out non-significant SNPs
            if pval <0.05 :
                parsed_data .append ({
                "Gene_ID":gene ['name'],
                "Trait":trait_name ,
                "Association_Type":"GWAS_SNP",
                "P_Value":pval ,
                "Effect_Size_or_Correlation":lead_snp ['BETA'],#Original effect size!
                "Evidence_Level":"High"if pval <1e-5 else "Medium"
                })

    return pd .DataFrame (parsed_data )

if __name__ =="__main__":
#Assume we have located these three genes approximately on the chromosome
    target_genes =[
    {"name":"OsNCED3","chr":3 ,"start":15000000 ,"end":15005000 },
    {"name":"OsLEA3","chr":5 ,"start":20000000 ,"end":20005000 },
    {"name":"SAPK2","chr":7 ,"start":12000000 ,"end":12005000 },
    {"name":"OsNAC-Novel","chr":1 ,"start":4999000 ,"end":5005000 }
    ]

    plink_file =str (PROCESSED_DIR /"mock_gwas_drought.assoc.linear")

    final_df =extract_gene_stats (plink_file ,"Drought Tolerance",target_genes )

    #Save as a standard table recognized by Agent
    final_df .to_csv (PROCESSED_DIR /"rice_gene_trait_stats.csv",index =False )
    print ("\n✅ Data transformation successful! The final table fed to Agent is as follows:")
    print (final_df .to_markdown (index =False ))
