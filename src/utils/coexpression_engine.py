import pandas as pd 
import numpy as np 
import os 
import time 
from src .project_paths import ANNOTATION_DIR 

class CoexpressionEngine :
    def __init__ (self ,tpm_file ):
        self .tpm_file =tpm_file 
        self .matrix =None 
        self .genes =[]

    def load_matrix (self ):
        """Load massive expression matrix into memory and prepare for normalization"""
        print (f"⏳ CurrentlyloadExpressionquantity matrix (this may take a few seconds)...")
        start_time =time .time ()

        #Read data, set first column as row names (Index)
        self .matrix =pd .read_csv (self .tpm_file ,index_col =0 )

        #Extremely important: filter out dead genes that barely express in any sample
        #e.g., require TPM > 1 in at least 3 samples
        mask =(self .matrix >1 ).sum (axis =1 )>=3 
        self .matrix =self .matrix [mask ]

        #Log transformation to reduce over-influence of highly expressed genes (TPM → log2(TPM+1))
        self .matrix =np .log2 (self .matrix +1 )

        self .genes =list (self .matrix .index )
        print (f"✅ loadCompleted！Time-consuming {time .time ()-start_time :.2f} second。co-retained {len (self .genes )} high-confidence activeGene。")

    def find_top_coexpressed_genes (self ,target_gene ,top_n =5 ):
        """Calculate Pearson correlation between target gene and entire gene set, return top N most similar 'comrades'"""
        if self .matrix is None :
            self .load_matrix ()

        if target_gene not in self .genes :
            print (f"⚠️ {target_gene } inExpressionFailed quality control or not detected in the matrix，Nomethod to calculate co-Expression。")
            return None 

        print (f"📡 CurrentlyfullGenesearch within the group {target_gene } their commonExpressioncomrade...")

        #Retrieve expression vector for target gene
        target_vector =self .matrix .loc [target_gene ]

        #Vectorized computation of Pearson correlation (ultra-fast!)
        #Pandas' corrwith calculates correlation between each row of the matrix and the target vector.
        correlations =self .matrix .corrwith (target_vector ,axis =1 )

        #Discard the 1.0 values from self-correlation and NaN entries.
        correlations =correlations .drop (index =target_gene ).dropna ()

        #Sort by correlation strength in descending order.
        sorted_corr =correlations .sort_values (ascending =False )

        #Extract the Top N results.
        top_genes =sorted_corr .head (top_n )

        results =[]
        for gene ,corr in top_genes .items ():
            results .append ({"Gene_ID":gene ,"Correlation":float (corr )})

        print (f"🎯 Successfullock {top_n } close combat comrades！highest correlation: {results [0 ]['Correlation']:.3f}")
        return pd .DataFrame (results )

if __name__ =="__main__":
    import argparse 
    import os 
    from rich .console import Console 
    from rich .markdown import Markdown 

    #Import the custom Gene System (Phase 1).
    try :
        from src .utils .gene_annotator import RiceGeneAnnotator 
    except ImportError :
    #Support running directly from this directory.
        import sys 
        sys .path .append (os .path .join (os .path .dirname (__file__ ),'../../'))
        from src .utils .gene_annotator import RiceGeneAnnotator 

    console =Console ()

    parser =argparse .ArgumentParser (description ="Legacy co-expression calculator for ad hoc matrix exploration.")
    parser .add_argument ("--matrix",required =True ,help ="Path to an explicit expression matrix file.")
    args =parser .parse_args ()

    #1. Start the Common Expression Engine (Phase 2).
    tpm_path =args .matrix 
    if not os .path .exists (tpm_path ):
        console .print (f"[red]❌ cannot findExpressionamountFile: {tpm_path }[/red]")
        exit ()

    engine =CoexpressionEngine (tpm_path )
    engine .load_matrix ()

    #2. Start the Local Gene System (Phase 1).
    GFF_FILE =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"locus.gff")
    MAP_FILE =str (ANNOTATION_DIR /"RAP-MSU.txt")
    annotator =RiceGeneAnnotator (GFF_FILE ,MAP_FILE )
    annotator .build_database ()

    #Begin Testing.
    test_gene ="Os01g0100100"
    result_df =engine .find_top_coexpressed_genes (test_gene ,top_n =5 )

    if result_df is not None :
        console .print ("\n[bold yellow]🔍 Currently querying the local Gene System to verify detailed background of these 5 genes...[/bold yellow]")

        detailed_results =[]
        for _ ,row in result_df .iterrows ():
            gene_id =row ['Gene_ID']
            corr =row ['Correlation']

            #Investigate thoroughly!
            info =annotator .query (gene_id )
            if info :
                desc =info .get ('Description','No Description')
                syns =", ".join (info .get ('Synonyms',[]))
                detailed_results .append ({
                "Gene_ID":gene_id ,
                "Corr":f"{corr :.3f}",
                "Description":desc [:60 ]+"..."if len (desc )>60 else desc ,#Truncate long outputs for better readability.
                "Synonyms":syns if syns else "-"
                })
            else :
                detailed_results .append ({
                "Gene_ID":gene_id ,
                "Corr":f"{corr :.3f}",
                "Description":"Unknown",
                "Synonyms":"-"
                })

        detailed_df =pd .DataFrame (detailed_results )

        console .print ("\n"+"="*80 )
        console .print (f"[bold cyan]🌟 {test_gene } core co-ExpressionNetwork (Top 5 detailed intelligence report)[/bold cyan]")
        console .print ("="*80 )
        console .print (Markdown (detailed_df .to_markdown (index =False )))
