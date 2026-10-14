import pandas as pd 
import numpy as np 
from scipy .stats import pearsonr 
from src .project_paths import NCRNA_DIR 

class LncRnaAgent :
    def __init__ (self ):
        self .lnc_matrix_file =str (NCRNA_DIR /"LNCandPHAS_matrix.coexp.txt")
        self .protein_matrix_file =str (NCRNA_DIR /"protein_coding_gene_matrix.coexp.txt")
        self .lnc_matrix =None 
        self .protein_matrix =None 
        self ._lnc_index_map ={}

    def load_matrices (self ):
        """Load co-expression matrix of lncRNA and protein-coding genes"""
        print ("Currently loading lncRNA co-expression matrix...")
        self .lnc_matrix =pd .read_csv (self .lnc_matrix_file ,sep ='\t',index_col =0 )
        self ._lnc_index_map ={str (idx ).strip ().lower ():idx for idx in self .lnc_matrix .index }
        print (f"✅ Hasload {len (self .lnc_matrix )} each lncRNA")

        print ("Currently loading protein-coding gene co-expression matrix...")
        self .protein_matrix =pd .read_csv (self .protein_matrix_file ,sep ='\t',index_col =0 )
        print (f"✅ Hasload {len (self .protein_matrix )} proteinGene")

    def get_top_coexpressed_genes (self ,lncrna_id ,top_n =5 ):
        """Find protein-coding genes highly co-expressed with the lncRNA"""
        if self .lnc_matrix is None or self .protein_matrix is None :
            self .load_matrices ()

            #Check if lncRNA exists
        query_key =str (lncrna_id ).strip ().lower ()
        real_id =self ._lnc_index_map .get (query_key )
        if real_id is None :
        #Fallback: Support minor variants like OSA-LNCNxxxx and osa-lncnxxxx
            query_key_alt =query_key .replace ("osa-","osa_").replace ("_","-")
            real_id =self ._lnc_index_map .get (query_key_alt )
        if real_id is None :
            return None ,f"Not Found lncRNA: {lncrna_id }"

            #Retrieve expression profile of the lncRNA
        lnc_expr =self .lnc_matrix .loc [real_id ]

        #Find common samples
        common_samples =list (set (self .lnc_matrix .columns )&set (self .protein_matrix .columns ))
        if len (common_samples )<10 :
            return None ,f"commonSampletoo few numbers: {len (common_samples )}"

        lnc_expr_common =lnc_expr [common_samples ]

        #If the number of non-zero samples is too low, Pearson correlation cannot be reliably calculated
        nonzero_count =int ((lnc_expr_common !=0 ).sum ())
        if nonzero_count <10 :
            return None ,f"this lncRNA in commonSampleonly inHas {nonzero_count } one non-zeroExpressionvalue，Samplenumber insufficient（need ≥10），Nomethod to calculate credible coExpressioncorrelation。"
        correlations =[]
        for gene_id in self .protein_matrix .index :
            protein_expr =self .protein_matrix .loc [gene_id ,common_samples ]

            #Filter out genes with all zero expression values
            if protein_expr .sum ()==0 or lnc_expr_common .sum ()==0 :
                continue 

            corr ,pval =pearsonr (lnc_expr_common ,protein_expr )
            if not np .isnan (corr ):
                correlations .append ({
                'Gene_ID':gene_id ,
                'Correlation':corr ,
                'P_value':pval 
                })

                #Convert to DataFrame and sort
        corr_df =pd .DataFrame (correlations )
        corr_df =corr_df .sort_values ('Correlation',ascending =False )

        #Return top N positively correlated genes
        top_genes =corr_df .head (top_n )

        return top_genes ,None 
