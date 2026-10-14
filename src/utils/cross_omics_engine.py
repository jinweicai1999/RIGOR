import pandas as pd 
import numpy as np 
import os 
from rich .console import Console 
from src .utils .rna_type_detector import RnaTypeDetector 
from src .project_paths import ANNOTATION_DIR 

class CrossOmicsEngine :
    def __init__ (self ,
    protein_matrix =str (ANNOTATION_DIR /"protein_coding_gene_matrix.coexp.txt"),
    mirna_matrix =str (ANNOTATION_DIR /"miRNA_matrix.coexp.txt"),
    lncrna_matrix =str (ANNOTATION_DIR /"LNCandPHAS_matrix.coexp.txt"),
    phasirna_matrix =str (ANNOTATION_DIR /"phasiRNA_matrix.coexp.txt")):
        self .console =Console ()
        self .files ={
        'protein':protein_matrix ,
        'miRNA':mirna_matrix ,
        'lncRNA':lncrna_matrix ,
        'phasiRNA':phasirna_matrix 
        }
        self .protein_df =None 
        self .nc_dfs ={}#Cache non-coding gene matrices.
        self .detector =RnaTypeDetector ()#Smart Type Identifier.

    def _load_protein_matrix (self ):
        if self .protein_df is None :
            self .console .print (" [cyan]⏳ Loading protein-coding gene matrix (37MB, initial load may take a few seconds)...[/cyan]")
            self .protein_df =pd .read_csv (self .files ['protein'],sep ='\t',index_col =0 )

            #Clean MSU IDs: Convert LOC_Os01g01010 to standard format Os01g01010.
            #Handle unconventional gene names found in files.
            self .protein_df .index =self .protein_df .index .astype (str ).str .replace ('LOC_','').str .split ('.').str [0 ]

    def _load_nc_matrix (self ,nc_type ):
        if nc_type not in self .nc_dfs :
            self .console .print (f"  [cyan]⏳ Currentlyload {nc_type } Matrix...[/cyan]")
            df =pd .read_csv (self .files [nc_type ],sep ='\t',index_col =0 )
            #miRNA names are case-insensitive for matching purposes.
            df .index =df .index .astype (str ).str .lower ()
            self .nc_dfs [nc_type ]=df 
        return self .nc_dfs [nc_type ]

    def _determine_nc_type (self ,nc_id ):
        """Smart Type Identifier (supports multiple ID formats and sequence validation)."""
        #Attempt to retrieve sequence from database.
        sequence =self .detector .get_sequence_from_id (nc_id )

        #Integrated identification.
        rna_type ,method =self .detector .detect (rna_id =nc_id ,sequence =sequence )

        if rna_type and rna_type !='miRNA_precursor':
            self .console .print (f"  [green]✓ identified as {rna_type } (Method: {method })[/green]")
            return rna_type 
        elif rna_type =='miRNA_precursor':
            self .console .print (f"  [yellow]⚠ detected miRNA ProgenitorSequence，will serve as miRNA Treatment[/yellow]")
            return 'miRNA'

        return None 

    def find_cross_omics_partners (self ,nc_id ,top_n =200 ):
        """
        Core logic: Compute genome-wide correlation.
        For miRNA/phasiRNA, return top N with strongest negative correlation.
        For lncRNA, return top N with strongest positive correlation.
        """
        nc_type =self ._determine_nc_type (nc_id )
        if not nc_type :
            self .console .print (f"  [red]❌ Nomethod to identify nonCoding RNA type: {nc_id }[/red]")
            return None 

        nc_df =self ._load_nc_matrix (nc_type )

        #Search for target gene lines.
        query_id =nc_id .lower ()
        #Since miRNA matrix indices include suffixes (e.g., _PmiREN), fuzzy matching is required.
        matched_index =[idx for idx in nc_df .index if query_id in idx ]

        if not matched_index :
            self .console .print (f"  [yellow]⚠️ Databasethis non- was not included inCoding RNA 'sExpressionData: {nc_id }[/yellow]")
            return None 

        exact_match =matched_index [0 ]
        nc_expr =nc_df .loc [exact_match ]

        self .console .print (f"  [green]✅ SuccessfulextractExpressionpolyline (Samplenumber: {len (nc_expr )})[/green]")

        #Ensure protein matrix is loaded.
        self ._load_protein_matrix ()

        #Core alignment: Extract shared sample columns.
        common_samples =nc_expr .index .intersection (self .protein_df .columns )
        if len (common_samples )<10 :
            self .console .print (" [red]❌ Insufficient overlapping samples; cannot compute correlation.[/red]")
            return None 

        aligned_nc =nc_expr [common_samples ]
        aligned_protein =self .protein_df [common_samples ]

        self .console .print (f"  [cyan]📡 Currentlywith fullGenehistoneCodingGenecompute Pearson correlation...[/cyan]")

        #Rapidly calculate full gene set correlation relationships.
        correlations =aligned_protein .corrwith (aligned_nc ,axis =1 )

        #Core branching: Apply biological rules to select positive or negative extremes.
        if nc_type in ['miRNA','phasiRNA']:
            self .console .print (" [cyan]🔍 Mechanism inference: Identifying strongly degraded/inhibited targets (selecting strong negative correlation)[/cyan]")
            #Filter out positive correlations and keep only negative ones, then sort by absolute value in descending order
            valid_corrs =correlations [correlations <-0.3 ].sort_values (ascending =True )

            #Since these are negative numbers, sorting from smallest to largest places the strongest negative correlations first. We'll take the absolute value of weights for easier uniform treatment later
            valid_corrs =valid_corrs .head (top_n ).abs ()

        else :# lncRNA
            self .console .print ("🔍 [cyan]Mechanism determination: Searching for co-transcribed/positively regulated partners (selecting strong positive correlations)[/cyan]")
            valid_corrs =correlations [correlations >0.3 ].sort_values (ascending =False ).head (top_n )

        if valid_corrs .empty :
            self .console .print ("[yellow]⚠️ No proteins with significant correlation (|r| > 0.3) were found.[/yellow]")
            return None 

            #Convert to DataFrame and standardize format for seamless integration with MacroNetworkAgent
        results =pd .DataFrame ({
        'Gene_ID':valid_corrs .index ,#Note: These are MSU format IDs (Os...)
        'Correlation':valid_corrs .values 
        })

        return results 
