import os 
import pandas as pd 
from src .project_paths import PROCESSED_DIR 

class NcRnaAgent :
    def __init__ (self ,clean_file =str (PROCESSED_DIR /"clean_mirna_targets.csv")):
        self .clean_file =clean_file 
        self .df =None 
        self .lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"
        self ._load_data ()

    def _txt (self ,zh ,en ):
        return en if self .lang =="en"else zh 

    def _load_data (self ):
        print (self ._txt ("Currently loading standardized miRNA target database...","Loading standardized miRNA target database..."))
        if not os .path .exists (self .clean_file ):
            print (self ._txt (f"cannot find cleanedDatabase: {self .clean_file }。Please firstRunning scripts/build_mirna_target_db.py",
            f"Clean target database not found: {self .clean_file }. Run scripts/build_mirna_target_db.py first."))
            return 
        self .df =pd .read_csv (self .clean_file )

    def find_top_targets (self ,mirna_id ,annotator ,is_star_gene_func =None ,top_n =3 ):
        """Core logic: prioritize returning target genes with HasAnnotation; candidate pool expands up to top_n*5 records if needed"""
        if self .df is None or self .df .empty :
            return []

            #1. Find targets
        target_df =self .df [self .df ['miRNA']==mirna_id .lower ().strip ()]

        if target_df .empty :
            print (self ._txt (f"Databasenot included {mirna_id }。",f"{mirna_id } is not found in database."))
            return []

        print (self ._txt (f"Successfullock {mirna_id } 's {len (target_df )} candidate targets。",f"Located {len (target_df )} candidate targets for {mirna_id }."))

        #2. Sort: directly use precomputed Has values (degradome-prioritized)
        target_df =target_df .sort_values (
        by =['Degradome_score','psRNAtarget_score','Lib_count'],
        ascending =[True ,True ,False ]
        )

        #3. Expand candidate pool ceiling, prioritizing enough HasAnnotation targets to reach top_n
        scan_limit =min (len (target_df ),top_n *5 )
        annotated_targets =[]
        unannotated_targets =[]

        for _ ,row in target_df .head (scan_limit ).iterrows ():
            base_id =str (row ['Target_ID']).strip ()
            target_info =annotator .query (base_id )
            if not target_info and "."not in base_id :
                target_info =annotator .query (f"{base_id }.1")

            entry ={
            "Degradome":row ['Degradome'],
            "Degradome_score":row ['Degradome_score'],
            "psRNAtarget":row ['psRNAtarget'],
            "psRNAtarget_score":row ['psRNAtarget_score'],
            "Libraries":row ['Libraries'],
            "Lib_count":row ['Lib_count'],
            }

            if target_info :
                entry .update ({
                "RAP_ID":target_info ['RAP_ID'],
                "Primary_Symbol":target_info .get ('Primary_Symbol',''),
                "Synonyms":target_info .get ('Synonyms',[]),
                "PMIDs":target_info .get ('PMIDs',[]),
                })
                annotated_targets .append (entry )
            else :
                entry .update ({
                "RAP_ID":base_id ,
                "Primary_Symbol":base_id ,
                "Synonyms":[],
                "PMIDs":[],
                })
                unannotated_targets .append (entry )

            if len (annotated_targets )>=top_n :
                break 

                #4. If insufficient HasAnnotation, fill with NoAnnotation entries
        valid_targets =annotated_targets [:top_n ]
        if len (valid_targets )<top_n :
            shortage =top_n -len (valid_targets )
            valid_targets +=unannotated_targets [:shortage ]
            if unannotated_targets :
                print (self ._txt (
                f"note：HasAnnotationtargetGeneInadequate {top_n } each，Hasuse {min (shortage ,len (unannotated_targets ))} eachNoAnnotationtargetGenefill in。",
                f"Note: fewer than {top_n } annotated targets found; {min (shortage ,len (unannotated_targets ))} unannotated target(s) added as fallback."
                ))

        return valid_targets 
