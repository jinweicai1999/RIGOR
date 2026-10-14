import json 
import os 
import pandas as pd 
from dotenv import load_dotenv 
from src .project_paths import EXPRESSION_DIR ,PROCESSED_DIR 
from src .utils .llm_client import chat_completion_text ,make_openai_client 
from src .utils .gwas_index import (
available_indexed_traits ,
default_index_dir ,
query_trait_chrom_significant ,
)

load_dotenv ()

class GwasAgent :
    def __init__ (self ,gwas_dir =None ,rna_dir =str (EXPRESSION_DIR )):
        self .gwas_dir =gwas_dir or os .getenv ("GWAS_DIR",str (PROCESSED_DIR ))
        self .gwas_index_dir =os .getenv ("GWAS_INDEX_DIR",str (default_index_dir (self .gwas_dir )))
        self .gwas_p_threshold =float (os .getenv ("GWAS_P_THRESHOLD","0.05"))
        self .gwas_trait_background_path =os .getenv ("GWAS_TRAIT_BACKGROUND",os .path .join (self .gwas_dir ,"gwas_trait_background.json"))
        self .gwas_max_top_fraction =float (os .getenv ("GWAS_MAX_TOP_FRACTION","0.05"))
        self .lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"
        self .client =make_openai_client (
        api_key =os .getenv ("LLM_API_KEY"),
        base_url =os .getenv ("LLM_BASE_URL","https://api.longcat.chat/openai").rstrip ("/"),
        )
        self .model_name =os .getenv ("LLM_MODEL_NAME","LongCat-Flash-Lite")#Ensure you fill in your model name

        self .traits =[
        "Heading_date","Plant_height","Num_panicles","Num_effective_panicles",
        "Yield","Grain_weight","Spikelet_length","Grain_length","Grain_width","Grain_thickness"
        ]
        self .trait_background =self ._load_trait_background ()

        #--- 🌟 New: Mount macro community dictionary for 8456 samples ---
        self .assign_file =os .path .join (rna_dir ,"module_assignment.csv")
        self .trait_file =os .path .join (rna_dir ,"module_trait_correlation.csv")
        self .gene_to_module ={}
        self .module_to_traits ={}
        self ._load_macro_modules ()

    def _txt (self ,zh ,en ):
        return en if self .lang =="en"else zh 

    def _load_trait_background (self ):
        if not os .path .exists (self .gwas_trait_background_path ):
            return {}
        try :
            with open (self .gwas_trait_background_path ,"r",encoding ="utf-8")as handle :
                payload =json .load (handle )
            return payload .get ("traits",{})
        except Exception :
            return {}

    def _trait_cutoff (self ,trait ,top_fraction =None ):
        stats =self .trait_background .get (trait ,{})
        cutoffs =stats .get ("p_cutoffs",{})
        top_fraction =self .gwas_max_top_fraction if top_fraction is None else top_fraction 
        key =f"top_{top_fraction:g}"
        if key in cutoffs :
            return float (cutoffs [key ])
        legacy_key ={0.001:"top_0.1_percent",0.01:"top_1_percent",0.05:"top_5_percent"}.get (top_fraction )
        if legacy_key and legacy_key in cutoffs :
            return float (cutoffs [legacy_key ])
        return self .gwas_p_threshold 

    def _trait_normalized_evidence (self ,trait ,p_value ):
        p_value =float (p_value )
        stats =self .trait_background .get (trait ,{})
        cutoffs =stats .get ("p_cutoffs",{})
        strong =float (cutoffs .get ("top_0.001",cutoffs .get ("top_0.1_percent",0.0))or 0.0 )
        moderate =float (cutoffs .get ("top_0.01",cutoffs .get ("top_1_percent",0.0))or 0.0 )
        weak =float (cutoffs .get ("top_0.05",cutoffs .get ("top_5_percent",self .gwas_p_threshold))or self .gwas_p_threshold )
        total_snps =stats .get ("total_snps")
        if strong and p_value <=strong :
            level ="Strong"
            top_fraction =0.001 
        elif moderate and p_value <=moderate :
            level ="Moderate"
            top_fraction =0.01 
        elif p_value <=weak :
            level ="Weak"
            top_fraction =0.05 
        else :
            level ="Filtered"
            top_fraction =None 
        return {
        "Trait_Normalized_Evidence":level ,
        "Trait_Top_Fraction_Max":top_fraction ,
        "Trait_Top_Percent_Max":None if top_fraction is None else top_fraction *100 ,
        "Trait_Background_Total_SNPs":total_snps ,
        "Trait_P_Cutoff_Top_0_1pct":strong or None ,
        "Trait_P_Cutoff_Top_1pct":moderate or None ,
        "Trait_P_Cutoff_Top_5pct":weak or None ,
        }

    def get_expected_gwas_files (self ):
        return [os .path .join (self .gwas_dir ,f"gwas_{trait }.assoc.linear")for trait in self .traits ]

    def inspect_gwas_status (self ):
        expected_files =self .get_expected_gwas_files ()
        existing_files =[path for path in expected_files if os .path .exists (path )]
        indexed_traits =available_indexed_traits (self .gwas_index_dir ,self .traits )
        return {
        "gwas_dir":self .gwas_dir ,
        "gwas_index_dir":self .gwas_index_dir ,
        "expected_count":len (expected_files ),
        "existing_count":len (existing_files ),
        "indexed_trait_count":len (indexed_traits ),
        "indexed_traits":indexed_traits ,
        "existing_files":existing_files ,
        "missing_files":[path for path in expected_files if path not in existing_files ],
        }

    def _load_macro_modules (self ):
        """Load macro qualitative prediction data"""
        print (self ._txt ("Currently loading macro co-expression community data for 8456 samples...","Loading macro co-expression module data (8456 samples)..."))
        if os .path .exists (self .assign_file ):
            df_a =pd .read_csv (self .assign_file )
            for _ ,row in df_a .iterrows ():
                gene_id =str (row ['Gene']).replace ("LOC_","").split ('.')[0 ]

                #🌟 Upgraded to list append mode
                if gene_id not in self .gene_to_module :
                    self .gene_to_module [gene_id ]=[]
                if row ['Module']not in self .gene_to_module [gene_id ]:
                    self .gene_to_module [gene_id ].append (row ['Module'])

        if os .path .exists (self .trait_file ):
            df_t =pd .read_csv (self .trait_file )
            for _ ,row in df_t .iterrows ():
                mod =row ['Module']
                trait =row ['Trait']
                if mod not in self .module_to_traits :
                    self .module_to_traits [mod ]=set ()
                self .module_to_traits [mod ].add (trait )

    def get_macro_evidence (self ,target_gene ,network_genes ,annotator ):
        """Extract macro qualitative evidence for this network"""
        evidence_lines =[]
        all_genes =list (set ([target_gene ]+network_genes ))
        modules_found ={}

        for g in all_genes :
            info =annotator .query (g )
            if not info :continue 

            #Use both MSU ID and RAP ID for dual search
            search_ids =[g ]+info .get ('MSU_IDs',[])
            for sid in search_ids :
            #Found this part:
                clean_id =sid .replace ("LOC_","").split ('.')[0 ]
                if clean_id in self .gene_to_module :
                    mods =self .gene_to_module [clean_id ]#Now this is a list
                    for mod in mods :#Iterate through each module it belongs to
                        if mod not in modules_found :
                            modules_found [mod ]=[]
                        modules_found [mod ].append (g )
                    break #Break upon finding

        for mod ,genes in modules_found .items ():
            traits =self .module_to_traits .get (mod ,set ())
            genes_str =", ".join (genes )
            if traits :
                traits_str =", ".join (traits )
                evidence_lines .append (self ._txt (
                f"- Gene [{genes_str }] belonging to macroNetworkmodule **{mod }**，this module in 8456 copySamplesignificantly enriched in traits：【{traits_str }】",
                f"- Genes [{genes_str }] belong to macro module **{mod }**, significantly enriched for traits in 8456 samples: [{traits_str }]"
                ))
            else :
                evidence_lines .append (self ._txt (
                f"- Gene [{genes_str }] belonging to macroNetworkmodule **{mod }**，But no significant enrichment of agronomic traits was found。",
                f"- Genes [{genes_str }] belong to macro module **{mod }**, with no significant agronomic-trait enrichment found."
                ))

        return "\n".join (evidence_lines )if evidence_lines else self ._txt ("Unable to find module affiliation or trait association for this community in the 8456-sample large network.","No module/trait association found in the 8456-sample macro network.")

    def scan_network_gwas (self ,gene_list ,annotator ):
        """Ultra-fast micro radar: Scan quantitative effect values at specific sites"""
        gene_coords ={}
        for gene in gene_list :
            info =annotator .query (gene )
            if info :
                loc =info ['Location']
                chr_str =loc .split (':')[0 ].replace ('chr','')
                try :
                    chrom =int (chr_str )
                    start =int (loc .split (':')[1 ].split ('-')[0 ])
                    end =int (loc .split ('-')[1 ].split (' ')[0 ])
                    gene_coords [gene ]={'chr':chrom ,'start':start -2000 ,'end':end +2000 }
                except :continue 

        if not gene_coords :return pd .DataFrame ()

        indexed_traits =available_indexed_traits (self .gwas_index_dir ,self .traits )
        if indexed_traits :
            indexed_result =self ._scan_network_gwas_indexed (gene_coords ,indexed_traits )
            if len (indexed_traits )==len (self .traits )or not indexed_result .empty :
                return indexed_result 

        all_signals =[]
        for trait in self .traits :
            file_path =os .path .join (self .gwas_dir ,f"gwas_{trait }.assoc.linear")
            if not os .path .exists (file_path ):continue 

            chunk_iter =pd .read_csv (file_path ,sep =r'\s+',chunksize =200000 ,usecols =['CHR','SNP','BP','BETA','P'])
            best_snps ={gene :{'P':1.0 ,'BETA':0 }for gene in gene_coords .keys ()}
            trait_p_cutoff =self ._trait_cutoff (trait )

            for chunk in chunk_iter :
                for gene ,coords in gene_coords .items ():
                    gene_snps =chunk [(chunk ['CHR']==coords ['chr'])&
                    (chunk ['BP']>=coords ['start'])&
                    (chunk ['BP']<=coords ['end'])]
                    if not gene_snps .empty :
                        min_idx =gene_snps ['P'].idxmin ()
                        best_p =gene_snps .loc [min_idx ,'P']
                        if best_p <best_snps [gene ]['P']:
                            best_snps [gene ]['P']=best_p 
                            best_snps [gene ]['BETA']=gene_snps .loc [min_idx ,'BETA']

            for gene ,data in best_snps .items ():
                if data ['P']<=trait_p_cutoff :
                    evidence =self ._trait_normalized_evidence (trait ,data ['P'])
                    if evidence .get ("Trait_Normalized_Evidence")=="Filtered":
                        continue 
                    record ={
                    "Gene":gene ,
                    "Trait":trait ,
                    "P_Value":data ['P'],
                    "Effect_Size(BETA)":data ['BETA']
                    }
                    record .update (evidence )
                    all_signals .append (record )

        return pd .DataFrame (all_signals )

    def _scan_network_gwas_indexed (self ,gene_coords ,traits ):
        """Fast GWAS scan using prebuilt trait/chromosome Parquet shards."""
        all_signals =[]
        genes_by_chr ={}
        for gene ,coords in gene_coords .items ():
            genes_by_chr .setdefault (coords ["chr"],[]).append ((gene ,coords ))
        for trait in traits :
            trait_p_cutoff =self ._trait_cutoff (trait )
            for chrom ,items in genes_by_chr .items ():
                try :
                    chrom_hits =query_trait_chrom_significant (
                    self .gwas_index_dir ,
                    trait ,
                    p_threshold =trait_p_cutoff ,
                    chrom =chrom ,
                    )
                except Exception :
                    continue 
                if chrom_hits .empty :
                    continue 
                for gene ,coords in items :
                    hits =chrom_hits [
                    (chrom_hits ["BP"]>=coords ["start"])&
                    (chrom_hits ["BP"]<=coords ["end"])
                    ]
                    if hits .empty :
                        continue 
                    min_idx =hits ["P"].idxmin ()
                    best =hits .loc [min_idx ]
                    evidence =self ._trait_normalized_evidence (trait ,best ["P"])
                    if evidence .get ("Trait_Normalized_Evidence")=="Filtered":
                        continue 
                    record ={
                    "Gene":gene ,
                    "Trait":trait ,
                    "P_Value":float (best ["P"]),
                    "Effect_Size(BETA)":float (best ["BETA"]),
                    "SNP":str (best ["SNP"]),
                    "CHR":int (best ["CHR"]),
                    "BP":int (best ["BP"]),
                    }
                    record .update (evidence )
                    all_signals .append (record )
        return pd .DataFrame (all_signals )

    def summarize (self ,target_gene ,target_ids ,gwas_df ,macro_evidence ):
        """Field effect and breeding application deep evaluation"""
        if gwas_df .empty :
            gwas_desc =self ._txt (
            "Local GWAS did not detect trait-normalized agronomic association signals for this network (target gene and allies) in natural populations.",
            "No trait-normalized agronomic GWAS signal was detected for this network (target + partners) in the local population-level scan."
            )
            data_markdown =self ._txt ("(No trait-normalized GWAS signal)","(No trait-normalized GWAS signal)")
        else :
            gwas_desc =self ._txt (
            "Local GWAS captured the following trait-normalized association signals in rice natural populations.",
            "The local GWAS scan captured the following trait-normalized associations in rice natural populations."
            )
            #🌟 Fixed: Use native to_csv(sep='\t') instead of to_markdown, completely eliminating dependency on tabulate!
            data_markdown =gwas_df .to_csv (sep ='\t',index =False )
        if self .lang =="en":
            prompt =f"""
            You are a rigorous rice quantitative genetics and breeding expert.
            Evaluate field-trait effects and breeding trade-offs for target gene '{target_gene }' and its partner network.

            [Background: macro-network phenotype context]
            {macro_evidence }

            [Core measured data: population GWAS quantitative effects]
            {gwas_desc }
            {data_markdown }
            Effect_Size(BETA) is retained as the raw regression coefficient from the GWAS model.
            Do not interpret the sign of BETA as positive or negative regulation, promotion, suppression,
            increase, or decrease. Use P_Value and -log10(P) to discuss association strength only.
            Trait_Normalized_Evidence is assigned within each trait's own genome-wide P-value distribution:
            Strong = top 0.1%, Moderate = top 1%, Weak = top 5%. Treat Strong/Moderate as the main GWAS evidence
            and Weak as suggestive background-level evidence.

            Requirements:
            Output concise, evidence-based Markdown with this structure:

            # {target_gene } GWAS Trait Association Assessment
            ## 1. Significant Trait Association Overview
            ## 2. Candidate Breeding Marker Value
            ## 3. Interpretation Limits and Validation Needs
            Use specific P_Value values and Trait_Normalized_Evidence levels. Do not infer trait direction.
            """
        else :
            prompt =f"""
            You are a rigorous and logically sound expert in rice quantitative genetics and breeding。
            youCurrentlytargetedGene '{target_gene }' its comradesNetworkConduct GWAS trait association analysis。
            
            【backgroundInformation：totalExpressionNetworkphenotypic characteristics of communities】
            {macro_evidence }
            
            【core measuredData：natural population GWAS quantifying effect】
            {gwas_desc }
            {data_markdown }
            *(note：Effect_Size(BETA) is the raw GWAS regression coefficient. Do not interpret its sign as positive/negative regulation, promotion/inhibition, increase/decrease, or yield gain/loss. Use P_Value and -log10(P) to discuss association strength only. If BETA is mentioned, use only its absolute magnitude as an effect-size reference.)*
            *(Trait_Normalized_Evidence is computed within each trait's own genome-wide P-value distribution: Strong = top 0.1%, Moderate = top 1%, Weak = top 5%. Only Strong/Moderate should be treated as core GWAS evidence; Weak is only suggestive.)*
            
            【task requirements】：
            please discard filler，strictly based on the aboveData，Outputfollowing format is exquisite Markdown report：

            # {target_gene } GWAS trait association assessment
            ## 1. significant trait association overview
            ## 2. candidate breeding marker value
            ## 3. interpretation limits and validation needs
            """

        try :
            return chat_completion_text (
            self .client ,
            self .model_name ,
            [{"role":"user","content":prompt }],
            temperature =0.1 
            )
        except Exception as e :
            return self ._txt (f"> large model reasoningInterfacedisconnect orReport Error: {e }",f"> LLM inference API failed: {e }")
    def summarize_ncrna (self ,mirna_id ,target_ids ,gwas_df ):
        """Summarize ncRNA target GWAS associations without directional inference."""
        if gwas_df .empty :
            return "⚠️ Local GWAS failed to detect significant association signal for the target gene."

            #🌟 Fix: Also replace with to_markdown
        data_markdown =gwas_df .to_csv (sep ='\t',index =False )

        prompt =f"""
        youCurrentlyfor non-Coding RNA '{mirna_id }' Evaluationits yield cost in the field。
        
        【microscopic GWAS empirical evidence】：
        weScandone '{mirna_id }' 's【targetGene】GWAS association signals in natural populations：
        {data_markdown }
        
        【important interpretation rule】：
        Effect_Size(BETA) is retained only as the raw regression coefficient from the GWAS model.
        Do not flip its sign and do not interpret it as positive/negative regulation, promotion/inhibition, increase/decrease, or yield gain/loss.
        Use P_Value and -log10(P) to discuss association strength only. If BETA is mentioned, use only the absolute magnitude.
        
        pleaseOutputreport：
        **1. target trait association summary**
        （briefly describe significant target-gene GWAS associations by P value）
        
        **2. {mirna_id } candidate marker value and validation needs**
        （do not infer phenotype direction; describe only associated traits and the need for validation）
        """
        return chat_completion_text (
        self .client ,
        self .model_name ,
        [{"role":"user","content":prompt }],
        temperature =0.1 
        )
