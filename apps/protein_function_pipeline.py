import argparse 
import base64 
import json 
import math 
import os 
import sys 
import re 
from pathlib import Path 
from collections import Counter 
from rich .console import Console 
from rich .panel import Panel 
from rich .markdown import Markdown 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .utils .gene_annotator import RiceGeneAnnotator 
from src .utils .epigenomic_analyzer import EpigenomicAnalyzer 
from src .utils .expression_pattern_service import ExpressionPatternService 
from src .utils .protein_feature_service import ProteinFeatureService 
from src .agents .macro_network_agent import MacroNetworkAgent 
from src .agents .literature_agent import LiteratureAgent 
from src .agents .gwas_agent import GwasAgent 
from src .project_paths import ANNOTATION_DIR 
from src .utils .protein_batch_analysis import SCENARIO_KEYWORDS 
from src .utils .llm_client import chat_completion_text 

PROTEIN_JSON_BEGIN ="PROTEIN_RESULT_JSON_BEGIN"
PROTEIN_JSON_END ="PROTEIN_RESULT_JSON_END"
PROTEIN_STEP_BEGIN ="PROTEIN_STEP_JSON_BEGIN"
PROTEIN_STEP_END ="PROTEIN_STEP_JSON_END"
PROTEIN_PROGRESS_PREFIX ="PROTEIN_PROGRESS_JSON:"

class EnhancedFunctionPredictor :
    def __init__ (self ):
        self .console =Console ()
        self .lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"

        #Initialize basic components
        self .gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")
        self .map_file =str (ANNOTATION_DIR /"RAP-MSU.txt")
        self .annotator =RiceGeneAnnotator (self .gff_file ,self .map_file )
        self .annotator .build_database ()

        self .macro_agent =MacroNetworkAgent ()
        self .lit_agent =LiteratureAgent ()
        self .epigenomic_analyzer =EpigenomicAnalyzer ()
        self .expression_pattern_service =ExpressionPatternService ()
        self .gwas_agent =GwasAgent ()
        self .protein_feature_service =ProteinFeatureService ()

    def _txt (self ,zh ,en ):
        return en if self .lang =="en"else zh 

    def _emit_step (self ,step ,title_zh ,title_en ,data ):
        """Output incremental step JSON for real-time dashboard rendering"""
        payload ={
        "step":step ,
        "title":self ._txt (title_zh ,title_en ),
        "status":"done",
        "data":data ,
        }
        print (PROTEIN_STEP_BEGIN )
        print (json .dumps (payload ,ensure_ascii =False ))
        print (PROTEIN_STEP_END )
        sys .stdout .flush ()

    def _emit_progress (self ,scope ,data ):
        payload ={"scope":scope ,**data }
        print (PROTEIN_PROGRESS_PREFIX +json .dumps (payload ,ensure_ascii =False ))
        sys .stdout .flush ()

    def _llm_call (self ,prompt ):
        """Unified LLM invocation entry point"""
        if os .getenv ("RIGOR_LLM_ENABLED","1")!="1"or not os .getenv ("LLM_API_KEY","").strip ():
            return ""
        try :
            return chat_completion_text (
            self .lit_agent .client ,
            model =self .lit_agent .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
        except Exception as e :
            self .console .print (self ._txt (f"LLM inferFailed: {e }",f"LLM inference failed: {e }"))
            return ""

    def _infer_homology (self ,gene_info ,abstracts ):
        """Infer gene function based on homologous literature"""
        if not abstracts or not abstracts .strip ():
            return self ._txt ("No homologous gene literature found; unable to perform inference in this dimension.","No homology literature found for this dimension.")
        target =gene_info ["RAP_ID"]
        note =gene_info .get ("Description","N/A")
        interpro =gene_info .get ("InterPro","N/A")
        if self .lang =="en":
            prompt =f"""You are a rigorous rice molecular biologist. Based ONLY on the homology literature below, infer the biological function of gene '{target }'.

[Gene info]
- RAP ID: {target }
- Description: {note }
- InterPro: {interpro }

[Homology literature]
{abstracts [:100000 ]}

Rules:
- Prioritize the gene's own Description (e.g. "Similar to XXX") as the strongest evidence.
- If cross-species literature is cited, clearly state the species.
- Output 3-5 sentences. Cite PMIDs where applicable.
- Do NOT fabricate titles, journals, or PMIDs.
"""
        else :
            prompt =f"""You are a rigorous expert in rice molecular biology。please**just**based on the following homologsGeneLiterature，InferenceGene '{target }' biologicalFunction。

【GeneInformation】
- RAP ID: {target }
- FunctionDescription: {note }
- InterProDomain: {interpro }

【homologousGeneLiterature】
{abstracts [:100000 ]}

【Inferencerule】
- highestPriority：Functionhomology describedInformation（if "Similar to XXX"）
- ifHascross-speciesLiterature，explicitly label species
- Output 3-5 sentence，reference PMID
- strictly prohibit fabricatingLiteraturetitle、journal name、PMID
"""
        self .console .print (self ._txt ("Homologous literature LLM inference in progress...","Homology literature LLM inference..."))
        return self ._llm_call (prompt )

    def _infer_coexpression (self ,gene_info ,part2_data ,abstracts ):
        """Infer Gene Function Based on Co-expression Network and Literature Independent Inference"""
        target =gene_info ["RAP_ID"]
        note =gene_info .get ("Description","N/A")
        interpro =part2_data .get ("interpro","N/A")
        go_enrichment =part2_data .get ("go_enrichment",[])
        go_text ="\n".join ([f"- {go } ({self ._txt (f'appear{c }time',f'count={c }')})"for go ,c in go_enrichment [:10 ]])
        if not abstracts or not abstracts .strip ():
            if not go_text :
                return self ._txt ("No co-expressed genes or literature, and no GO enrichment data.","No co-expression literature and no GO enrichment data.")
        if self .lang =="en":
            prompt =f"""You are a rigorous rice molecular biologist. Based on the GO enrichment analysis and co-expression gene literature below, infer the function of gene '{target }'.

[Gene info]
- RAP ID: {target }
- Description: {note }
- InterPro: {interpro }

[GO enrichment (Top 10)]
{go_text or "No significant GO terms"}

[Co-expression gene literature]
{abstracts [:100000 ]if abstracts else "No co-expression literature"}

Rules:
- First summarize partner gene functions, then infer target gene function.
- Use Description and InterPro as primary evidence; co-expression as supporting.
- Output 3-5 sentences. Cite PMIDs where applicable.
- Do NOT fabricate titles, journals, or PMIDs.
"""
        else :
            prompt =f"""You are a rigorous expert in rice molecular biology。based on the following GO EnrichmentAnalysiswith commonExpressioncomradeGeneLiterature，InferenceGene '{target }' 'sFunction。

【GeneInformation】
- RAP ID: {target }
- FunctionDescription: {note }
- InterProDomain: {interpro }

【GO EnrichmentAnalysis (Top 10)】
{go_text or "No significant GO terms found."}

【totalExpressionGeneLiterature】
{abstracts [:100000 ]if abstracts else "No co-expressed genes or literature."}

【Inferencerule】
- first briefly inventory comradesGene'sFunction，againInferencetargetGeneFunction
- Functiondescription and InterPro forPrimaryEvidence，totalExpressionforAuxiliaryValidation
- Output 3-5 sentence，reference PMID
- strictly prohibit fabricatingLiteraturetitle、journal name、PMID
"""
        self .console .print (self ._txt ("Co-expression network LLM reasoning in progress...","Co-expression network LLM inference..."))
        return self ._llm_call (prompt )

    def _infer_epigenomic (self ,gene_info ,part4_data ):
        """Infer Regulation Status Based on Start Sub-epigenetic Regulation Data Independent Inference"""
        epigenomic_summary =part4_data .get ("summary_text","")
        if not epigenomic_summary or not epigenomic_summary .strip ():
            return self ._txt ("No start sub-epigenomic evidence found.","No promoter epigenomic evidence.")
        target =gene_info ["RAP_ID"]
        if self .lang =="en":
            prompt =f"""You are a rice epigenetics expert. Based on the promoter peak scan data below, assess the regulatory state of gene '{target }' across tissues.

[Promoter epigenomic data]
{epigenomic_summary }

[Interpretation guide]
- ATAC/FAIRE/MNase open peaks → chromatin accessibility
- H3K4me3, H3K27ac, H3K9ac → active promoter marks
- H3K27me3, H3K9me2 → repressive marks
- Co-existence of active and repressive marks → bivalent/dynamic state

Output: In which tissues is the promoter likely active, repressed, or bivalent? 2-3 sentences.
"""
        else :
            prompt =f"""You are an expert in rice epigenetics。please refer to the followingStartatom peak ScanData，JudgmentGene '{target }' differentTissueinRegulationstate。

【Startsub-epigenomeData】
{epigenomic_summary }

【explanation prompt】
- ATAC/FAIRE/MNase open peak → chromatin open
- H3K4me3, H3K27ac, H3K9ac → ActiveStartsub marker
- H3K27me3, H3K9me2 → Inhibitionmarker
- ActivitywithInhibitionmarkers coexist → bivalent/dynamic state

Output：thisGeneStartin which subTissuemoderately active、biasInhibitionor bivalent？2-3 sentence。
"""
        self .console .print (self ._txt ("Epigenetic regulation LLM reasoning in progress...","Epigenomic LLM inference..."))
        return self ._llm_call (prompt )

    def _split_report_sections (self ,report_text ):
        if not report_text :
            return []
        sections =[]
        current_title =None 
        current_lines =[]

        def flush ():
            nonlocal current_title ,current_lines 
            body ="\n".join (current_lines ).strip ()
            if current_title or body :
                sections .append ({
                "title":current_title or self ._txt ("Unnamed Module","Untitled Section"),
                "body":body ,
                })
            current_title =None 
            current_lines =[]

        for raw_line in report_text .splitlines ():
            line =raw_line .strip ()
            title =None 
            if line .startswith ("**")and line .endswith ("**")and len (line )>4 :
                title =line .strip ("*").strip ()
            elif re .match (r"^#{1,3}\s+",line ):
                title =re .sub (r"^#{1,3}\s+","",line ).strip ()
            elif re .match (r"^\d+\s+\S+",line ):
                title =line 
            elif line .lower ().startswith ("final conclusion"):
                title =line 
            elif line .startswith ("Final Conclusion"):
                title =line 

            if title :
                flush ()
                current_title =title 
            else :
                current_lines .append (raw_line )
        flush ()
        return [item for item in sections if item .get ("title")or item .get ("body")]

    def _infer_theme (self ,text_blob ):
        text =(text_blob or "").lower ()
        best_name ="development_architecture"
        best_score =-1 
        for name ,keywords in SCENARIO_KEYWORDS .items ():
            score =sum (1 for kw in keywords if kw in text )
            if score >best_score :
                best_name =name 
                best_score =score 
        return best_name 

    def _build_visualization_payload (self ,gene_info ,part1_abstracts ,part2_data ,part3_abstracts ,part4_data ,part5_data ,part6_data ):
        refs1 =self ._extract_references_from_abstracts (part1_abstracts )
        refs3 =self ._extract_references_from_abstracts (part3_abstracts )
        tissue_scores =part4_data .get ("tissue_scores",[])
        gwas_df =part6_data ["gwas_df"]if part6_data is not None else None 
        target_gene =gene_info .get ("RAP_ID","")
        target_gwas_df =gwas_df [gwas_df ["Gene"]==target_gene ].copy ()if gwas_df is not None and not gwas_df .empty and "Gene"in gwas_df .columns else None 
        if target_gwas_df is not None and not target_gwas_df .empty :
            target_gwas_df ["neg_log10_p"]=target_gwas_df ["P_Value"].apply (lambda x :round (-math .log10 (max (float (x ),1e-300 )),4 ))
            core_target_gwas_df =target_gwas_df [target_gwas_df ["Trait_Normalized_Evidence"].isin (["Strong","Moderate"])]if "Trait_Normalized_Evidence"in target_gwas_df .columns else target_gwas_df 
            target_gwas_records =core_target_gwas_df .to_dict ("records")
        else :
            target_gwas_records =[]

        lit_count =len ({item ["pmid"]for item in refs1 +refs3 })
        interpro =gene_info .get ("InterPro","")
        go_count =len (gene_info .get ("GO",[]))
        network_genes =part2_data .get ("network_genes",[])
        epigen_hit_count =len (part4_data .get ("chromatin_accessibility",[]))+len (part4_data .get ("histone_modification",[]))
        top_tissue_score =max ([float (item .get ("score",0.0 ))for item in tissue_scores ],default =0.0 )
        target_trait_count =len (target_gwas_records )
        gwas_level_weights ={"Strong":20.0 ,"Moderate":12.0 ,"Weak":4.0 }
        gwas_score =sum (gwas_level_weights .get (str (item .get ("Trait_Normalized_Evidence","Weak")),4.0 )for item in target_gwas_records )

        radar_scores ={
        "literature":round (min (100.0 ,lit_count *12.0 ),3 ),
        "domain_go":round (min (100.0 ,(30.0 if interpro not in ("","N/A","No")else 0.0 )+min (70.0 ,go_count *8.0 )),3 ),
        "coexpression":round (min (100.0 ,len (network_genes )*7.5 +(25.0 if part3_abstracts .strip ()else 0.0 )),3 ),
        "epigenomics":round (min (100.0 ,epigen_hit_count *4.0 +top_tissue_score *2.0 ),3 ),
        "gwas":round (min (100.0 ,gwas_score ),3 ),
        }

        heatmap_rows =[]
        for item in tissue_scores :
            heatmap_rows .append ({
            "tissue":item .get ("tissue","unknown"),
            "open_count":len (item .get ("open_evidence",[])),
            "active_count":len (item .get ("active_marks",[])),
            "repressive_count":len (item .get ("repressive_marks",[])),
            "poised_count":len (item .get ("poised_marks",[])),
            "context_count":len (item .get ("context_marks",[])),
            "score":float (item .get ("score",0.0 )),
            "state":item .get ("state","unknown"),
            "open_marks":item .get ("open_evidence",[]),
            "active_marks":item .get ("active_marks",[]),
            "repressive_marks":item .get ("repressive_marks",[]),
            })

        network_nodes =[{
        "id":target_gene ,
        "label":gene_info .get ("Primary_Symbol")or target_gene ,
        "role":"target",
        "theme":self ._infer_theme (" ".join ([gene_info .get ("Description",""),interpro ," ".join (gene_info .get ("GO",[]))])),
        "correlation":1.0 ,
        }]
        network_edges =[]
        for friend in network_genes :
            label =friend .get ("Primary_Symbol")or friend .get ("RAP_ID","")
            theme =self ._infer_theme (" ".join ([friend .get ("Note","")," ".join (friend .get ("GO",[])),label ]))
            corr =float (friend .get ("Correlation",0.0 ))
            network_nodes .append ({
            "id":friend .get ("RAP_ID",""),
            "label":label ,
            "role":"partner",
            "theme":theme ,
            "correlation":corr ,
            })
            network_edges .append ({
            "source":target_gene ,
            "target":friend .get ("RAP_ID",""),
            "weight":corr ,
            })

        return {
        "evidence_radar":{
        "scores":radar_scores ,
        "meta":{
        "literature_count":lit_count ,
        "go_count":go_count ,
        "network_gene_count":len (network_genes ),
        "epigenomic_hit_count":epigen_hit_count ,
        "target_trait_count":target_trait_count ,
        },
        },
        "promoter_heatmap":heatmap_rows ,
        "trait_effect_plot":target_gwas_records ,
        "coexpression_network":{
        "nodes":network_nodes ,
        "edges":network_edges ,
        },
        }

    def _emit_structured_payload (self ,gene_info ,is_coordinate_mode ,final_report ,gwas_report =None ,visualizations =None ):
        payload ={
        "report_type":"coordinate"if is_coordinate_mode else "protein_gene",
        "target":gene_info .get ("RAP_ID",""),
        "symbol":gene_info .get ("Primary_Symbol",""),
        "description":gene_info .get ("Description",""),
        "location":gene_info .get ("Location",""),
        "final_sections":self ._split_report_sections (final_report ),
        "gwas_sections":self ._split_report_sections (gwas_report or ""),
        "visualizations":visualizations or {},
        }
        print (PROTEIN_JSON_BEGIN )
        print (json .dumps (payload ,ensure_ascii =False ))
        print (PROTEIN_JSON_END )

    def part0_protein_features (self ,gene_info ):
        """Supplement protein features: UniProt subcellular localization, InterPro/Pfam domains, and plotting data via web query"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Protein Feature Supplement Annotation','Protein Feature Annotation')}[/bold cyan]")

        uniprot_id =gene_info .get ("UniProt_ID","")
        transcript_ids =gene_info .get ("Transcript_IDs",[])
        if transcript_ids :
            self .console .print (self ._txt (f"Transcriptnumber: {len (transcript_ids )}",f"Transcript count: {len (transcript_ids )}"))
            self .console .print (self ._txt (f"TranscriptList: {', '.join (transcript_ids )}",f"Transcript IDs: {', '.join (transcript_ids )}"))
        if not uniprot_id :
            self .console .print (self ._txt ("Local annotation does not provide UniProt ID; skipping web-based protein feature query","UniProt ID missing in local annotation. Skip online protein feature query."))
            return {}

        self .console .print (f"UniProt ID: {uniprot_id }")
        features =self .protein_feature_service .collect (gene_info )

        if features .get ("protein_length"):
            self .console .print (self ._txt (f"protein length: {features ['protein_length']} aa",f"Protein length: {features ['protein_length']} aa"))
        if features .get ("subcellular_locations"):
            self .console .print (self ._txt (f"subcellular localization: {', '.join (features ['subcellular_locations'])}",f"Subcellular localization: {', '.join (features ['subcellular_locations'])}"))
        if features .get ("interpro_domains"):
            self .console .print (self ._txt (f"InterPro Domainsegment: {len (features ['interpro_domains'])} each",f"InterPro domain fragments: {len (features ['interpro_domains'])}"))
        if features .get ("pfam_domains"):
            self .console .print (self ._txt (f"Pfam drawDomain: {len (features ['pfam_domains'])} each",f"Pfam plotting domains: {len (features ['pfam_domains'])}"))
        for warning in features .get ("warnings",[]):
            self .console .print (f"{warning }")

        return features 

    def part1_homology_literature (self ,gene_info ):
        """Part 1: Homologous Gene Literature Analysis (Rice → Plants → Mammals)"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Part 1: Homologous Gene Literature Analysis','Part 1: Homology Literature Analysis')}[/bold cyan]")

        #Extract keywords
        note =gene_info .get ('Description','')
        primary_sym =gene_info .get ('Primary_Symbol','')

        search_terms =[]
        if note and len (note .strip ())>5 :
            llm_keywords =self .lit_agent ._extract_protein_keywords_with_llm (note )
            if llm_keywords :
                self .console .print (self ._txt (f"LLMidentify keywords: {', '.join (llm_keywords )}",f"LLM keywords: {', '.join (llm_keywords )}"))
                search_terms .extend (llm_keywords )

        if primary_sym and len (primary_sym )>2 :
            search_terms .append (primary_sym )

        if not search_terms :
            self .console .print (self ._txt ("Unable to extract effective keywords; skipping homologous literature search","No valid keyword extracted. Skip homology literature search."))
            return ""

        all_pmids =[]
        rice_pmids =[]
        plant_pmids =[]
        mammal_pmids =[]

        #1. Search rice literature (up to 5 articles)
        self .console .print (self ._txt ("Currently searching for homologous gene literature in rice...","Searching rice homology literature..."))

        #Prioritize direct search using Gene Symbol
        if primary_sym and len (primary_sym )>2 :
            query =f'"{primary_sym }"[Title/Abstract] AND (rice[Title/Abstract] OR Oryza[Title/Abstract])'
            self .console .print (self ._txt (f"directlySearchGenesymbol: {primary_sym }",f"Search by primary symbol: {primary_sym }"))
            try :
                rice_pmids .extend (self .lit_agent .search_pubmed_ids (query ,retmax =5 ,sort ="relevance"))
            except :
                pass 

                #Then use keywords extracted by LLM for search
        for term in search_terms [:2 ]:
            if term ==primary_sym :#Avoid duplicate searches
                continue 
            query =f'"{term }"[Title/Abstract] AND (rice[Title/Abstract] OR Oryza[Title/Abstract])'
            try :
                rice_pmids .extend (self .lit_agent .search_pubmed_ids (query ,retmax =5 ,sort ="relevance"))
            except :
                pass 

        rice_pmids =list (set (rice_pmids ))[:5 ]
        all_pmids .extend (rice_pmids )

        #2. If rice literature is insufficient, search other plants (10 articles)
        if len (all_pmids )<5 :
            self .console .print (self ._txt ("Insufficient rice literature; currently searching for homologous genes in other plants...","Rice literature is insufficient. Expanding to other plant homologs..."))
            main_keyword =search_terms [0 ]
            query =f'"{main_keyword }"[Title/Abstract] AND (Arabidopsis[Title/Abstract] OR maize[Title/Abstract] OR wheat[Title/Abstract] OR plant[Title/Abstract])'
            try :
                plant_pmids =self .lit_agent .search_pubmed_ids (query ,retmax =10 ,sort ="relevance")
                all_pmids .extend (plant_pmids )
            except :
                pass 

                #3. If plant literature is also insufficient, search mammals (10 articles)
        if len (all_pmids )<5 :
            self .console .print (self ._txt ("Insufficient plant literature; currently searching for homologous genes in mammals...","Plant literature is insufficient. Expanding to mammalian homologs..."))
            main_keyword =search_terms [0 ]
            query =f'"{main_keyword }"[Title/Abstract] AND (human[Title/Abstract] OR mouse[Title/Abstract] OR mammal[Title/Abstract])'
            try :
                mammal_pmids =self .lit_agent .search_pubmed_ids (query ,retmax =10 ,sort ="relevance")
                all_pmids .extend (mammal_pmids )
            except :
                pass 

                #Deduplicate and fetch literature
        all_pmids =list (set (all_pmids ))[:15 ]

        if not all_pmids :
            self .console .print (self ._txt ("Homologous gene literature not found","No homology literature found."))
            return ""

        self .console .print (self ._txt (f"Successfulfound {len (all_pmids )} homologous articleGeneLiterature",f"Found {len (all_pmids )} homology papers."))
        abstracts =self .lit_agent ._fetch_medline_by_ids (all_pmids )

        return abstracts 

    def part2_domain_and_go_enrichment (self ,gene_info ,target_gene ):
        """Part 2: Domain Analysis + GO Enrichment Analysis"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Part 2: Domain and GO Enrichment Analysis','Part 2: Domain and GO Enrichment')}[/bold cyan]")

        # 2.1 DomainAnalysis
        interpro =gene_info .get ('InterPro','No')
        self .console .print (self ._txt (f"InterProDomain: {interpro }",f"InterPro domains: {interpro }"))

        #2.2 Find 100 co-expressed genes (for GO enrichment, use as many as available).
        self .console .print (self ._txt ("Currently searching for co-expressed genes...","Searching co-expression genes..."))
        _ ,_ ,_ ,network_genes =self .macro_agent .get_top_hub_partners (
        target_gene ,gene_info ,self .annotator ,top_n =100 
        )

        if not network_genes :
            self .console .print (self ._txt ("No co-expressed genes found; using only the target gene's own GO annotations.","No co-expression genes found. Use target gene GO only."))
            all_go_terms =gene_info .get ('GO',[])
        else :
            self .console .print (self ._txt (f"found {len (network_genes )} commonExpressionGene",f"Found {len (network_genes )} co-expression genes."))

            #2.3 Perform GO enrichment analysis on 11 genes.
            all_go_terms =gene_info .get ('GO',[]).copy ()
            for g in network_genes :
                all_go_terms .extend (g .get ('GO',[]))

                #GO Enrichment Statistics
        go_counter =Counter (all_go_terms )
        top_go =go_counter .most_common (10 )

        self .console .print (self ._txt ("GO Enrichment Analysis (Top 10):","GO enrichment (Top 10):"))
        for go_term ,count in top_go :
            self .console .print (self ._txt (f"    - {go_term } (appear{count }time)",f"    - {go_term } (count={count })"))

        return {
        'interpro':interpro ,
        'network_genes':network_genes ,
        'go_enrichment':top_go 
        }

    def part3_coexpression_literature (self ,network_genes ,gene_info ):
        """Part Three: Co-expression Gene Literature Analysis"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Part Three: Co-expression Gene Literature Analysis','Part 3: Co-expression Literature Analysis')}[/bold cyan]")

        if not network_genes :
            self .console .print (self ._txt ("No co-expressed genes; skipping this section.","No co-expression genes. Skip this part."))
            return ""

            #Using existing literature search logic.
        def _progress_callback (
        done ,
        total ,
        gene_data ,
        retrieved_for_gene ,
        cumulative ,
        maximum ,
        rate_limited =False ,
        rate_limit_stage ="",
        rate_limit_detail ="",
        ):
            current_gene =gene_data or {}
            self ._emit_progress ("coexpression_literature",{
            "done":done ,
            "total":total ,
            "gene_id":current_gene .get ("RAP_ID",""),
            "symbol":current_gene .get ("Primary_Symbol",""),
            "retrieved_for_gene":retrieved_for_gene ,
            "cumulative_literature":cumulative ,
            "max_literature":maximum ,
            "rate_limited":rate_limited ,
            "rate_limit_stage":rate_limit_stage ,
            "rate_limit_detail":rate_limit_detail ,
            })

        abstracts =self .lit_agent .fetch_abstracts_route_b (
        network_genes ,
        target_gene_info =gene_info ,
        progress_callback =_progress_callback ,
        )

        if abstracts .strip ():
            self .console .print (f"SuccessfulFetchtotalExpressionGeneLiterature")
        else :
            self .console .print (self ._txt ("No literature found for co-expressed genes.","No co-expression literature found."))

        return abstracts 

    def part4_epigenomic_promoter_evidence (self ,gene_info =None ,coord =None ,strand =None ):
        """Part Four: Promoter Epigenetic Regulation Evidence"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Part Four: Promoter Epigenetic Regulation Evidence','Part 4: Promoter Epigenomic Evidence')}[/bold cyan]")

        epigenomic_data =self .epigenomic_analyzer .analyze (
        gene_info =gene_info ,
        coord =coord ,
        strand =strand 
        )
        promoter =epigenomic_data ["promoter_region"]

        self .console .print (self ._txt (
        f"Startsubregion: {promoter ['promoter_label']} (TSS={promoter ['tss']}, strand={promoter ['strand']})",
        f"Promoter region: {promoter ['promoter_label']} (TSS={promoter ['tss']}, strand={promoter ['strand']})"
        ))
        self .console .print (self ._txt (
        f"chromatin openness: {len (epigenomic_data ['chromatin_accessibility'])} eachDatacollect",
        f"Chromatin accessibility datasets: {len (epigenomic_data ['chromatin_accessibility'])}"
        ))
        self .console .print (self ._txt (
        f"histone modification: {len (epigenomic_data ['histone_modification'])} eachDatacollect",
        f"Histone modification datasets: {len (epigenomic_data ['histone_modification'])}"
        ))

        return epigenomic_data 

    def part5_expression_patterns (self ,gene_info ):
        """Part Five: Multi-Transcript Expression Pattern Analysis"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Part Five: Multi-Transcript Expression Patterns','Part 5: Multi-transcript Expression Patterns')}[/bold cyan]")
        expression_data =self .expression_pattern_service .collect (gene_info )
        if expression_data .get ("available"):
            count =len (expression_data .get ('transcript_results',[]))
            self .console .print (self ._txt (f"Hasgenerate {count } eachTranscript/model'sExpressionmode summary",f"Generated expression summaries for {count } transcript/model entries."))
        else :
            self .console .print (self ._txt ("No usable expression pattern summary generated at this time.","No usable expression summary generated."))
        for warning in expression_data .get ("warnings",[]):
            self .console .print (f"{warning }")
        return expression_data 

    def part6_gwas_trait_assessment (self ,gene_info ,network_genes ):
        """Part Six: GWAS Trait and Breeding Evaluation"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Part Six: GWAS Traits and Breeding Evaluation','Part 6: GWAS Trait Assessment')}[/bold cyan]")

        target_gene =gene_info ["RAP_ID"]
        network_gene_ids =[gene ["RAP_ID"]for gene in network_genes ]if network_genes else []
        gwas_status =self .gwas_agent .inspect_gwas_status ()

        macro_evidence =self .gwas_agent .get_macro_evidence (target_gene ,network_gene_ids ,self .annotator )
        self .console .print (self ._txt ("Macro-level trait evidence compiled.","Macro-module trait evidence prepared."))
        self .console .print (
        self ._txt (
        f"GWAS Directory: {gwas_status ['gwas_dir']} | Hasdiscover {gwas_status ['existing_count']}/{gwas_status ['expected_count']} standardResultFile",
        f"GWAS directory: {gwas_status ['gwas_dir']} | Found {gwas_status ['existing_count']}/{gwas_status ['expected_count']} standard result files"
        )
        )

        all_scan_targets =list (dict .fromkeys ([target_gene ]+network_gene_ids ))
        gwas_df =self .gwas_agent .scan_network_gwas (all_scan_targets ,self .annotator )

        if gwas_df .empty :
            self .console .print (self ._txt ("No significant GWAS quantitative signals detected, or local GWAS file missing.","No significant GWAS signals found, or local GWAS files are missing."))
        else :
            self .console .print (self ._txt (f"Scanreached {len (gwas_df )} significant GWAS Signal",f"Detected {len (gwas_df )} significant GWAS signals."))

        return {
        "macro_evidence":macro_evidence ,
        "gwas_df":gwas_df ,
        "scan_targets":all_scan_targets ,
        "gwas_status":gwas_status ,
        }

    def _extract_references_from_abstracts (self ,abstracts_text ):
        if not abstracts_text :
            return []
        refs =[]
        blocks =[blk .strip ()for blk in abstracts_text .split ("\n\n")if "PMID:"in blk and "Title:"in blk ]
        for block in blocks :
            pmid_match =block .split ("PMID:",1 )[1 ].split ("\n",1 )[0 ].strip ()if "PMID:"in block else ""
            title_match =block .split ("Title:",1 )[1 ].split ("\n",1 )[0 ].strip ()if "Title:"in block else "N/A"
            journal_match =block .split ("Journal:",1 )[1 ].split ("\n",1 )[0 ].strip ()if "Journal:"in block else "N/A"
            journal_match =re .sub (r"\s*\(DOI:.*?\)\s*$","",journal_match ).strip ()
            if pmid_match :
                refs .append ({"pmid":pmid_match ,"title":title_match ,"journal":journal_match })
        unique ={}
        for item in refs :
            unique [item ["pmid"]]=item 
        return list (unique .values ())

    def _format_reference_list (self ,refs ):
        if not refs :
            return "No citable references"if self .lang =="en"else "No citable literature."
        if self .lang =="en":
            return "\n".join ([
            f"- *{item ['title']} - {item ['journal']} (PMID: {item ['pmid']})*"
            for item in refs 
            ])
        return "\n".join ([
        f"- *《{item ['title']}》 - {item ['journal']} (PMID: {item ['pmid']})*"
        for item in refs 
        ])

    def synthesize_results (self ,gene_info ,part1_abstracts ,part2_data ,part3_abstracts ,part4_data ,part5_data ):
        """Synthesize functional evidence to generate final prediction report."""
        self .console .print (f"\n[bold cyan]{self ._txt ('Comprehensive Analysis and Functional Prediction','Integrated Analysis and Function Prediction')}[/bold cyan]")

        target_gene =gene_info ['RAP_ID']
        note =gene_info .get ('Description','No')
        interpro =part2_data ['interpro']
        go_enrichment =part2_data ['go_enrichment']
        epigenomic_summary =part4_data .get ('summary_text',self ._txt ('No promoter epigenetic evidence.','No promoter epigenomic evidence'))
        expression_summary =part5_data .get ('summary_text',self ._txt ('No expression pattern evidence.','No expression pattern evidence'))
        refs =self ._extract_references_from_abstracts (part1_abstracts )+self ._extract_references_from_abstracts (part3_abstracts )
        refs ={item ["pmid"]:item for item in refs }.values ()
        references_for_prompt =self ._format_reference_list (list (refs ))

        #Constructing prompt...
        if self .lang =="en":
            prompt =f"""You are a rigorous rice molecular biologist. Infer the biological function of rice gene '{target_gene }' based on five independent evidence sources, and assess potential promoter regulatory states across tissues.

[Gene info]
- RAP ID: {target_gene }
- Description: {note }
- InterPro: {interpro }

[Evidence 1: Homology literature, highest priority]
{part1_abstracts [:100000 ]if part1_abstracts else "No homology literature"}

[Evidence 2: GO enrichment]
{chr (10 ).join ([f"- {go } (count={count })"for go ,count in go_enrichment [:10 ]])}

[Evidence 3: Co-expression literature]
{part3_abstracts [:100000 ]if part3_abstracts else "No co-expression literature"}

[Evidence 4: Promoter epigenomic evidence]
Peak scans in promoter region (TSS upstream 2000 bp, downstream 500 bp):
{epigenomic_summary }

[Evidence 5: Multi-transcript expression]
{expression_summary }

[Citable references from system retrieval]
{references_for_prompt }

Output format:
1) Function inference from homology literature
2) Function inference from domain and GO enrichment
3) Function inference from co-expression network
4) Promoter regulatory-state judgment from epigenomic evidence
5) Judgment from multi-transcript expression patterns
Final conclusion in 2-4 sentences.

Do not create an extra references section inside the generated text. Do not fabricate title, journal, or PMID.
"""
        else :
            prompt =f"""You are a rigorous expert in rice molecular biology。Please base your analysis on the following five independent evidence sources，comprehensiveInferencericeGene '{target_gene }' biologicalFunction，and judge itsStartsub in differentTissuepotential inRegulationstate。

【GeneBasicInformation】
- RAP ID: {target_gene }
- FunctionDescription: {note }
- InterProDomain: {interpro }

【Evidence1：homologousGeneLiterature（highestPriority）】
following is theGeneStudies of homologous proteins in different speciesLiterature：
{part1_abstracts [:100000 ]if part1_abstracts else "No homologous gene literature."}

【Evidence2：GOEnrichmentAnalysis】
based on targetGeneits10commonExpressionGene'sGOAnnotationstatistics，discover following high-frequencyGO terms：
{chr (10 ).join ([f"- {go } (appear{count }time)"for go ,count in go_enrichment [:10 ]])}

【Evidence3：totalExpressionGeneLiterature】
following are target-relatedGenetotalExpressioncomradeGene'sFunctionLiterature：
{part3_abstracts [:100000 ]if part3_abstracts else "No co-expressed genes or literature."}

【Evidence4：Startsub epigeneticRegulationEvidence】
following is based on targetGene TSS upstream 2000 bp、downstream 500 bp Startsubregion of peak ScanResult：
{epigenomic_summary }

【Evidence5：manyTranscriptExpressionpattern】
following are targetsGenein organ/period and differentTreatmentConditionunderExpressionmode summary：
{expression_summary }

【citableLiteraturelist（SystemSearch）】
{references_for_prompt }

【Inferencerule】
1. **highestPriority**：Evidence1（homologousGeneLiterature）- this is the most directFunctionEvidence
2. **timePriority**：Evidence2（Domain+GOEnrichment）- provide molecularFunctionstatistical support
 3. **AuxiliaryValidation**：Evidence3（totalExpressionGene）- used forValidationand supplement the aboveInference
 4. **epigeneticRegulationlayer**：Evidence4（promoter peak）- used to judge thisGenedifferentTissuepotential activation in、Inhibitionor dynamicRegulationstate
 5. **Expressionpattern**：Evidence5（organ/Period + TreatmentExpressionpattern）- used to judge differentTranscript/Expressionmodel in organ、period andTreatmentConditionunderExpressiondifference

【epigenome interpretation hints】
- ATAC/FAIRE/MNase Open peaks generally support open chromatin
- H3K4me3、H3K27ac、H3K9ac、H3K4ac usually support activityStartatom
- H3K27me3、H3K9me2、H3K9me1 commonly indicatesInhibitionstate
- ifActivitywithInhibitionmarkers coexist，should be consideredTissuespecificity、stage differences or bivalent/dynamic state
- promoter state score higher，more support thisTissuemiddle promoter preferential activation；if state for dynamic_or_bivalent，explain different markers coexist

please follow the followingStructureOutputreport：

**Result1：homology-basedGene'sFunctionInference**
evidence-based1'sLiterature，SummarythisGeneHomologous proteins in other speciesHasknownFunction。ifHascross-speciesLiterature，please clearly label species。

**Result2：basedDomainandGOenrichedFunctionInference**
CombinationInterProDomainandGOEnrichmentAnalysis，from molecularFunctionangleInferencethisGenePossible biological processes involved。

**Result3：based on co-ExpressionNetwork'sFunctionInference**
combining evidence3，from co-Expressioncomrade'sFunctionInferencetargetGenepotentialFunction。byResult1forPrimarybase comprehensive judgment on。

**Result4：basedStartevidence from offspring epigenomeRegulationstate judgment**
combining evidence4，this indicatesGene promoter in whichTissuemore like open activation、Inhibitionsilence or dynamic switch state，and explain supporting your peak type。

**Result5：based on multipleTranscriptExpressionpattern judgment**
SummarythisGenedifferentTranscript/Expressionmodel in organ、period andTreatmentConditionunderExpressionpattern，point out whether there is obviousTissuePreference、Differences in developmental stage preference or stress response。

**finalConclusion**
integrating fiveResult，gives thisGeneMost probable biologicalFunction，and supplement its possibleTissuespecificityRegulationcharacteristics andExpressionmode features（2-4sentence）。

【reference rule】
can mention in main text PMID as argument source，But do not create additional content at the end of the main text“referenceLiterature/Appendix”chapter。
strictly prohibit fabricatingLiteraturetitle、journal name、PMID。
"""

        try :
            if os .getenv ("RIGOR_LLM_ENABLED","1")!="1"or not os .getenv ("LLM_API_KEY","").strip ():
                return ""
            self .console .print (self ._txt ("Currently performing comprehensive reasoning by large language model...","Running integrated LLM inference..."))
            llm_text =chat_completion_text (
            self .lit_agent .client ,
            model =self .lit_agent .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
            tail_title ="## References (System Retrieved)"if self .lang =="en"else "## Reference Literature (System Search)"
            return llm_text +f"\n\n---\n\n{tail_title }\n"+self ._format_reference_list (list (refs ))
        except Exception as e :
            self .console .print (self ._txt (f"large model reasoningFailed: {e }",f"LLM inference failed: {e }"))
            return ""

    def synthesize_coordinate_results (self ,coord ,strand ,part4_data ):
        """In coordinate mode, only perform regional-level promoter epigenetic interpretation, do not pretend to predict gene function."""
        self .console .print (f"\n[bold cyan]{self ._txt ('Regional-Level Epigenetic Regulation Interpretation','Regional Epigenomic Interpretation')}[/bold cyan]")

        promoter =part4_data ["promoter_region"]
        tissue_scores =part4_data .get ("tissue_scores",[])
        epigenomic_summary =part4_data .get ("summary_text",self ._txt ("No promoter epigenetic evidence.","No promoter epigenomic evidence"))

        if self .lang =="en":
            prompt =f"""You are a rigorous rice epigenetics expert. The input is a genomic coordinate region, not a confirmed gene.

[Input region]
- Original coordinate: {coord }
- Strand: {strand or '+'}
- Estimated TSS: {promoter ['tss']}
- Promoter region: {promoter ['promoter_label']}

[Constraints]
1. No gene ID/domain/GO/co-expression/literature is provided.
2. Do not treat this region as a known gene for functional prediction.
3. Only infer potential regulatory state from promoter-region epigenomic peaks.
4. Do not fabricate gene names, literature, PMIDs, domains, or expression conclusions.

[Promoter epigenomic evidence]
{epigenomic_summary }

[Tissue score summary]
{chr (10 ).join ([f"- tissue={item ['tissue']}, score={item ['score']}, state={item ['state']}"for item in tissue_scores ])if tissue_scores else "- No tissue score"}

Output format:
1) Region-level epigenomic observations
2) Region-level regulatory-state judgment
3) What can and cannot be concluded
Final conclusion in 2-4 sentences. Subject must be "this region" or "this promoter region".
"""
        else :
            prompt =f"""You are a rigorous expert in rice epigenetics。now given to you is notHasknownGene，but rather aGenegroup coordinate region。

【Inputregion】
- original coordinates: {coord }
- chain direction: {strand or '+'}
- predicted TSS: {promoter ['tss']}
- promoter region: {promoter ['promoter_label']}

【important limitation】
1. here noHasprovide truthGene ID、GeneAnnotation、Domain、GO、totalExpressionorLiteratureInformation。
2. You cannot directly treat this coordinate region asHasknownGeneTo doFunctionpredict。
3. you can only base on promoter regional epigenome peak，judge the differences in this region acrossTissuepotential inRegulationstate。
4. strictly prohibit fabricatingGenename、Literature、PMID、DomainorExpressionConclusion。

【promoter epigenomic evidence】
{epigenomic_summary }

【Tissuescoring summary】
{chr (10 ).join ([f"- tissue={item ['tissue']}, score={item ['score']}, state={item ['state']}"for item in tissue_scores ])if tissue_scores else "- No tissue score"}

Please strictly follow the followingStructureOutput：

**Result1：regional level epigenome observation**
ObjectiveSummarythis promoter regions in differentTissueOpen chromatin and histone modification features。

**Result2：regional levelRegulationstate judgment**
judge the differences in this region acrossTissuemore like in active、poised_or_weakly_active、repressed Or dynamic_or_bivalent，and explain basis。

**Result3：can and cannot be drawnConclusion**
clarify which onesConclusioncan be derived from this regionDataget，whichConclusioncannot obtain。
for example：can determine promoter Activitytrend；cannot directly assert specificGeneFunction。

**finalConclusion**
use 2-4 sentenceSummaryThe most reasonable explanation for this coordinate region。subject must be“this region”or“this promoter region”，cannot be written as“thisGene”。
"""

        try :
            if os .getenv ("RIGOR_LLM_ENABLED","1")!="1"or not os .getenv ("LLM_API_KEY","").strip ():
                return ""
            self .console .print (self ._txt ("Currently, regional-level epigenetic regulation is being explained by the large language model...","Running LLM inference for regional epigenomic interpretation..."))
            return chat_completion_text (
            self .lit_agent .client ,
            model =self .lit_agent .model_name ,
            messages =[{"role":"user","content":prompt }],
            temperature =0.1 
            )
        except Exception as e :
            self .console .print (self ._txt (f"regional level large model reasoningFailed: {e }",f"Regional LLM inference failed: {e }"))
            return ""

    def render_gwas_results (self ,gene_info ,part6_data ):
        """Output GWAS trait and breeding evaluation report"""
        self .console .print (f"\n[bold cyan]{self ._txt ('Comprehensive GWAS result evaluation','GWAS Summary Assessment')}[/bold cyan]")

        target_gene =gene_info ["RAP_ID"]
        gwas_df =part6_data ["gwas_df"]
        macro_evidence =part6_data ["macro_evidence"]
        gwas_status =part6_data ["gwas_status"]

        if gwas_df .empty :
            self .console .print (self ._txt ("No significant GWAS quantitative traits generated; attempting to output a simplified explanation","No significant GWAS quantitative table generated. Fallback summary will be used."))
        else :
            self .console .print (self ._txt (f"will {len (gwas_df )} line GWAS quantitativeResultgive to large modelSummary",f"Submitting {len (gwas_df )} GWAS quantitative records for LLM summarization."))

        if gwas_status ["existing_count"]==0 :
            self .console .print (self ._txt ("Currently only retaining the GWAS analysis framework; once server results are mounted or synchronized, full field evaluation will be automatically generated","Only GWAS analysis framework is currently available. Full field assessment requires mounted/synced server results."))

        return self .gwas_agent .summarize (
        target_gene =target_gene ,
        target_ids =part6_data ["scan_targets"],
        gwas_df =gwas_df [gwas_df ["Trait_Normalized_Evidence"].isin (["Strong","Moderate"])]if gwas_df is not None and not gwas_df .empty and "Trait_Normalized_Evidence"in gwas_df .columns else gwas_df ,
        macro_evidence =macro_evidence 
        )

    def run (self ,target_gene =None ,coord =None ,strand =None ):
        """Running step-by-step analysis process to generate results"""
        self .console .print (Panel (
        self ._txt (
        f"[bold cyan]GeneFunctionfour-dimensional predictionSystem[/bold cyan]\nInputtarget: {target_gene or coord }",
        f"[bold cyan]Four-dimensional Gene Function Prediction[/bold cyan]\nInput target: {target_gene or coord }",
        ),
        border_style ="cyan"
        ))

        #── Phase 1: Query Gene Information ──
        self .console .print (self ._txt ("\n[bold green][Phase 1] Gene Information query[/bold green]","\n[bold green][Phase 1] Gene Information Query[/bold green]"))
        gene_info =None 
        if target_gene :
            gene_info =self .annotator .query (target_gene )

        if target_gene and not gene_info :
            self .console .print (self ._txt (f"[bold red]Not FoundGene '{target_gene }'[/bold red]",f"[bold red]Gene '{target_gene }' not found[/bold red]"))
            sys .exit (1 )

        is_coordinate_mode =bool (coord and not target_gene )

        if not gene_info and coord :
            gene_info ={
            'RAP_ID':coord ,
            'Primary_Symbol':'N/A',
            'Description':self ._txt ('Annotation not provided; currently performing coordinate-based start sub-epigenome analysis','No annotation provided. Coordinate mode runs promoter epigenomic analysis only.'),
            'InterPro':'N/A',
            'GO':[]
            }

            #── Step 1: Gene Information + protein features (combined output) ──
        protein_features =self .part0_protein_features (gene_info )if target_gene else {}
        plot_b64 =""
        if protein_features .get ("plot_png"):
            raw =protein_features ["plot_png"]
            if isinstance (raw ,bytes ):
                plot_b64 =base64 .b64encode (raw ).decode ("ascii")
            elif hasattr (raw ,"getvalue"):
                plot_b64 =base64 .b64encode (raw .getvalue ()).decode ("ascii")
        self ._emit_step ("basic_info","Gene information and protein features","Gene Info & Protein Features",{
        "target":gene_info .get ("RAP_ID",""),
        "symbol":gene_info .get ("Primary_Symbol","N/A"),
        "description":gene_info .get ("Description","N/A"),
        "location":gene_info .get ("Location","N/A"),
        "uniprot":gene_info .get ("UniProt_ID",""),
        "protein_length":protein_features .get ("protein_length"),
        "subcellular_locations":protein_features .get ("subcellular_locations",[]),
        "interpro_domains":protein_features .get ("interpro_domains",[]),
        "pfam_domains":protein_features .get ("pfam_domains",[]),
        "plot_png_b64":plot_b64 ,
        })

        #Collecting reports from each step for final assembly
        step_reports =[]

        if is_coordinate_mode :
        #── Coordinate mode: only performing epigenetic regulation ──
            part4_data =self .part4_epigenomic_promoter_evidence (coord =coord ,strand =strand )
            final_report =self .synthesize_coordinate_results (coord ,strand ,part4_data )
            self ._emit_step ("epigenomic","Start sub-epigenetic regulation","Promoter Epigenomics",{
            "report":final_report ,
            "tissue_scores":[
            {"tissue":ts ["tissue"],"score":ts .get ("score",0 ),"state":ts .get ("state","")}
            for ts in part4_data .get ("tissue_scores",[])
            ],
            "heatmap_data":self ._build_promoter_heatmap_data (part4_data ),
            })
            step_reports .append (final_report )
            gwas_report =""
            part1_abstracts =""
            part2_data ={'interpro':'N/A','network_genes':[],'go_enrichment':[]}
            part3_abstracts =""
            part5_data ={}
            part6_data =None 
        else :
        #── Step 2: Homologous literature → instant LLM summary ──
            has_symbol =gene_info .get ('Primary_Symbol','')!='N/A'
            if target_gene and has_symbol :
                part1_abstracts =self .part1_homology_literature (gene_info )
                homology_report =self ._infer_homology (gene_info ,part1_abstracts )
            else :
                part1_abstracts =""
                homology_report =self ._txt ("This gene has no known symbol; skipping homologous literature analysis.","Gene has no known symbol. Skipping homology literature analysis.")
            self ._emit_step ("homology","Homologous literature function inference","Homology Literature Inference",{
            "report":homology_report ,
            "references":self ._parse_references (part1_abstracts ),
            })
            step_reports .append (homology_report )

            #── Step 3: GO enrichment and co-expression network first; companion-gene literature is deferred. ──
            if target_gene and has_symbol :
                part2_data =self .part2_domain_and_go_enrichment (gene_info ,target_gene )
                part3_abstracts =""
                coexp_report =""
            else :
                part2_data ={'interpro':gene_info .get ('InterPro','N/A'),'network_genes':[],'go_enrichment':[]}
                part3_abstracts =""
                coexp_report =self ._txt ("This gene has no known symbol; skipping co-expression analysis.","Gene has no known symbol. Skipping co-expression analysis.")

                #Constructing co-expression network visualization data
            network_nodes =[]
            network_edges =[]
            for ng in part2_data .get ("network_genes",[]):
                nid =ng .get ("RAP_ID","")
                corr =float (ng .get ("Correlation",0.0 ))
                network_nodes .append ({
                "id":nid ,
                "label":ng .get ("Primary_Symbol")or nid ,
                "role":"partner",
                "theme":self ._infer_theme (" ".join ([ng .get ("Note","")," ".join (ng .get ("GO",[])),ng .get ("Primary_Symbol","")])),
                "correlation":corr ,
                })
                network_edges .append ({
                "source":gene_info .get ("RAP_ID",""),
                "target":nid ,
                "weight":corr ,
                })
            if network_nodes :
                network_nodes .insert (0 ,{
                "id":gene_info .get ("RAP_ID",""),
                "label":gene_info .get ("Primary_Symbol")or gene_info .get ("RAP_ID",""),
                "role":"target",
                "theme":self ._infer_theme (" ".join ([gene_info .get ("Description",""),gene_info .get ("InterPro","")," ".join (gene_info .get ("GO",[]))])),
                "correlation":1.0 ,
                })

            self ._emit_step ("coexpression","Co-expression network function inference","Co-expression Network Inference",{
            "pending":True ,
            "message":self ._txt (
            "共表达伙伴文献检索通常耗时较长。RIGOR 会先继续完成表观组、表达模式和 GWAS 模块，随后回到本模块补全文献推理结果，请耐心等待。",
            "Co-expression companion-gene literature retrieval can take longer than other modules. RIGOR will first continue with epigenomics, expression, and GWAS, then return here with the literature-based network inference. Please wait patiently."
            ),
            "report":"",
            "references":[],
            "go_terms":[{"term":t ,"count":c }for t ,c in part2_data .get ("go_enrichment",[])[:10 ]],
            "network_gene_count":len (part2_data .get ("network_genes",[])),
            "network_viz":{"nodes":network_nodes ,"edges":network_edges },
            })

            #── Step 4: Start sub-epigenetic regulation → instant LLM summary ──
            part4_data =self .part4_epigenomic_promoter_evidence (gene_info =gene_info )
            epigenomic_report =self ._infer_epigenomic (gene_info ,part4_data )
            self ._emit_step ("epigenomic","Start sub-epigenetic regulation","Promoter Epigenomics",{
            "report":epigenomic_report ,
            "tissue_scores":[
            {"tissue":ts ["tissue"],"score":ts .get ("score",0 ),"state":ts .get ("state","")}
            for ts in part4_data .get ("tissue_scores",[])
            ],
            "heatmap_data":self ._build_promoter_heatmap_data (part4_data ),
            })
            step_reports .append (epigenomic_report )

            #── Step 5: Multi-transcript expression pattern (no LLM) ──
            part5_data =self .part5_expression_patterns (gene_info )if target_gene else {}
            expr_summary =part5_data .get ("summary_text",self ._txt ("No expression pattern data.","No expression pattern data."))
            #Coding expression chart as base64 sent to frontend
            expr_transcripts =[]
            for tr in part5_data .get ("transcript_results",[]):
                entry ={
                "transcript_id":tr .get ("transcript_id",""),
                "expression_model_id":tr .get ("expression_model_id",""),
                "developmental_summary":tr .get ("developmental_summary",[]),
                "treatment_summary":tr .get ("treatment_summary",[]),
                "warnings":tr .get ("warnings",[]),
                "developmental_plot_b64":"",
                "treatment_plot_b64":"",
                }
                for key_src ,key_dst in [("developmental_plot_png","developmental_plot_b64"),
                ("treatment_plot_png","treatment_plot_b64")]:
                    raw =tr .get (key_src )
                    if raw :
                        if isinstance (raw ,bytes ):
                            entry [key_dst ]=base64 .b64encode (raw ).decode ("ascii")
                        elif hasattr (raw ,"getvalue"):
                            entry [key_dst ]=base64 .b64encode (raw .getvalue ()).decode ("ascii")
                expr_transcripts .append (entry )
            self ._emit_step ("expression","Multi-transcript expression pattern","Expression Patterns",{
            "summary":expr_summary ,
            "transcript_ids":part5_data .get ("transcript_ids",[]),
            "expression_model_ids":part5_data .get ("expression_model_ids",[]),
            "transcript_results":expr_transcripts ,
            "warnings":part5_data .get ("warnings",[]),
            })
            step_reports .append (expr_summary )

            # ── Step 6: GWAS Evaluation → gwas_agent.summarize ──
            part6_data =self .part6_gwas_trait_assessment (gene_info ,part2_data ["network_genes"])if target_gene else None 
            if part6_data is not None :
                gwas_report =self .render_gwas_results (gene_info ,part6_data )
                gwas_df =part6_data .get ("gwas_df")
                target_gene_id =gene_info .get ("RAP_ID","")
                target_gwas_df =gwas_df [gwas_df ["Gene"]==target_gene_id ].copy ()if gwas_df is not None and not gwas_df .empty and "Gene"in gwas_df .columns else None 
                trait_effect_data =[]
                suggestive_trait_effect_data =[]
                if target_gwas_df is not None and not target_gwas_df .empty :
                    target_gwas_df ["neg_log10_p"]=target_gwas_df ["P_Value"].apply (lambda x :round (-math .log10 (max (float (x ),1e-300 )),4 ))
                    if "Trait_Normalized_Evidence"in target_gwas_df .columns :
                        trait_effect_data =target_gwas_df [target_gwas_df ["Trait_Normalized_Evidence"].isin (["Strong","Moderate"])].to_dict ("records")
                        suggestive_trait_effect_data =target_gwas_df [target_gwas_df ["Trait_Normalized_Evidence"]=="Weak"].to_dict ("records")
                    else :
                        trait_effect_data =target_gwas_df .to_dict ("records")
                self ._emit_step ("gwas","GWAS trait and breeding evaluation","GWAS Trait & Breeding Assessment",{
                "report_sections":self ._split_report_sections (gwas_report ),
                "trait_effect_data":trait_effect_data ,
                "suggestive_trait_effect_data":suggestive_trait_effect_data ,
                })
            else :
                gwas_report =""

            #── Deferred Step 3: Companion-gene literature retrieval and co-expression interpretation ──
            if target_gene and has_symbol :
                part3_abstracts =self .part3_coexpression_literature (
                part2_data ['network_genes'],gene_info 
                )
                coexp_report =self ._infer_coexpression (gene_info ,part2_data ,part3_abstracts )

            self ._emit_step ("coexpression","Co-expression network function inference","Co-expression Network Inference",{
            "report":coexp_report ,
            "references":self ._parse_references (part3_abstracts ),
            "go_terms":[{"term":t ,"count":c }for t ,c in part2_data .get ("go_enrichment",[])[:10 ]],
            "network_gene_count":len (part2_data .get ("network_genes",[])),
            "network_viz":{"nodes":network_nodes ,"edges":network_edges },
            })

                #── Final: Output complete JSON (maintain compatibility for review) ──
        if not is_coordinate_mode and coexp_report :
            step_reports .insert (1 ,coexp_report )
        final_report ="\n\n---\n\n".join (step_reports )
        visualizations ={}
        if not is_coordinate_mode :
            visualizations =self ._build_visualization_payload (
            gene_info ,part1_abstracts ,part2_data ,part3_abstracts ,
            part4_data ,part5_data ,part6_data ,
            )
        self ._emit_structured_payload (gene_info ,is_coordinate_mode ,final_report ,gwas_report ,visualizations )
        self .console .print (self ._txt ("[bold green]AnalysisCompleted。[/bold green]","[bold green]Analysis completed.[/bold green]"))

    def _parse_references (self ,abstracts_text ):
        """Extracting structured literature list from Medline text"""
        if not abstracts_text :
            return []
        refs =[]
        for block in re .split (r'\n(?=- PMID:)',abstracts_text .strip ()):
            if not block .strip ():
                continue 
            pmid_m =re .search (r'PMID:\s*(\d+)',block )
            title_m =re .search (r'Title:\s*(.+)',block )
            journal_m =re .search (r'Journal:\s*(.+)',block )
            if pmid_m :
                raw_journal =journal_m .group (1 ).strip ()if journal_m else ""
                clean_journal =re .sub (r"\s*\(DOI:[^)]*\)\s*$","",raw_journal ).strip ()
                refs .append ({
                "pmid":pmid_m .group (1 ),
                "title":title_m .group (1 ).strip ()if title_m else "",
                "journal":clean_journal ,
                })
        return refs 

    def _build_promoter_heatmap_data (self ,part4_data ):
        """Extracting start sub-heatmap data from part4_data"""
        tissue_scores =part4_data .get ("tissue_scores",[])
        return [
        {
        "tissue":ts .get ("tissue",""),
        "open_count":len (ts .get ("open_evidence",[])),
        "active_count":len (ts .get ("active_marks",[])),
        "repressive_count":len (ts .get ("repressive_marks",[])),
        "open_marks":ts .get ("open_evidence",[]),
        "active_marks":ts .get ("active_marks",[]),
        "repressive_marks":ts .get ("repressive_marks",[]),
        "score":ts .get ("score",0 ),
        "state":ts .get ("state",""),
        }
        for ts in tissue_scores 
        ]


if __name__ =="__main__":
    parser =argparse .ArgumentParser (description ="GeneFunction Four-Dimensional Prediction System")
    parser .add_argument ("--gene",type =str ,help ="Target Gene ID")
    parser .add_argument ("--coord",type =str ,help ="Target Gene coordinates, e.g., chr04:18172160-18177058")
    parser .add_argument ("--strand",type =str ,choices =["+","-"],help ="Strand direction in coordinate mode")
    parser .add_argument ("--llm-base-url",type =str ,help ="LLM Interface address")
    parser .add_argument ("--llm-api-key",type =str ,help ="LLM API Key")
    parser .add_argument ("--pubmed-max-date",type =str ,help ="PubMed Literature cutoff date, format YYYY/MM/DD, e.g., 2025/12/31")
    args =parser .parse_args ()

    if not args .gene and not args .coord :
        parser .error ("At least one of --gene or --coord must be provided")

    if args .llm_base_url :
        os .environ ["LLM_BASE_URL"]=args .llm_base_url 
    if args .llm_api_key :
        os .environ ["LLM_API_KEY"]=args .llm_api_key 
    if args .pubmed_max_date :
        os .environ ["PUBMED_MAX_DATE"]=args .pubmed_max_date 

    predictor =EnhancedFunctionPredictor ()
    predictor .run (target_gene =args .gene ,coord =args .coord ,strand =args .strand )
