from __future__ import annotations 

import json 
import math 
import os 
import socket 
from dataclasses import dataclass 
from datetime import datetime 
from pathlib import Path 
from typing import Any 

import pandas as pd 

from src .agents .gwas_agent import GwasAgent 
from src .agents .literature_agent import LiteratureAgent 
from src .agents .macro_network_agent import MacroNetworkAgent 
from src .project_paths import ANNOTATION_DIR ,PROCESSED_DIR 
from src .utils .gene_annotator import RiceGeneAnnotator 


SCENARIO_KEYWORDS ={
"development_architecture":[
"development",
"flower",
"panicle",
"seed",
"embryo",
"organ",
"meristem",
"growth",
],
"hormone_signaling":[
"auxin",
"aba",
"gibberellin",
"ethylene",
"jasmon",
"cytokinin",
"signaling",
"kinase",
],
"membrane_vesicle_transport":[
"transport",
"vesicle",
"clathrin",
"endocyt",
"golgi",
"membrane",
"rab",
"coatomer",
],
"stress_defense":[
"stress",
"defense",
"immune",
"resistance",
"pathogen",
"cold",
"drought",
"salt",
"heat",
],
"metabolism_energy":[
"metabolic",
"biosynthesis",
"catabolic",
"lipid",
"carbon",
"nitrogen",
"energy",
"oxidation",
],
"epigenetic_transcriptional_regulation":[
"transcription",
"chromatin",
"histone",
"methyl",
"acetyl",
"rna polymerase",
"transcription factor",
"dna-binding",
],
}


@dataclass 
class BatchScoringConfig :
    func_weight :float =0.38 
    network_weight :float =0.31 
    gwas_weight :float =0.19 
    completeness_weight :float =0.12 
    top_per_scenario :int =10 
    max_genes :int =300 
    selected_scenarios :list [str ]|None =None 


class ProteinBatchAnalyzer :
    def __init__ (self ,config :BatchScoringConfig |None =None ):
        self .config =config or BatchScoringConfig ()
        self .lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"

        gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")
        map_file =str (ANNOTATION_DIR /"RAP-MSU.txt")
        self .annotator =RiceGeneAnnotator (gff_file ,map_file )
        self .annotator .build_database ()

        self .macro_agent =MacroNetworkAgent ()
        self .gwas_agent =GwasAgent ()
        self .lit_agent =LiteratureAgent ()
        self .enable_literature =os .getenv ("BATCH_ENABLE_LITERATURE","1")=="1"
        self .enable_llm_summary =os .getenv ("BATCH_ENABLE_LLM_SUMMARY","1")=="1"
        self ._network_ok =self ._probe_network ()
        if not self ._network_ok :
            self .enable_literature =False 
            self .enable_llm_summary =False 

        self ._degree_map =self ._build_degree_map ()
        self ._active_scenarios =self ._resolve_scenarios ()

    def _txt (self ,zh :str ,en :str )->str :
        return en if self .lang =="en"else zh 

    def _build_degree_map (self )->dict [str ,int ]:
        df =self .macro_agent .edges_df 
        a =df ["Gene_A"].value_counts ()
        b =df ["Gene_B"].value_counts ()
        degree =(a .add (b ,fill_value =0 )).astype (int )
        return degree .to_dict ()

    def _resolve_scenarios (self )->dict [str ,list [str ]]:
        selected =self .config .selected_scenarios or []
        if not selected :
            return SCENARIO_KEYWORDS 
        keys =[k for k in selected if k in SCENARIO_KEYWORDS ]
        return {k :SCENARIO_KEYWORDS [k ]for k in keys }if keys else SCENARIO_KEYWORDS 

    @staticmethod 
    def _keyword_hit_count (blob :str ,keywords :list [str ])->int :
        text =(blob or "").lower ()
        return sum (1 for kw in keywords if kw in text )

    @staticmethod 
    def _probe_network ()->bool :
        try :
            socket .setdefaulttimeout (3.0 )
            host =socket .gethostbyname ("eutils.ncbi.nlm.nih.gov")
            return bool (host )
        except Exception :
            return False 

    def _normalize_ids (self ,ids :list [str ])->dict [str ,Any ]:
        cleaned =[str (x ).strip ()for x in ids if str (x ).strip ()]
        dedup =list (dict .fromkeys (cleaned ))
        recognized =[]
        unrecognized =[]
        for raw in dedup :
            info =self .annotator .query (raw )
            if info :
                recognized .append (info )
            else :
                unrecognized .append (raw )
        return {
        "input_count":len (ids ),
        "non_empty_count":len (cleaned ),
        "unique_count":len (dedup ),
        "duplicate_count":len (cleaned )-len (dedup ),
        "recognized_count":len (recognized ),
        "unrecognized_count":len (unrecognized ),
        "unrecognized_ids":unrecognized ,
        "recognized_infos":recognized [:self .config .max_genes ],
        "truncated":max (0 ,len (recognized )-self .config .max_genes ),
        }

    @staticmethod 
    def _safe_len (x :Any )->int :
        if x is None :
            return 0 
        if isinstance (x ,(list ,tuple ,set )):
            return len (x )
        if isinstance (x ,str ):
            return 1 if x .strip ()else 0 
        return 0 

    def _scenario_scores (self ,gene_info :dict [str ,Any ])->dict [str ,float ]:
        return self ._scenario_scores_from_evidence (gene_info )

    def _scenario_scores_from_evidence (
    self ,
    gene_info :dict [str ,Any ],
    partner_genes :list [dict [str ,Any ]]|None =None ,
    traits :list [str ]|None =None ,
    literature_boost :float =0.0 ,
    )->dict [str ,float ]:
        desc_blob =" ".join (
        [
        gene_info .get ("Description",""),
        gene_info .get ("Primary_Symbol",""),
        " ".join (gene_info .get ("GO",[])),
        ]
        ).lower ()
        domain_blob =str (gene_info .get ("InterPro","")).lower ()
        partner_blob =" ".join (
        " ".join (
        [
        str (p .get ("Primary_Symbol","")),
        str (p .get ("Note","")),
        " ".join (p .get ("GO",[])),
        ]
        )
        for p in (partner_genes or [])
        ).lower ()
        trait_blob =" ".join (traits or []).lower ()

        scores ={}
        for scenario ,keywords in self ._active_scenarios .items ():
            ann_hits =self ._keyword_hit_count (desc_blob ,keywords )
            dom_hits =self ._keyword_hit_count (domain_blob ,keywords )
            partner_hits =self ._keyword_hit_count (partner_blob ,keywords )
            trait_hits =self ._keyword_hit_count (trait_blob ,keywords )

            ann_score =min (40.0 ,ann_hits *12.0 )
            dom_score =min (20.0 ,dom_hits *10.0 )
            partner_score =min (25.0 ,partner_hits *8.0 )
            trait_score =min (10.0 ,trait_hits *10.0 )
            lit_score =min (10.0 ,literature_boost *0.5 )if (ann_hits +dom_hits +partner_hits +trait_hits )>0 else 0.0 
            scores [scenario ]=round (min (100.0 ,ann_score +dom_score +partner_score +trait_score +lit_score ),3 )
        return scores 

    def _scenario_evidence_flags (
    self ,
    gene_info :dict [str ,Any ],
    partner_genes :list [dict [str ,Any ]],
    traits :list [str ],
    literature_boost :float ,
    scenario :str ,
    )->dict [str ,Any ]:
        keywords =self ._active_scenarios .get (scenario ,[])
        desc_blob =" ".join ([gene_info .get ("Description","")," ".join (gene_info .get ("GO",[]))]).lower ()
        domain_blob =str (gene_info .get ("InterPro","")).lower ()
        partner_blob =" ".join (
        " ".join ([str (p .get ("Primary_Symbol","")),str (p .get ("Note",""))," ".join (p .get ("GO",[]))])
        for p in partner_genes 
        ).lower ()
        trait_blob =" ".join (traits or []).lower ()
        annotation_hits =self ._keyword_hit_count (desc_blob ,keywords )
        domain_hits =self ._keyword_hit_count (domain_blob ,keywords )
        partner_hits =self ._keyword_hit_count (partner_blob ,keywords )
        trait_hits =self ._keyword_hit_count (trait_blob ,keywords )
        flags ={
        "annotation":annotation_hits >0 ,
        "domain_go":domain_hits >0 or bool (gene_info .get ("GO")),
        "coexpression":partner_hits >0 or bool (partner_genes ),
        "gwas_trait":trait_hits >0 or bool (traits ),
        "literature":literature_boost >0 ,
        }
        evidence_labels =[k for k ,v in flags .items ()if v ]
        return {
        "flags":flags ,
        "evidence_labels":evidence_labels ,
        "annotation_hits":annotation_hits ,
        "domain_hits":domain_hits ,
        "partner_hits":partner_hits ,
        "trait_hits":trait_hits ,
        }

    def _func_score (self ,gene_info :dict [str ,Any ])->float :
        go_count =len (gene_info .get ("GO",[]))
        interpro =gene_info .get ("InterPro","")
        base =min (60.0 ,go_count *6.0 )
        if interpro and interpro not in ("N/A","No"):
            base +=30.0 
        if gene_info .get ("Description"):
            base +=10.0 
        return round (min (100.0 ,base ),3 )

    def _network_score (self ,gene_info :dict [str ,Any ])->tuple [float ,list [str ]]:
        rap =gene_info .get ("RAP_ID","")
        msu_candidates =[f"LOC_{x }"if not str (x ).startswith ("LOC_")else str (x )for x in gene_info .get ("MSU_IDs",[])]
        msu_candidates +=[self .macro_agent .rap2msu .get (rap ,"")]
        msu_candidates =[x for x in msu_candidates if x ]
        modules =[]
        degree =0 
        for msu in msu_candidates :
            mods =self .macro_agent .gene2mod .get (msu ,[])
            if mods :
                modules =mods 
                degree =int (self ._degree_map .get (msu ,0 ))
                break 
        if not modules :
            return 0.0 ,[]
        deg_score =min (50.0 ,math .log1p (degree )*10.0 )
        mod_score =min (50.0 ,len (modules )*20.0 )
        return round (mod_score +deg_score ,3 ),modules 

    def _gwas_proxy_score (self ,modules :list [str ])->tuple [float ,list [str ]]:
        if not modules :
            return 0.0 ,[]
        traits =set ()
        for mod in modules :
            traits |=set (self .gwas_agent .module_to_traits .get (mod ,set ()))
        score =min (100.0 ,len (traits )*12.0 +len (modules )*5.0 )
        return round (score ,3 ),sorted (traits )

    def _completeness_score (self ,gene_info :dict [str ,Any ],modules :list [str ])->float :
        score =0.0 
        if gene_info .get ("Description"):
            score +=25 
        if self ._safe_len (gene_info .get ("GO")):
            score +=25 
        if gene_info .get ("InterPro")not in ("","N/A","No",None ):
            score +=25 
        if modules :
            score +=25 
        return float (score )

    def _risk_penalty (self ,gene_info :dict [str ,Any ],modules :list [str ])->tuple [float ,str ]:
        missing =0 
        if not gene_info .get ("Description"):
            missing +=1 
        if not gene_info .get ("GO"):
            missing +=1 
        if gene_info .get ("InterPro")in ("","N/A","No",None ):
            missing +=1 
        if not modules :
            missing +=1 
        penalty =missing *3.0 
        if missing <=1 :
            tag ="low"
        elif missing <=3 :
            tag ="medium"
        else :
            tag ="high"
        return penalty ,tag 

    def _compute_base_rows (self ,recognized_infos :list [dict [str ,Any ]])->list [dict [str ,Any ]]:
        rows =[]
        for info in recognized_infos :
            rap_id =info ["RAP_ID"]
            func_score =self ._func_score (info )
            network_score ,modules =self ._network_score (info )
            gwas_score ,traits =self ._gwas_proxy_score (modules )
            completeness_score =self ._completeness_score (info ,modules )
            risk_penalty ,risk_tag =self ._risk_penalty (info ,modules )
            target_mod ,total_same_mod ,same_mod_found ,partner_genes =self .macro_agent .get_top_hub_partners (
            rap_id ,info ,self .annotator ,top_n =5 
            )
            scenario_scores =self ._scenario_scores_from_evidence (info ,partner_genes =partner_genes ,traits =traits ,literature_boost =0.0 )

            base_score =(
            self .config .func_weight *func_score 
            +self .config .network_weight *network_score 
            +self .config .gwas_weight *gwas_score 
            +self .config .completeness_weight *completeness_score 
            )
            base_score =round (max (0.0 ,min (100.0 ,base_score -risk_penalty )),3 )
            main_scenario =max (scenario_scores ,key =scenario_scores .get )

            rows .append (
            {
            "gene_id":rap_id ,
            "primary_symbol":info .get ("Primary_Symbol",""),
            "description":info .get ("Description",""),
            "interpro":info .get ("InterPro",""),
            "go_terms":info .get ("GO",[]),
            "base_score":base_score ,
            "func_score":func_score ,
            "network_score":network_score ,
            "gwas_score":gwas_score ,
            "completeness_score":completeness_score ,
            "risk_penalty":risk_penalty ,
            "risk_tag":risk_tag ,
            "main_scenario":main_scenario ,
            "scenario_scores":scenario_scores ,
            "modules":modules ,
            "primary_module":target_mod ,
            "same_module_partner_count":total_same_mod ,
            "top_same_module_partner_count":same_mod_found ,
            "partner_genes":partner_genes ,
            "traits":traits ,
            "pmids":[],
            "literature_boost_score":0.0 ,
            }
            )
        return rows 

    def _collect_literature_boost (self ,rows :list [dict [str ,Any ]])->None :
        if not rows or not self .enable_literature :
            return 
        by_scenario :dict [str ,list [dict [str ,Any ]]]={}
        for row in rows :
            by_scenario .setdefault (row ["main_scenario"],[]).append (row )

        for scenario ,subset in by_scenario .items ():
            top_rows =sorted (subset ,key =lambda x :x ["base_score"],reverse =True )[:self .config .top_per_scenario ]
            for row in top_rows :
                gene_info =self .annotator .query (row ["gene_id"])
                if not gene_info :
                    continue 
                try :
                    pmids ,has_cross =self .lit_agent ._search_target_gene_literature (gene_info ,use_symbol_search =True )
                except Exception :
                    pmids ,has_cross =[],False 
                pmids =list (dict .fromkeys (pmids ))[:5 ]
                boost =min (20.0 ,len (pmids )*3.5 +(2.5 if has_cross else 0.0 ))
                row ["literature_boost_score"]=round (boost ,3 )
                row ["pmids"]=pmids 
                gene_info ["Description"]=gene_info .get ("Description","")
                row ["scenario_scores"]=self ._scenario_scores_from_evidence (
                gene_info ,
                partner_genes =row .get ("partner_genes",[]),
                traits =row .get ("traits",[]),
                literature_boost =row ["literature_boost_score"],
                )

    def _build_rank_df (self ,rows :list [dict [str ,Any ]])->pd .DataFrame :
        if not rows :
            return pd .DataFrame ()
        for row in rows :
            total =row ["base_score"]+row .get ("literature_boost_score",0.0 )
            row ["total_score"]=round (min (100.0 ,total ),3 )
        rank_df =pd .DataFrame (rows )
        rank_df =rank_df .sort_values ("total_score",ascending =False ).reset_index (drop =True )
        rank_df ["rank"]=rank_df .index +1 
        return rank_df 

    def _build_scenario_df (self ,rank_df :pd .DataFrame )->pd .DataFrame :
        if rank_df .empty :
            return pd .DataFrame ()
        records =[]
        for scenario ,sub in rank_df .groupby ("main_scenario"):
            top =sub .sort_values ("total_score",ascending =False ).head (10 )
            records .append (
            {
            "scenario":scenario ,
            "gene_count":len (sub ),
            "scenario_score":round (float (sub ["total_score"].mean ()),3 ),
            "representative_genes":", ".join (top ["gene_id"].tolist ()),
            "support_pmids":", ".join (sorted ({pid for p in top ["pmids"]for pid in p })[:20 ]),
            }
            )
        return pd .DataFrame (records ).sort_values ("scenario_score",ascending =False )

    def _build_assignment_matrix_df (self ,rank_df :pd .DataFrame )->pd .DataFrame :
        if rank_df .empty :
            return pd .DataFrame ()
        records =[]
        scenario_names =list (self ._active_scenarios .keys ())
        for _ ,row in rank_df .iterrows ():
            record ={
            "gene_id":row ["gene_id"],
            "primary_symbol":row ["primary_symbol"],
            "main_scenario":row ["main_scenario"],
            "total_score":row ["total_score"],
            }
            for scen in scenario_names :
                record [scen ]=float (row ["scenario_scores"].get (scen ,0.0 ))
            records .append (record )
        matrix_df =pd .DataFrame (records )
        return matrix_df .sort_values ("total_score",ascending =False ).reset_index (drop =True )

    def _assignment_tier (self ,assignment_score :float )->str :
        if assignment_score >=70 :
            return "Tier 1"
        if assignment_score >=45 :
            return "Tier 2"
        if assignment_score >=25 :
            return "Tier 3"
        return "Background"

    def _build_selected_scenario_tables (self ,rank_df :pd .DataFrame )->dict [str ,pd .DataFrame ]:
        if rank_df .empty :
            return {}
        out :dict [str ,pd .DataFrame ]={}
        for scenario in self ._active_scenarios .keys ():
            records =[]
            for _ ,row in rank_df .iterrows ():
                assignment_score =float (row ["scenario_scores"].get (scenario ,0.0 ))
                if assignment_score <25 :
                    continue 
                evidence =self ._scenario_evidence_flags (
                {
                "Description":row .get ("description",""),
                "InterPro":row .get ("interpro",""),
                "GO":row .get ("go_terms",[]),
                },
                row .get ("partner_genes",[]),
                row .get ("traits",[]),
                float (row .get ("literature_boost_score",0.0 )),
                scenario ,
                )
                records .append (
                {
                "gene_id":row ["gene_id"],
                "primary_symbol":row ["primary_symbol"],
                "assignment_score":round (assignment_score ,3 ),
                "priority_score":round (float (row ["total_score"]),3 ),
                "tier":self ._assignment_tier (assignment_score ),
                "evidence_types":", ".join (evidence ["evidence_labels"])if evidence ["evidence_labels"]else "annotation",
                "traits":", ".join (row .get ("traits",[]))if row .get ("traits")else "",
                "modules":", ".join (row .get ("modules",[]))if row .get ("modules")else "",
                "pmid_count":len (row .get ("pmids",[])),
                "risk_tag":row .get ("risk_tag",""),
                "description":row .get ("description",""),
                }
                )
            out [scenario ]=(
            pd .DataFrame (records ).sort_values (["assignment_score","priority_score"],ascending =False ).reset_index (drop =True )
            if records 
            else pd .DataFrame ()
            )
        return out 

    def _build_gene_cards (self ,rank_df :pd .DataFrame )->list [dict [str ,Any ]]:
        cards =[]
        for _ ,row in rank_df .iterrows ():
            scen =row ["main_scenario"]
            cards .append (
            {
            "gene_id":row ["gene_id"],
            "symbol":row ["primary_symbol"],
            "rank":int (row ["rank"]),
            "total_score":float (row ["total_score"]),
            "main_scenario":scen ,
            "sub_scores":{
            "func":float (row ["func_score"]),
            "network":float (row ["network_score"]),
            "gwas_proxy":float (row ["gwas_score"]),
            "completeness":float (row ["completeness_score"]),
            "literature_boost":float (row ["literature_boost_score"]),
            "risk_penalty":float (row ["risk_penalty"]),
            },
            "risk_tag":row ["risk_tag"],
            "traits":row ["traits"],
            "modules":row ["modules"],
            "partner_genes":row .get ("partner_genes",[]),
            "pmids":row ["pmids"],
            "evidence_note":self ._txt (
            "This score is used for experimental priority screening, not causal conclusion.",
            "This score is for experimental prioritization, not causal proof.",
            ),
            }
            )
        return cards 

    def _cluster_name_label (self ,scenario :str )->str :
        label_map ={
        "development_architecture":"Development and Architecture",
        "hormone_signaling":"Hormone Signaling",
        "membrane_vesicle_transport":"Membrane and Vesicle Transport",
        "stress_defense":"Stress and Defense",
        "metabolism_energy":"Metabolism and Energy",
        "epigenetic_transcriptional_regulation":"Epigenetic and Transcriptional Regulation",
        }
        return label_map .get (scenario ,scenario .replace ("_"," ").title ())

    def _build_deg_clusters (self ,rank_df :pd .DataFrame )->list [dict [str ,Any ]]:
        if rank_df .empty :
            return []

        clusters =[]
        for scenario in self ._active_scenarios .keys ():
            sub =rank_df [rank_df ["main_scenario"]==scenario ].copy ()
            if sub .empty :
                continue 
            sub =sub .sort_values (["total_score","base_score"],ascending =False )
            gene_infos =[]
            for _ ,row in sub .head (max (self .config .top_per_scenario ,12 )).iterrows ():
                info =self .annotator .query (row ["gene_id"])
                if not info :
                    continue 
                gene_infos .append (info )
            if not gene_infos :
                continue 
            clusters .append (
            {
            "cluster_key":scenario ,
            "cluster_name":self ._cluster_name_label (scenario ),
            "gene_infos":gene_infos ,
            }
            )
        return clusters 

    def _run_deg_cluster_summaries (self ,rank_df :pd .DataFrame )->tuple [str ,list [dict [str ,Any ]]]:
        if rank_df .empty :
            return self ._txt ("No summary result available.","No analyzable results."),[]
        if not self .enable_llm_summary :
            return "LLM summary skipped (network unavailable or disabled).",[]

        reports =[]
        rendered_sections =[]
        for cluster in self ._build_deg_clusters (rank_df ):
            cluster_name =cluster ["cluster_name"]
            gene_infos =cluster ["gene_infos"]
            abstracts_text =self .lit_agent .fetch_abstracts_route_c (gene_infos )
            summary =self .lit_agent .summarize_deg_cluster (cluster_name ,gene_infos ,abstracts_text )
            pmid_count =abstracts_text .count ("- PMID:")
            member_labels =[g .get ("Primary_Symbol")or g .get ("RAP_ID")for g in gene_infos ]
            reports .append (
            {
            "cluster_key":cluster ["cluster_key"],
            "cluster_name":cluster_name ,
            "gene_count":len (gene_infos ),
            "genes":member_labels ,
            "abstracts_pmid_count":pmid_count ,
            "summary":summary ,
            }
            )
            rendered_sections .append (
            "\n\n".join (
            [
            f"## {cluster_name }",
            f"Cluster genes ({len (member_labels )}): {', '.join (member_labels [:20 ])}",
            f"Rice literature records retrieved: {pmid_count }",
            summary ,
            ]
            )
            )
        if not rendered_sections :
            return self ._txt ("No summary result available.","No analyzable results."),reports 
        return "\n\n---\n\n".join (rendered_sections ),reports 

    def _run_builtin_benchmark (self ,rank_df :pd .DataFrame )->dict [str ,Any ]:
        """
        Optional offline benchmark:
        If data/processed/batch_benchmark_sets.json exists, evaluate top-k hit rate.
        JSON format:
        { "stress_defense": ["Os01g...", ...], ... }
        """
        bench_path =PROCESSED_DIR /"batch_benchmark_sets.json"
        if not bench_path .exists ()or rank_df .empty :
            return {"available":False ,"message":"benchmark set not found"}
        try :
            benchmark =json .loads (bench_path .read_text (encoding ="utf-8"))
            top20 =set (rank_df .head (20 )["gene_id"].tolist ())
            hit_map ={}
            for name ,genes in benchmark .items ():
                gs =set (genes )
                if not gs :
                    continue 
                hit =len (top20 &gs )
                hit_map [name ]={"top20_hit":hit ,"top20_recall":round (hit /len (gs ),4 )}
            return {"available":True ,"metrics":hit_map }
        except Exception as exc :
            return {"available":False ,"message":str (exc )}

    def analyze_gene_ids (self ,gene_ids :list [str ])->dict [str ,Any ]:
        qc =self ._normalize_ids (gene_ids )
        recognized_infos =qc ["recognized_infos"]
        rows =self ._compute_base_rows (recognized_infos )
        self ._collect_literature_boost (rows )
        rank_df =self ._build_rank_df (rows )
        scenario_df =self ._build_scenario_df (rank_df )
        assignment_matrix_df =self ._build_assignment_matrix_df (rank_df )
        selected_scenario_tables =self ._build_selected_scenario_tables (rank_df )
        gene_cards =self ._build_gene_cards (rank_df )
        llm_summary ,cluster_reports =self ._run_deg_cluster_summaries (rank_df )
        benchmark =self ._run_builtin_benchmark (rank_df )
        return {
        "qc_report":qc ,
        "rank_df":rank_df ,
        "scenario_df":scenario_df ,
        "assignment_matrix_df":assignment_matrix_df ,
        "selected_scenario_tables":selected_scenario_tables ,
        "gene_cards":gene_cards ,
        "llm_summary":llm_summary ,
        "cluster_reports":cluster_reports ,
        "benchmark":benchmark ,
        }

    def analyze_from_csv (self ,input_csv :str ,gene_column :str ="gene_id")->dict [str ,Any ]:
        df =pd .read_csv (input_csv )
        if gene_column not in df .columns :
            raise ValueError (f"Column '{gene_column }' not found. Available columns: {', '.join (df .columns )}")
        gene_ids =[str (x )for x in df [gene_column ].tolist ()]
        return self .analyze_gene_ids (gene_ids )

    def dump_results (self ,result :dict [str ,Any ],output_dir :str |None =None )->dict [str ,str ]:
        ts =datetime .now ().strftime ("%Y%m%d_%H%M%S")
        root =Path (output_dir )if output_dir else (PROCESSED_DIR /"batch_results"/ts )
        root .mkdir (parents =True ,exist_ok =True )

        rank_path =root /"rank_table.csv"
        scenario_path =root /"scenario_table.csv"
        cards_path =root /"gene_cards.json"
        clusters_path =root /"cluster_reports.json"
        qc_path =root /"qc_report.json"
        llm_path =root /"llm_summary.md"
        bench_path =root /"benchmark_report.json"

        rank_df :pd .DataFrame =result ["rank_df"]
        scenario_df :pd .DataFrame =result ["scenario_df"]
        rank_export =rank_df .drop (columns =["scenario_scores"],errors ="ignore")
        rank_export .to_csv (rank_path ,index =False )
        scenario_df .to_csv (scenario_path ,index =False )
        cards_path .write_text (json .dumps (result ["gene_cards"],ensure_ascii =False ,indent =2 ),encoding ="utf-8")
        clusters_path .write_text (json .dumps (result .get ("cluster_reports",[]),ensure_ascii =False ,indent =2 ),encoding ="utf-8")
        qc_path .write_text (json .dumps (result ["qc_report"],ensure_ascii =False ,indent =2 ),encoding ="utf-8")
        llm_path .write_text (result ["llm_summary"],encoding ="utf-8")
        bench_path .write_text (json .dumps (result .get ("benchmark",{}),ensure_ascii =False ,indent =2 ),encoding ="utf-8")

        return {
        "output_dir":str (root ),
        "rank_table":str (rank_path ),
        "scenario_table":str (scenario_path ),
        "gene_cards":str (cards_path ),
        "cluster_reports":str (clusters_path ),
        "qc_report":str (qc_path ),
        "llm_summary":str (llm_path ),
        "benchmark_report":str (bench_path ),
        }
