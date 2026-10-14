#!/usr/bin/env python3

import argparse 
import json 
import os 
import sys 
from pathlib import Path 

from rich .console import Console 
from rich .markdown import Markdown 
from rich .panel import Panel 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .utils .protein_batch_analysis import BatchScoringConfig ,ProteinBatchAnalyzer 

BATCH_JSON_BEGIN ="BATCH_RESULT_JSON_BEGIN"
BATCH_JSON_END ="BATCH_RESULT_JSON_END"


def parse_weights (weights :str |None ,default :BatchScoringConfig )->BatchScoringConfig :
    if not weights :
        return default 
    try :
        payload =json .loads (weights )
    except Exception as exc :
        raise ValueError (f"Invalid --weights JSON: {exc }")from exc 
    return BatchScoringConfig (
    func_weight =float (payload .get ("func_weight",default .func_weight )),
    network_weight =float (payload .get ("network_weight",default .network_weight )),
    gwas_weight =float (payload .get ("gwas_weight",default .gwas_weight )),
    completeness_weight =float (payload .get ("completeness_weight",default .completeness_weight )),
    top_per_scenario =int (payload .get ("top_per_scenario",default .top_per_scenario )),
    max_genes =int (payload .get ("max_genes",default .max_genes )),
    selected_scenarios =payload .get ("selected_scenarios",default .selected_scenarios ),
    )


def main ():
    parser =argparse .ArgumentParser (description ="Batch protein-coding gene analysis (scenario + prioritization)")
    parser .add_argument ("--input-csv",required =True ,help ="Input CSV path")
    parser .add_argument ("--gene-column",default ="gene_id",help ="Gene ID column name")
    parser .add_argument ("--output-dir",default ="",help ="Optional output directory")
    parser .add_argument ("--weights",default ="",help ="JSON string of scoring weights")
    args =parser .parse_args ()

    console =Console ()
    lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"

    default_cfg =BatchScoringConfig ()
    cfg =parse_weights (args .weights or None ,default_cfg )
    analyzer =ProteinBatchAnalyzer (config =cfg )

    console .print (
    Panel (
    "Batch Protein Analysis\nTwo-stage workflow: local scoring + literature boost + LLM synthesis"
    if lang =="en"
    else "Protein-coding gene batch analysis\nTwo-stage process: local scoring + literature enhancement + LLM comprehensive summary",
    border_style ="cyan",
    )
    )
    console .print (f"Input file: {args .input_csv }")
    console .print (f"Gene column: {args .gene_column }")

    result =analyzer .analyze_from_csv (args .input_csv ,args .gene_column )
    qc =result ["qc_report"]
    console .print ("\nQC Summary")
    console .print (
    f"input={qc ['input_count']}, unique={qc ['unique_count']}, "
    f"recognized={qc ['recognized_count']}, unrecognized={qc ['unrecognized_count']}, duplicates={qc ['duplicate_count']}"
    )
    if qc .get ("unrecognized_ids"):
        console .print (f"Unrecognized IDs (top 20): {', '.join (qc ['unrecognized_ids'][:20 ])}")

    rank_df =result ["rank_df"]
    scenario_df =result ["scenario_df"]
    assignment_matrix_df =result .get ("assignment_matrix_df")
    selected_scenario_tables =result .get ("selected_scenario_tables",{})

    console .print (f"\nRank table rows: {len (rank_df )}")
    if not rank_df .empty :
        show_cols =["rank","gene_id","primary_symbol","total_score","main_scenario","risk_tag"]
        console .print (rank_df [show_cols ].head (20 ).to_string (index =False ))

    console .print (f"\nScenario table rows: {len (scenario_df )}")
    if not scenario_df .empty :
        console .print (scenario_df .head (20 ).to_string (index =False ))

    if assignment_matrix_df is not None and not assignment_matrix_df .empty :
        console .print (f"\nScenario assignment matrix rows: {len (assignment_matrix_df )}")
        console .print (assignment_matrix_df .head (20 ).to_string (index =False ))

    for scenario ,df in selected_scenario_tables .items ():
        console .print (f"\nSelected scenario candidates: {scenario } ({len (df )})")
        if not df .empty :
            console .print (df .head (20 ).to_string (index =False ))

    console .print ("\nLLM Integrated Summary")
    console .print (Markdown (result ["llm_summary"]))

    bench =result .get ("benchmark",{})
    if bench .get ("available"):
        console .print ("\nBenchmark metrics")
        console .print (json .dumps (bench .get ("metrics",{}),ensure_ascii =False ,indent =2 ))
    else :
        console .print (f"\nBenchmark: {bench .get ('message','unavailable')}")

    rank_export =rank_df .drop (columns =["scenario_scores"],errors ="ignore")
    qc_payload ={
    "input_count":qc .get ("input_count",0 ),
    "unique_count":qc .get ("unique_count",0 ),
    "recognized_count":qc .get ("recognized_count",0 ),
    "unrecognized_count":qc .get ("unrecognized_count",0 ),
    "duplicate_count":qc .get ("duplicate_count",0 ),
    "unrecognized_ids":qc .get ("unrecognized_ids",[]),
    }
    payload ={
    "qc_report":qc_payload ,
    "rank_table":rank_export .to_dict ("records"),
    "scenario_table":scenario_df .to_dict ("records"),
    "assignment_matrix":assignment_matrix_df .to_dict ("records")if assignment_matrix_df is not None else [],
    "selected_scenario_tables":{k :v .to_dict ("records")for k ,v in selected_scenario_tables .items ()},
    "llm_summary":result ["llm_summary"],
    "cluster_reports":result .get ("cluster_reports",[]),
    "benchmark":bench ,
    }
    # Emit machine-readable payload without Rich formatting, otherwise wrapping
    # corrupts JSON extraction in the dashboard.
    print (BATCH_JSON_BEGIN )
    print (json .dumps (payload ,ensure_ascii =False ))
    print (BATCH_JSON_END )


if __name__ =="__main__":
    main ()
