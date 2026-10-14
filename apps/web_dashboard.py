import os 
import sys 
import subprocess 
import time 
import re 
import json 
import tempfile 
import math 
from pathlib import Path 

import pandas as pd 
import streamlit as st 
try :
    import plotly .graph_objects as go 
except Exception :
    go =None 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .utils .gene_annotator import RiceGeneAnnotator 
from src .utils .expression_pattern_service import ExpressionPatternService 
from src .utils .protein_feature_service import ProteinFeatureService 
from src .project_paths import ANNOTATION_DIR 


st .set_page_config (
page_title ="RIGOR",
page_icon ="RIGOR",
layout ="wide",
)


APP_MAP ={
"Protein Function Prediction":"apps/protein_function_pipeline.py",
"Protein Batch Analysis":"apps/protein_batch_pipeline.py",
"sRNA":"apps/mirna_pipeline.py",
"lncRNA":"apps/lncrna_pipeline.py",
"RNA Type Identification":"apps/sequence_identifier.py",
}

PAGE_OPTIONS =[
("protein","Protein Function Prediction"),
("batch","Protein Batch Analysis"),
("srna","sRNA"),
("lncrna","lncRNA"),
("identifier","RNA Type Identification"),
("help","Help"),
]

PAGE_LABEL_BY_KEY ={key :label for key ,label in PAGE_OPTIONS }

PROTEIN_NAV_SECTIONS =[
("Gene info","step-basic-info"),
("Homology","step-homology"),
("Co-expression","step-coexpression"),
("Epigenomics","step-epigenomic"),
("Expression","step-expression"),
("GWAS","step-gwas"),
("Radar","evidence-radar"),
]

BATCH_NAV_SECTIONS =[
("QC","batch-qc"),
("Scenario","batch-scenario-table"),
("Matrix","batch-assignment-matrix"),
("Candidates","batch-candidate-tables"),
("Rank","batch-rank-table"),
("Summary","batch-integrated-summary"),
]

SRNA_NAV_SECTIONS =[
("Targets","step-mirna-targets"),
("Function","step-mirna-function"),
("GWAS","step-mirna-gwas"),
("phasi targets","step-phasirna-targets"),
("phasi func","step-phasirna-function"),
("phasi GWAS","step-phasirna-gwas"),
]

LNCRNA_NAV_SECTIONS =[
("Basic","step-lncrna-basic"),
("Evidence","step-lncrna-supplement"),
("Co-expression","step-lncrna-coexpression"),
("GWAS","step-lncrna-gwas"),
]

IDENTIFIER_NAV_SECTIONS =[
]

HELP_NAV_SECTIONS =[
("Overview","help-overview"),
("Modules","help-modules"),
("Usage","help-usage"),
("Results","help-results"),
("FAQ","help-faq"),
]

NAV_SECTIONS_BY_PAGE ={
"protein":PROTEIN_NAV_SECTIONS ,
"batch":BATCH_NAV_SECTIONS ,
"srna":SRNA_NAV_SECTIONS ,
"lncrna":LNCRNA_NAV_SECTIONS ,
"identifier":IDENTIFIER_NAV_SECTIONS ,
"help":HELP_NAV_SECTIONS ,
}

ANSI_RE =re .compile (r"\x1b\[[0-9;?]*[A-Za-z]")
EMOJI_RE =re .compile (
"["
"\U0001F300-\U0001FAFF"
"\U00002700-\U000027BF"
"\U00002600-\U000026FF"
"\U0001F1E6-\U0001F1FF"
"]+",
flags =re .UNICODE ,
)
ERROR_HINT_RE =re .compile (r"(traceback|error|exception|failed|module not found|not found|report error)",re .IGNORECASE )
API_ERROR_HINT_RE =re .compile (
r"(401|403|unauthorized|invalid api key|incorrect api key|authentication|invalid[_\s-]*key|invalid[_\s-]*token|permission denied|insufficient_quota|api key)",
re .IGNORECASE ,
)
IMPORTANT_RUNTIME_RE =re .compile (
r"(\[Phase|\bPhase\b|Traceback|loading|scan|fetch|analysis|predict|done|failed|error)",
re .IGNORECASE ,
)
BATCH_JSON_BEGIN ="BATCH_RESULT_JSON_BEGIN"
BATCH_JSON_END ="BATCH_RESULT_JSON_END"
PROTEIN_JSON_BEGIN ="PROTEIN_RESULT_JSON_BEGIN"
PROTEIN_JSON_END ="PROTEIN_RESULT_JSON_END"
PROTEIN_STEP_BEGIN ="PROTEIN_STEP_JSON_BEGIN"
PROTEIN_STEP_END ="PROTEIN_STEP_JSON_END"
PROTEIN_PROGRESS_PREFIX ="PROTEIN_PROGRESS_JSON:"
NAV_PLACEHOLDER =None 
PLOT_RENDER_COUNTER =0 
BUILTIN_LLM_BASE_URL ="https://api.deepseek.com"
BUILTIN_LLM_API_KEY = os.getenv("DEEPSEEK_API_KEY", "")
BUILTIN_LLM_MODEL_NAME ="deepseek-v4-flash"


def _anchor_id (name ):
    safe =re .sub (r"[^a-zA-Z0-9_-]+","-",str (name )).strip ("-").lower ()
    return f"rigor-{safe }"


def render_nav_anchor (name ):
    st .markdown (f'<div id="{_anchor_id (name )}" class="rigor-anchor"></div>',unsafe_allow_html =True )


def _nav_links_html (page_key ):
    sections =NAV_SECTIONS_BY_PAGE .get (page_key ,[])
    return "".join (
    f'<a class="rigor-local-nav-link" href="#{_anchor_id (anchor )}">{label}</a>'
    for label ,anchor in sections 
    )


def _has_any_analysis_result ():
    result_keys =[
    "protein_result",
    "protein_batch_result",
    "srna_result",
    "lncrna_result",
    "identifier_result",
    ]
    return bool (st .session_state .get ("rigor_nav_auto_open"))or any (st .session_state .get (key )for key in result_keys )


def open_result_navigation ():
    st .session_state ["rigor_nav_auto_open"]=True 


def render_global_fixed_nav (placeholder =None ,force_open =None ):
    tab_html ="".join (
    f'<span class="rigor-nav-tab rigor-nav-tab-{key}">{label}</span>'
    for key ,label in PAGE_OPTIONS 
    )
    group_html ="".join (
    f"""
<div class="rigor-nav-group rigor-nav-group-{key}">
  <div class="rigor-local-nav-title">{PAGE_LABEL_BY_KEY .get (key ,"RIGOR")}</div>
  {_nav_links_html (key)}
</div>
"""
    for key ,_label in PAGE_OPTIONS 
    )
    has_result =_has_any_analysis_result ()if force_open is None else bool (force_open )
    checked_attr =" checked"if has_result else ""
    toggle_id ="rigor-nav-toggle-open"if has_result else "rigor-nav-toggle-closed"
    nav_html =f"""
<div class="rigor-nav-drawer">
  <input id="{toggle_id}" class="rigor-nav-toggle" type="checkbox"{checked_attr}>
  <label class="rigor-nav-toggle-label" for="{toggle_id}" title="Toggle navigation"></label>
  <div class="rigor-local-nav">
    <div class="rigor-nav-tabs">{tab_html}</div>
    {group_html}
  </div>
</div>
"""
    target =placeholder if placeholder is not None else st 
    target .markdown (nav_html ,unsafe_allow_html =True )


SCENARIO_OPTIONS =[
"development_architecture",
"hormone_signaling",
"membrane_vesicle_transport",
"stress_defense",
"metabolism_energy",
"epigenetic_transcriptional_regulation",
]

THEME_LABELS ={
"development_architecture":"Development",
"hormone_signaling":"Hormone/Signaling",
"membrane_vesicle_transport":"Membrane/Transport",
"stress_defense":"Stress/Defense",
"metabolism_energy":"Metabolism/Energy",
"epigenetic_transcriptional_regulation":"Epigenetic/Transcription",
}

THEME_COLORS ={
"development_architecture":"#2D6A4F",
"hormone_signaling":"#B56576",
"membrane_vesicle_transport":"#1D4ED8",
"stress_defense":"#C2410C",
"metabolism_energy":"#A16207",
"epigenetic_transcriptional_regulation":"#6D28D9",
}

ZH_EN_REPLACEMENTS =[
("Start sub-region","promoter region"),
("Start sub-area","promoter region"),
("Currently loading","Loading"),
("Currently loading","Loading"),
("Database has loaded","Database loaded"),
("Gene annotation database has loaded","Gene annotation database loaded"),
("Total of","Indexed"),
("rice genes collected","rice genes"),
("Standardized miRNA target database","standardized miRNA target database"),
("Macro co-expression community data","macro co-expression module data"),
("miRNA intelligent analyzer started","miRNA analysis started"),
("lncRNA intelligent analyzer started","lncRNA analysis started"),
("phasiRNA trait gene analysis","phasiRNA trait-gene analysis"),
("Target","Target"),
("Mode","Mode"),
("Database query","Database mode"),
("Sequence mode","Sequence mode"),
("Currently querying miRNA target gene","Querying miRNA target genes"),
("Successfully found","Found"),
("target genes","target gene(s)"),
("Currently fetching target gene literature","Fetching target-gene literature"),
("Currently performing miRNA function reverse inference","Running reverse functional inference for miRNA"),
("No literature found, skipping function inference","No literature found. Functional inference skipped"),
("Currently scanning GWAS signal of target gene","Scanning GWAS signals for target genes"),
("No significant GWAS signal found","No significant GWAS signals detected"),
("AnalysisCompleted","Analysis completed"),
("Currently analyzing lncRNA co-expression network","Analyzing lncRNA co-expression network"),
("Error","Error"),
("Not Found lncRNA","lncRNA not found"),
("Cleaned database not found","Clean database not found"),
("Successfully locked","Located"),
("candidate targets","candidate target(s)"),
("SearchFailed","Search failed"),
("Finally obtained","Retrieved"),
("literature blocks","literature block(s)"),
("Currently performing comprehensive reasoning by large model","Running integrated LLM inference"),
("Large model inference failed","LLM inference failed"),
("Macro Network positioning","Macro-network mapping"),
("Running status","Run Status"),
("Task has ended","Task completed"),
("Task failed","Task failed"),
("Exit code","exit code"),
("Detected suspected report error/anomalous information","Potential error/exception entries detected"),
("No obvious report error detected","No clear errors detected"),
("Running command","Command"),
("Running prompt","Runtime log"),
("Exception summary","Error summary"),
("Full report","Full report"),
("ErrorOutput","Error output"),
("AnalysisCompleted","Analysis completed"),
("Currently running, please wait","Running, please wait"),
("Task has been submitted, waiting for initial logs","Task submitted, waiting for logs"),
("Protein coding gene function prediction","Protein Function Prediction"),
("RNA type identification","RNA Type Identification"),
("Not Found","Not found"),
("Failed","failed"),
("Completed","completed"),
]

I18N ={
"zh":{
"app_title":"RiceTraitAgent Web Console",
"app_caption":"Unified entry point for protein genes, sRNAs, lncRNAs, and RNA type identification",
"language":"Language",
"runtime_config":"RunningConfiguration",
"llm_base_help":"Please input the complete interface address, e.g., https://api.deepseek.com",
"llm_key_help":"Please input the model service API key",
"runtime_notes":(
"Running instructions:\n"
"- This page reuses the main CLI entry point\n"
"- Automatically degrades when GWAS file is missing\n"
"- Enable LLM by filling in the sidebar configuration\n"
"- For first-time execution, run: pip install -r requirements.txt"
),
"tab_protein":"Protein Function Prediction",
"tab_batch":"Protein Batch Analysis",
"tab_srna":"sRNA Analysis",
"tab_lncrna":"lncRNA",
"tab_identifier":"RNA Type Identification",
"tab_help":"Help",
"runtime_hint":"Running log",
"runtime_pending":"Task has been submitted, waiting for initial logs...",
"run_status":"Running status",
"task_done":"Task has ended",
"task_failed":"Task failed, exit code {code}",
"error_detected":"Detected suspected report error/anomalous information {count} items",
"no_error":"No obvious report error detected",
"cmd":"Running command",
"runtime_log":"Running log",
"error_summary":"Exception summary",
"full_report":"Full report",
"error_output":"ErrorOutput",
},
"en":{
"app_title":"RiceTraitAgent Web Dashboard",
"app_caption":"Unified entry for protein, sRNA, lncRNA and RNA type identification",
"language":"Language",
"runtime_config":"Runtime Configuration",
"llm_base_help":"Enter full endpoint URL, e.g. https://api.deepseek.com",
"llm_key_help":"Enter API key for model service",
"runtime_notes":(
"Notes:\n"
"- All pages reuse CLI entry scripts\n"
"- Automatically degrades when GWAS files are missing\n"
"- Fill sidebar config to enable LLM\n"
"- First run: pip install -r requirements.txt"
),
"tab_protein":"Protein Function",
"tab_batch":"Batch Protein",
"tab_srna":"sRNA",
"tab_lncrna":"lncRNA",
"tab_identifier":"RNA Type ID",
"tab_help":"Help",
"running_wait":"Running, please wait...",
"runtime_hint":"Runtime Log",
"runtime_pending":"Task submitted. Waiting for logs...",
"run_status":"Run Status",
"task_done":"Task completed",
"task_failed":"Task failed, exit code {code}",
"error_detected":"Potential error/exception entries detected: {count}",
"no_error":"No clear errors detected",
"cmd":"Command",
"runtime_log":"Runtime Log",
"error_summary":"Error Summary",
"full_report":"Full Report",
"error_output":"Error Output",
},
}




def get_lang ():
    return "en"


def t (key ,**kwargs ):
    lang =get_lang ()
    text =I18N .get (lang ,I18N ["en"]).get (key ,key )
    return text .format (**kwargs )if kwargs else text 


def sanitize_output (text ):
    if not text :
        return ""
    text =text .replace ("\r","\n")
    text =ANSI_RE .sub ("",text )
    text =EMOJI_RE .sub ("",text )
    text ="".join (ch for ch in text if ch =="\n"or ch =="\t"or ord (ch )>=32 )
    text =text .replace ("Start sub-region","Start sub-area")
    return text 


def _is_structured_runtime_line (line ):
    if not line :
        return False 
    if line in (BATCH_JSON_BEGIN ,BATCH_JSON_END ,PROTEIN_JSON_BEGIN ,PROTEIN_JSON_END ,
    PROTEIN_STEP_BEGIN ,PROTEIN_STEP_END ):
        return True 
    if line .startswith (PROTEIN_PROGRESS_PREFIX ):
        return True 
    if line .startswith ('{"step"')or line .startswith ('{"step":'):
        return True 
    if line .startswith ('{"qc_report"')or line .startswith ('{"qc_report":'):
        return True 
    if line .startswith ('{"report_type"')or line .startswith ('{"report_type":'):
        return True 
    return False 


def strip_incremental_step_blocks (text ):
    if not text :
        return ""
    cleaned_lines =[]
    in_step_block =False 
    in_structured_result_block =False 
    for raw_line in text .splitlines ():
        line =raw_line .strip ()
        if line in (PROTEIN_JSON_BEGIN ,BATCH_JSON_BEGIN ):
            in_structured_result_block =True 
            cleaned_lines .append (raw_line )
            continue 
        if line in (PROTEIN_JSON_END ,BATCH_JSON_END ):
            in_structured_result_block =False 
            cleaned_lines .append (raw_line )
            continue 
        if in_structured_result_block :
            cleaned_lines .append (raw_line )
            continue 
        if line ==PROTEIN_STEP_BEGIN :
            in_step_block =True 
            continue 
        if line ==PROTEIN_STEP_END :
            in_step_block =False 
            continue 
        if in_step_block :
            continue 
        if _is_structured_runtime_line (line ):
            continue 
        cleaned_lines .append (raw_line )
    return "\n".join (cleaned_lines )


def run_command (script_rel_path ,args ,extra_env =None ):
    env =os .environ .copy ()
    env .update (extra_env or {})
    cmd =[sys .executable ,str (PROJECT_ROOT /script_rel_path ),*args ]
    result =subprocess .run (
    cmd ,
    cwd =PROJECT_ROOT ,
    env =env ,
    capture_output =True ,
    text =True ,
    )
    return cmd ,result 


LOGO_PATH =PROJECT_ROOT /"apps"/"rigor_logo.svg"
DARK_LOGO_PATH =PROJECT_ROOT /"apps"/"rigor_logo_dark.svg"


def render_logo_row ():
    with st .container (border =False ):
        left ,right =st .columns ([3.4 ,4.6 ])
        with left :
            use_dark =bool (st .session_state .get ("ui_dark_mode",False ))
            logo_to_use =DARK_LOGO_PATH if use_dark and DARK_LOGO_PATH .exists ()else LOGO_PATH
            if logo_to_use .exists ():
                st .image (str (logo_to_use ),width =460 )
        with right :
            st .empty ()


def render_app_footer ():
    st .markdown (
    """
<div class="rigor-app-footer">
  <div>© 2026, School of Life Sciences, Sun Yat-sen University, China. All Rights Reserved.</div>
  <div>Contact: <a href="mailto:jinwc@mail2.sysu.edu.cn">jinwc@mail2.sysu.edu.cn</a></div>
</div>
""",
    unsafe_allow_html =True ,
    )


def loading_animation_html (message ):
    return f"""
<div class="rigor-loader-wrap">
  <div class="rigor-loader-stage">
    <div class="rigor-loader-orb"></div>
  </div>
  <div class="rigor-loader-text">{message}</div>
</div>
"""


def radar_wait_message (lang ):
    return (
    "The radar chart will be generated after all evidence modules finish loading. Please wait..."
    if lang =="en"
    else "雷达图需要等待所有模块加载完毕之后生成，请耐心等候。"
    )


def llm_enabled_for_ui ():
    has_key =bool (get_effective_llm_config ()["api_key"])
    invalid =bool (st .session_state .get ("ui_llm_api_invalid",False ))
    return has_key and not invalid


def get_effective_llm_config ():
    return {
    "base_url":BUILTIN_LLM_BASE_URL ,
    "api_key":BUILTIN_LLM_API_KEY ,
    "model_name":BUILTIN_LLM_MODEL_NAME ,
    "source":"builtin",
    "ui_has_key":False ,
    "ui_has_model":False ,
    }


def render_api_required_hint (lang ):
    st .warning (
    "The built-in LLM service is temporarily unavailable. Please try again later."
    if lang =="en"
    else "The built-in LLM service is temporarily unavailable. Please try again later."
    )


def detect_api_error_message (stdout_text ,stderr_text ):
    text ="\n".join (part for part in [stdout_text ,stderr_text ]if part )
    if not text :
        return ""
    for line in text .splitlines ():
        line_clean =line .strip ()
        if _is_structured_runtime_line (line_clean ):
            continue 
        if line_clean and API_ERROR_HINT_RE .search (line_clean ):
            return line_clean 
    return ""


def validate_user_llm_config_once (llm_cfg ):
    if llm_cfg .get ("source")!="input"or not llm_cfg .get ("api_key"):
        return True ,""
    cache =st .session_state .setdefault ("llm_validation_cache",{})
    sig =(llm_cfg .get ("base_url",""),llm_cfg .get ("model_name",""),llm_cfg .get ("api_key",""))
    if sig in cache :
        return cache [sig ]
    try :
        from src .utils .llm_client import chat_completion_text ,make_openai_client 
        client =make_openai_client (
        api_key =llm_cfg .get ("api_key",""),
        base_url =llm_cfg .get ("base_url",""),
        timeout =20.0 ,
        )
        chat_completion_text (
        client ,
        llm_cfg .get ("model_name","deepseek-v4-flash"),
        [{"role":"user","content":"Reply with OK."}],
        temperature =0.0 ,
        max_tokens =64 ,
        )
        cache [sig ]=(True ,"")
        return True ,""
    except Exception as exc :
        msg =str (exc ).strip ()or "Authentication failed."
        cache [sig ]=(False ,msg )
        return False ,msg 


def translate_text_to_en (text ):
    out =text 
    for zh ,en in ZH_EN_REPLACEMENTS :
        out =out .replace (zh ,en )
    return out 


def normalize_display_text (text ,lang ="en"):
    cleaned =sanitize_output (text )
    if lang =="en":
        cleaned =translate_text_to_en (cleaned )
    return cleaned 


def extract_important_runtime_lines (text ,limit =10 ):
    lines =[line .strip ()for line in sanitize_output (text ).splitlines ()if line .strip ()]
    noise_patterns =[
    r"^={10,}$",
    r"^[-_]{6,}$",
    r"^Rich .*",
    ]
    filtered =[]
    for line in lines :
        if _is_structured_runtime_line (line ):
            continue 
        if any (re .search (pat ,line )for pat in noise_patterns ):
            continue 
        filtered .append (line )
    lines =filtered 
    important =[line for line in lines if IMPORTANT_RUNTIME_RE .search (line )]
    if not important :
        important =lines [-limit :]
    out ="\n".join (important [-limit :])
    return normalize_display_text (out ,get_lang ())


def _parse_step_json (buffer_lines ):
    """Attempting to parse STEP JSON from buffered lines"""
    try :
        text ="".join (buffer_lines ).strip ()
        return json .loads (text )
    except Exception :
        return None 


def _plot_key (name ):
    return f"{name }_{PLOT_RENDER_COUNTER }"


def _parse_progress_line (line ):
    if not line .startswith (PROTEIN_PROGRESS_PREFIX ):
        return None 
    try :
        return json .loads (line [len (PROTEIN_PROGRESS_PREFIX ):].strip ())
    except Exception :
        return None 


def render_coexpression_literature_progress (payload ,placeholder ):
    if not payload or payload .get ("scope")!="coexpression_literature":
        return 
    done =int (payload .get ("done",0 )or 0 )
    total =max (1 ,int (payload .get ("total",1 )or 1 ))
    cumulative =int (payload .get ("cumulative_literature",0 )or 0 )
    maximum =max (1 ,int (payload .get ("max_literature",total *10 )or total *10 ))
    gene_id =payload .get ("gene_id","")
    symbol =payload .get ("symbol","")
    gene_label =""
    if gene_id or symbol :
        gene_label =f" Current gene: {gene_id } {f'({symbol })'if symbol else ''}"
    fraction =min (1.0 ,max (0.0 ,done /total ))
    with placeholder .container ():
        with st .container (border =True ):
            st .caption ("Co-expression literature retrieval")
            if payload .get ("rate_limited"):
                st .warning (
                "PubMed is temporarily rate-limiting requests. Please wait patiently; RIGOR will continue once the remote service responds. "
                "If this takes too long, you can stop this run and rerun later to preserve full literature coverage."
                )
            st .progress (
            fraction ,
            text =f"{done }/{total } companion genes scanned · {cumulative }/{maximum } literature records retrieved.{gene_label }",
            )


def _render_step_card (step_data ,lang ):
    """Render complete AnalysisResult card based on step type"""
    step =step_data .get ("step","")
    title =step_data .get ("title","")
    data =step_data .get ("data",{})
    llm_enabled =llm_enabled_for_ui ()
    render_nav_anchor (f"step-{step .replace ('_','-')}")

    if step =="basic_info":
        with st .container (border =True ):
            st .subheader (title )
            col1 ,col2 ,col3 =st .columns (3 )
            col1 .metric ("Target"if lang =="en"else "Target",data .get ("target","N/A"))
            col2 .metric ("Symbol",data .get ("symbol","N/A"))
            if data .get ("uniprot"):
                col3 .metric ("UniProt",data ["uniprot"])
            else :
                col3 .metric ("Type"if lang =="en"else "Type","protein_gene")
            if data .get ("description"):
                st .caption ("Description"if lang =="en"else "Description")
                st .write (data ["description"])
            if data .get ("location"):
                st .caption ("Location"if lang =="en"else "Coordinates")
                st .write (data ["location"])
                #Protein features
            plen =data .get ("protein_length")
            subcell =data .get ("subcellular_locations",[])
            interpro =data .get ("interpro_domains",[])
            pfam =data .get ("pfam_domains",[])
            if plen or subcell :
                st .divider ()
                st .caption ("Protein Features"if lang =="en"else "Protein features")
                pcols =st .columns (3 )
                pcols [0 ].metric ("Length"if lang =="en"else "Protein length",
                f"{plen } aa"if plen else "N/A")
                pcols [1 ].metric ("InterPro",len (interpro ))
                pcols [2 ].metric ("Pfam",len (pfam ))
                #Display full Information for subcellular localization in a separate line
                st .caption ("Subcellular Localization"if lang =="en"else "Subcellular localization")
                st .write (", ".join (subcell )if subcell else "N/A")
                #InterPro Domain table
            if interpro :
                st .caption ("InterPro Domains"if lang =="en"else "InterPro Domain")
                st .dataframe (pd .DataFrame (interpro ),use_container_width =True ,hide_index =True )
                #Pfam Core Domain table
            if pfam :
                st .caption ("Pfam Core Domains"if lang =="en"else "Pfam Core Domain")
                st .dataframe (pd .DataFrame (pfam ),use_container_width =True ,hide_index =True )
                #Pfam Domain plot
            plot_b64 =data .get ("plot_png_b64","")
            if plot_b64 :
                import base64 as _b64 
                st .caption ("Pfam Domain Plot"if lang =="en"else "Pfam Domain plot")
                st .image (_b64 .b64decode (plot_b64 ),use_container_width =True )

    elif step =="homology":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .info ("No homology literature found."if lang =="en"else "No homologous literature found.")
                #Literature list
            refs =data .get ("references",[])
            if refs :
                st .divider ()
                st .caption (f"References ({len (refs )})"if lang =="en"else f"referenceLiterature（{len (refs )} piece）")
                for ref in refs :
                    pmid =ref .get ("pmid","")
                    ref_title =ref .get ("title","")
                    journal =ref .get ("journal","")
                    st .markdown (f"- **{ref_title }** — {journal } ([PMID:{pmid }](https://pubmed.ncbi.nlm.nih.gov/{pmid }/))")

    elif step =="coexpression":
        with st .container (border =True ):
            st .subheader (title )
            #GO Enrichment bar chart (Top 10 displayed, based on Top 100 co-expressed genes)
            go_terms =data .get ("go_terms",[])
            net_count =data .get ("network_gene_count",0 )
            if go_terms and go is not None :
                st .caption (
                f"GO Enrichment Top 10 (from {net_count } co-expressed genes)"
                if lang =="en"
                else f"GO Enrichment Top 10（based {net_count } commonExpressionGenestatistics）"
                )
                terms =[item ["term"]for item in go_terms ]
                counts =[item ["count"]for item in go_terms ]
                fig_go =go .Figure (go .Bar (
                x =counts [::-1 ],
                y =terms [::-1 ],
                orientation ="h",
                marker_color ="#0F766E",
                ))
                fig_go .update_layout (
                margin =dict (l =20 ,r =20 ,t =10 ,b =20 ),
                height =max (200 ,len (go_terms )*32 ),
                xaxis_title ="Count",
                yaxis =dict (automargin =True ),
                )
                st .plotly_chart (fig_go ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("coexpression_go"))
            elif go_terms :
            #Fallback to text list when Plotly is not installed
                st .caption (
                f"GO Enrichment Top 10 (from {net_count } co-expressed genes)"
                if lang =="en"
                else f"GO Enrichment Top 10（based {net_count } commonExpressionGenestatistics）"
                )
                for item in go_terms :
                    st .write (f"- {item ['term']} ({item ['count']})")
            if data .get ("pending"):
                st .divider ()
                pending_message =data .get ("message")or (
                "Co-expression companion-gene literature retrieval can take longer than other modules. Please wait patiently."
                if lang =="en"
                else "共表达伙伴文献检索通常耗时较长，请耐心等候。"
                )
                st .markdown (loading_animation_html (pending_message ),unsafe_allow_html =True )
                progress =data .get ("progress")or {}
                if progress :
                    total =max (1 ,int (progress .get ("total",1 )or 1 ))
                    done =int (progress .get ("done",0 )or 0 )
                    cumulative =int (progress .get ("cumulative_literature",0 )or 0 )
                    maximum =max (1 ,int (progress .get ("max_literature",total *10 )or total *10 ))
                    gene_id =progress .get ("gene_id","")
                    symbol =progress .get ("symbol","")
                    gene_label =""
                    if gene_id or symbol :
                        gene_label =(
                        f" Current gene: {gene_id } {f'({symbol })'if symbol else ''}"
                        if lang =="en"
                        else f" 当前基因：{gene_id } {f'({symbol })'if symbol else ''}"
                        )
                    fraction =min (1.0 ,max (0.0 ,done /total ))
                    if progress .get ("rate_limited"):
                        st .warning (
                        "PubMed is temporarily rate-limiting requests. RIGOR will continue once the remote service responds."
                        if lang =="en"
                        else "PubMed 当前可能存在访问限速。RIGOR 会在远程服务响应后继续执行，请耐心等候。"
                        )
                    st .progress (
                    fraction ,
                    text =(
                    f"{done }/{total } companion genes scanned · {cumulative }/{maximum } literature records retrieved.{gene_label }"
                    if lang =="en"
                    else f"已检索 {done }/{total } 个共表达伙伴基因 · 已获取 {cumulative }/{maximum } 条文献记录。{gene_label }"
                    ),
                    )
                else :
                    st .progress (
                    0 ,
                    text =(
                    "Waiting for companion-gene literature retrieval..."
                    if lang =="en"
                    else "正在等待共表达伙伴文献检索..."
                    ),
                    )
                return 
            #LLM Report
            report =data .get ("report","")
            if not llm_enabled :
                st .divider ()
                render_api_required_hint (lang )
            elif report :
                st .divider ()
                st .markdown (normalize_display_text (report ,lang ))
                #Literature list
            refs =data .get ("references",[])
            if refs :
                st .divider ()
                st .caption (f"References ({len (refs )})"if lang =="en"else f"referenceLiterature（{len (refs )} piece）")
                for ref in refs :
                    pmid =ref .get ("pmid","")
                    ref_title =ref .get ("title","")
                    journal =ref .get ("journal","")
                    st .markdown (f"- **{ref_title }** — {journal } ([PMID:{pmid }](https://pubmed.ncbi.nlm.nih.gov/{pmid }/))")

    elif step =="epigenomic":
        with st .container (border =True ):
            st .subheader (title )
            #Start site epigenomic plot
            heatmap_data =data .get ("heatmap_data",[])
            if heatmap_data :
                render_promoter_chart (heatmap_data ,lang )
                #LLM Report
            report =data .get ("report","")
            if not llm_enabled :
                st .divider ()
                render_api_required_hint (lang )
            elif report :
                st .divider ()
                st .markdown (normalize_display_text (report ,lang ))

    elif step =="expression":
        with st .container (border =True ):
            st .subheader (title )
            summary =data .get ("summary","")
            if summary :
                st .markdown (normalize_display_text (summary ,lang ))
                # WarningInformation
            for warning in data .get ("warnings",[]):
                st .warning (normalize_display_text (warning ,lang ))
                #Transcript expression plot
            import base64 as _b64 
            for tr in data .get ("transcript_results",[]):
                with st .container (border =True ):
                    st .markdown (
                    f"**Transcript {tr ['transcript_id']}**"
                    if lang =="en"
                    else f"**Transcript {tr ['transcript_id']}**"
                    )
                    st .caption (
                    f"Expression model: {tr .get ('expression_model_id','N/A')}"
                    if lang =="en"
                    else f"Expressionmatrix model: {tr .get ('expression_model_id','N/A')}"
                    )
                    for w in tr .get ("warnings",[]):
                        st .warning (normalize_display_text (w ,lang ))
                    col1 ,col2 =st .columns (2 )
                    with col1 :
                        st .caption ("Organ/time expression"if lang =="en"else "Organism/temporal expression pattern")
                        if tr .get ("developmental_plot_b64"):
                            st .image (_b64 .b64decode (tr ["developmental_plot_b64"]),use_container_width =True )
                        else :
                            st .info ("No organ/time expression plot."if lang =="en"else "No organism/temporal expression plot generated currently.")
                        if tr .get ("developmental_summary"):
                            st .dataframe (pd .DataFrame (tr ["developmental_summary"]),use_container_width =True ,hide_index =True )
                    with col2 :
                        st .caption ("Treatment expression"if lang =="en"else "Treatment expression pattern")
                        if tr .get ("treatment_plot_b64"):
                            st .image (_b64 .b64decode (tr ["treatment_plot_b64"]),use_container_width =True )
                        else :
                            st .info ("No treatment expression plot."if lang =="en"else "No treatment expression plot generated currently.")
                        if tr .get ("treatment_summary"):
                            st .dataframe (pd .DataFrame (tr ["treatment_summary"]),use_container_width =True ,hide_index =True )
            if not data .get ("transcript_results")and not summary :
                st .info ("No expression data."if lang =="en"else "NoExpressionData。")

    elif step =="gwas":
        with st .container (border =True ):
            st .subheader (title )
            #Trait effect plot
            trait_data =data .get ("trait_effect_data",[])
            suggestive_trait_data =data .get ("suggestive_trait_effect_data",[])
            if trait_data or suggestive_trait_data :
                render_trait_effect_plot ({"trait_effect_plot":trait_data ,"suggestive_trait_effect_plot":suggestive_trait_data },lang )
            elif data .get ("suggestive_trait_effect_data"):
                st .info ("No Strong/Moderate GWAS evidence for the target gene. Weak trait-normalized signals are listed below as suggestive evidence only."if lang =="en"else "目标基因未检出 Strong/Moderate 级别的核心 GWAS 证据。Weak 信号仅作为提示列在下方。")
            if suggestive_trait_data :
                with st .expander ("Suggestive GWAS signals (Weak)"if lang =="en"else "提示性 GWAS 信号（Weak）",expanded =False ):
                    st .caption ("Weak signals are within the trait-specific top 5% only and should not be used as core predictions."if lang =="en"else "Weak 信号只进入该性状全基因组背景的 top 5%，不应作为核心预测结论。")
                    st .dataframe (pd .DataFrame (suggestive_trait_data ),use_container_width =True ,hide_index =True )
                #GWAS LLM Report
            sections =data .get ("report_sections",[])
            if not llm_enabled :
                render_api_required_hint (lang )
            else :
                for section in sections :
                    sec_title =section .get ("title","")
                    body =normalize_display_text (section .get ("body","").strip (),lang )
                    if sec_title :
                        st .markdown (f"### {sec_title }")
                    if body :
                        st .markdown (body )

    elif step =="lncrna_basic":
        with st .container (border =True ):
            st .subheader (title )
            lncrna_id =data .get ("lncrna_id","")
            chrom =data .get ("chromosome","")
            start =data .get ("start",0 )
            end =data .get ("end",0 )
            strand =data .get ("strand",".")
            length =data .get ("length",0 )
            coord_hit =data .get ("coord_hit")

            #Basic positional information
            col1 ,col2 ,col3 =st .columns (3 )
            col1 .metric ("lncRNA ID",lncrna_id )
            col2 .metric ("Location"if lang =="en"else "Position",f"{chrom }:{start }-{end }"if chrom else "N/A")
            col3 .metric ("Length (bp)"if lang =="en"else "Length (bp)",f"{length :,}"if length else "N/A")

            col4 ,col5 =st .columns (2 )
            col4 .metric ("Strand"if lang =="en"else "Strand Direction",strand )
            col5 .metric (
            "Cis / Trans",
            f"Cis {data .get ('cis_count',0 )} / Trans {data .get ('trans_count',0 )}"
            )

            #Coordinate Search Hit Information
            if coord_hit :
                st .caption ("Coordinate Match"if lang =="en"else "Coordinate Match")
                st .info (
                f"Query → **{coord_hit ['id']}**  similarity {coord_hit ['similarity']*100 :.1f}%"
                if lang =="en"
                else f"query hit → **{coord_hit ['id']}**  similarity {coord_hit ['similarity']*100 :.1f}%"
                )

                #Shared Expression Gene Table
            coexp =data .get ("coexp_genes",[])
            if coexp :
                st .divider ()
                st .caption (
                f"Top Co-expressed Protein Genes ({len (coexp )})"
                if lang =="en"
                else f"totalExpressionproteinCodingGene（{len (coexp )} each）"
                )
                df =pd .DataFrame (coexp )[["rap_id","symbol","correlation","regulation","description"]]
                df .columns =(
                ["RAP ID","Symbol","Pearson r","Regulation","Description"]
                if lang =="en"
                else ["RAP ID","Symbol","Pearson r","Regulation Type","Description"]
                )
                st .dataframe (df ,use_container_width =True ,hide_index =True )

    elif step =="lncrna_coexpression":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .info ("No functional inference available."if lang =="en"else "No Function Inference Result.")
            refs =data .get ("references",[])
            if refs :
                st .divider ()
                st .caption (f"References ({len (refs )})"if lang =="en"else f"referenceLiterature（{len (refs )} piece）")
                for ref in refs :
                    pmid =ref .get ("pmid","")
                    ref_title =ref .get ("title","")
                    journal =ref .get ("journal","")
                    st .markdown (f"- **{ref_title }** — {journal } ([PMID:{pmid }](https://pubmed.ncbi.nlm.nih.gov/{pmid }/))")

    elif step =="lncrna_gwas":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                has_signal =data .get ("has_signal",False )
                st .info (
                ("No significant GWAS signals detected in co-expressed gene network."
                if lang =="en"
                else "No significant GWAS Signal detected in the shared Expression Gene Network.")
                if not has_signal 
                else ("GWAS signals detected but summary unavailable."
                if lang =="en"
                else "GWAS Signal detected, but summary generation failed.")
                )

    elif step =="mirna_targets":
        with st .container (border =True ):
            st .subheader (title )
            mirna_id =data .get ("mirna_id","")
            mode =data .get ("mode","")
            target_count =data .get ("target_count",0 )
            annotated_count =data .get ("annotated_count",0 )
            cols =st .columns (3 )
            cols [0 ].metric ("miRNA"if lang =="en"else "miRNA ID",mirna_id )
            cols [1 ].metric ("Mode"if lang =="en"else "Mode",mode )
            cols [2 ].metric ("Targets"if lang =="en"else "Number of Target Genes",f"{annotated_count }/{target_count }")
            targets =data .get ("targets",[])
            if targets :
                st .caption ("Target Genes"if lang =="en"else "Target Gene List")
                st .dataframe (
                pd .DataFrame (targets ).rename (columns ={
                "rap_id":"RAP ID",
                "symbol":"Symbol",
                "evidence":"Evidence",
                "score":"Score",
                "seed_matches":"Seed Matches",
                "has_pmid":"Has Papers",
                }),
                use_container_width =True ,
                hide_index =True ,
                )

    elif step =="mirna_function":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .info ("No literature found for target genes. Functional inference skipped."
                if lang =="en"
                else "No literature related to target genes found; skipping Function Inference.")

    elif step =="mirna_gwas":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            gwas_signals =data .get ("gwas_signals",[])
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .info ("No significant GWAS signals detected for target genes."
                if lang =="en"
                else "No significant GWAS Signal detected in target genes.")
            if gwas_signals :
                st .divider ()
                st .caption (f"GWAS Signals ({len (gwas_signals )})"if lang =="en"else f"GWAS Signal（{len (gwas_signals )} line）")
                st .dataframe (pd .DataFrame (gwas_signals ),use_container_width =True ,hide_index =True )

    elif step =="phasirna_targets":
        with st .container (border =True ):
            st .subheader (title )
            phasirna_id =data .get ("phasirna_id","")
            sequence =data .get ("sequence","")
            evidence_mode =data .get ("evidence_mode","")
            target_count =data .get ("target_count",0 )
            cols =st .columns (3 )
            cols [0 ].metric ("phasiRNA ID",phasirna_id )
            cols [1 ].metric ("Sequence"if lang =="en"else "Sequence",f"{sequence } ({len (sequence )}nt)"if sequence else "N/A")
            cols [2 ].metric ("Targets"if lang =="en"else "Number of Target Genes",target_count )
            if evidence_mode :
                mode_label ={
                "expression_anticorrelation":("Expression Anti-correlation"if lang =="en"else "Negative Correlation with Expression"),
                "heuristic_complementarity":("Heuristic Complementarity"if lang =="en"else "Heuristic Complementarity"),
                }.get (evidence_mode ,evidence_mode )
                st .caption (f"{'Evidence Mode'if lang =='en'else 'Evidence Pattern'}: {mode_label }")
            targets =data .get ("targets",[])
            if targets :
                st .dataframe (
                pd .DataFrame (targets ).rename (columns ={
                "rap_id":"RAP ID","symbol":"Symbol",
                "anti_corr":"Anti-Corr","score":"Score",
                "seed_matches":"Seed Matches","has_pmid":"Has Papers",
                }),
                use_container_width =True ,hide_index =True ,
                )

    elif step =="phasirna_function":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .info ("No literature found for target genes. Functional inference skipped."
                if lang =="en"
                else "No literature related to target genes found; skipping Function Inference.")

    elif step =="phasirna_gwas":
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            gwas_signals =data .get ("gwas_signals",[])
            if not llm_enabled :
                render_api_required_hint (lang )
            elif report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .info ("No significant GWAS signals detected for target genes."
                if lang =="en"
                else "No significant GWAS Signal detected in target genes.")
            if gwas_signals :
                st .divider ()
                st .caption (f"GWAS Signals ({len (gwas_signals )})"if lang =="en"else f"GWAS Signal（{len (gwas_signals )} line）")
                st .dataframe (pd .DataFrame (gwas_signals ),use_container_width =True ,hide_index =True )

    else :
        with st .container (border =True ):
            st .subheader (title )
            report =data .get ("report","")
            if report :
                st .markdown (normalize_display_text (report ,lang ))
            else :
                st .json (data )


                #─── Step Name Mapping (for progress indication, key matches pipeline _emit_step) ───
_STEP_NAMES_ZH ={
"basic_info":"Gene information and protein features",
"homology":"Homologous literature function inference",
"coexpression":"Co-expression network function inference",
"epigenomic":"Start sub-epigenetic regulation",
"expression":"Multi-transcript expression pattern",
"gwas":"GWAS trait and breeding evaluation",
"lncrna_basic":"lncRNA Basic Information",
"lncrna_supplement":"Supplement Additional Evidence",
"lncrna_coexpression":"Co-expression Function Inference",
"lncrna_gwas":"GWAS Breeding Effect Evaluation",
"mirna_targets":"Target Gene Query/Prediction",
"mirna_function":"Function reverse inference",
"mirna_gwas":"GWAS Breeding Effect Evaluation",
"phasirna_targets":"Candidate target gene identification",
"phasirna_function":"Function reverse inference",
"phasirna_gwas":"GWAS breeding value evaluation",
}
_STEP_NAMES_EN ={
"basic_info":"Gene Info & Protein Features",
"homology":"Homology Literature Inference",
"coexpression":"Co-expression Network Inference",
"epigenomic":"Promoter Epigenomics",
"expression":"Expression Patterns",
"gwas":"GWAS Trait & Breeding Assessment",
"lncrna_basic":"lncRNA Basic Info",
"lncrna_supplement":"Supplement Evidence",
"lncrna_coexpression":"Co-expression Inference",
"lncrna_gwas":"GWAS Breeding Assessment",
"mirna_targets":"Target Gene Query/Prediction",
"mirna_function":"Functional Reverse Inference",
"mirna_gwas":"GWAS Breeding Assessment",
"phasirna_targets":"Target Gene Identification",
"phasirna_function":"Functional Reverse Inference",
"phasirna_gwas":"GWAS Breeding Assessment",
}

STEP_RENDER_ORDER =[
"basic_info",
"homology",
"coexpression",
"epigenomic",
"expression",
"gwas",
"lncrna_basic",
"lncrna_supplement",
"lncrna_coexpression",
"lncrna_gwas",
"mirna_targets",
"mirna_function",
"mirna_gwas",
"phasirna_targets",
"phasirna_function",
"phasirna_gwas",
]


def _render_ordered_step_cards (step_payloads ,lang ):
    global PLOT_RENDER_COUNTER 
    PLOT_RENDER_COUNTER +=1 
    ordered_steps =[step for step in STEP_RENDER_ORDER if step in step_payloads ]
    remaining_steps =[step for step in step_payloads if step not in STEP_RENDER_ORDER ]
    for step in ordered_steps +remaining_steps :
        _render_step_card (step_payloads [step ],lang )


def execute_and_store (result_key ,script_rel_path ,args ,extra_env =None ):
    st .session_state ["rigor_nav_auto_open"]=True 
    if NAV_PLACEHOLDER is not None :
        render_global_fixed_nav (NAV_PLACEHOLDER ,force_open =True )

    # Render immediate feedback before API validation or subprocess startup.
    status_placeholder =st .empty ()
    step_container =st .empty ()#Step Card Rendering Area
    runtime_box =st .empty ()
    coexpression_literature_progress_box =st .empty ()
    combined_output =[]
    step_payloads ={}

    runtime_box .markdown (
    loading_animation_html ("Preparing runtime..."),
    unsafe_allow_html =True ,
    )

    env =os .environ .copy ()
    env .update (extra_env or {})
    llm_cfg =get_effective_llm_config ()
    llm_ok ,llm_err_msg =validate_user_llm_config_once (llm_cfg )
    if not llm_ok :
        st .session_state ["ui_llm_api_invalid"]=True 
        st .session_state ["ui_llm_api_invalid_msg"]=llm_err_msg 
        st .error (f"Built-in LLM service is unavailable or has insufficient permission.\\n\\nDetected: {llm_err_msg }")
        env ["LLM_API_KEY"]=""
        env ["BATCH_ENABLE_LLM_SUMMARY"]="0"
        env ["RIGOR_LLM_ENABLED"]="0"
    else :
        st .session_state ["ui_llm_api_invalid"]=False 
        st .session_state ["ui_llm_api_invalid_msg"]=""
    cmd =[sys .executable ,str (PROJECT_ROOT /script_rel_path ),*args ]
    lang =get_lang ()
    is_protein ="protein_function_pipeline"in script_rel_path 
    is_lncrna ="lncrna_pipeline"in script_rel_path 
    is_mirna ="mirna_pipeline"in script_rel_path 
    is_step_pipeline =is_protein or is_lncrna or is_mirna 

    runtime_box .markdown (
    loading_animation_html ("RIGOR is reasoning across evidence modules..."),
    unsafe_allow_html =True ,
    )

    process =subprocess .Popen (
    cmd ,
    cwd =PROJECT_ROOT ,
    env =env ,
    stdout =subprocess .PIPE ,
    stderr =subprocess .STDOUT ,
    text =True ,
    bufsize =1 ,
    )

    #── Real-time Read Loop ──
    step_buffer =[]
    in_step_block =False 

    if process .stdout is not None :
        for line in iter (process .stdout .readline ,""):
            if not line :
                break 
            combined_output .append (line )
            stripped =line .strip ()

            progress_payload =_parse_progress_line (stripped )
            if progress_payload :
                if is_protein and progress_payload .get ("scope")=="coexpression_literature":
                    runtime_box .markdown (
                    loading_animation_html (radar_wait_message (lang )),
                    unsafe_allow_html =True ,
                    )
                coexpression_literature_progress_box .empty ()
                existing =step_payloads .get ("coexpression",{
                "step":"coexpression",
                "title":"Co-expression Network Inference",
                "status":"running",
                "data":{},
                })
                existing_data =existing .setdefault ("data",{})
                existing_data ["pending"]=True 
                existing_data ["progress"]=progress_payload 
                existing_data ["message"]=(
                "Co-expression companion-gene literature retrieval is running. This step can take longer than other modules; please wait patiently."
                if lang =="en"
                else "共表达伙伴文献检索正在运行，这一步耗时通常较长，请耐心等候。"
                )
                step_payloads ["coexpression"]=existing 
                with step_container .container ():
                    _render_ordered_step_cards (step_payloads ,lang )
                continue 

            #Detecting STEP JSON Marker
            if is_step_pipeline and stripped ==PROTEIN_STEP_BEGIN :
                in_step_block =True 
                step_buffer =[]
                continue 
            elif is_step_pipeline and stripped ==PROTEIN_STEP_END and in_step_block :
                in_step_block =False 
                step_data =_parse_step_json (step_buffer )
                if step_data :
                    if step_data .get ("step")=="coexpression":
                        coexpression_literature_progress_box .empty ()
                    step_payloads [step_data .get ("step","")]=step_data 
                    with step_container .container ():
                        _render_ordered_step_cards (step_payloads ,lang )
                step_buffer =[]
                continue 
            elif in_step_block :
                step_buffer .append (line )
                continue 
            time .sleep (0.01 )
        process .stdout .close ()

    returncode =process .wait ()
    status_placeholder .empty ()
    stdout_text =strip_incremental_step_blocks (sanitize_output ("".join (combined_output )))
    runtime_box .empty ()

    st .session_state [result_key ]={
    "cmd":cmd ,
    "returncode":returncode ,
    "stdout":stdout_text ,
    "stderr":"",
    }


def extract_runtime_summary (stdout ,stderr ,returncode ):
    combined ="\n".join (part for part in [stdout ,stderr ]if part ).strip ()
    lines =[line for line in combined .splitlines ()if line .strip ()]
    error_lines =[line for line in lines if ERROR_HINT_RE .search (line )]
    recent_lines =lines [-12 :]
    return {
    "error_count":len (error_lines ),
    "recent_lines":"\n".join (recent_lines ),
    "error_lines":"\n".join (error_lines [:20 ]),
    "has_errors":returncode !=0 or bool (error_lines ),
    }


def render_runtime_panel (cmd ,result ):
    return 


def render_full_report (result ):
    lang =get_lang ()
    if result .stderr .strip ():
        with st .container (border =True ):
            st .subheader (t ("error_output"))
            st .code (normalize_display_text (result .stderr ,lang ))


def render_run_result (cmd ,result ):
    render_runtime_panel (cmd ,result )
    render_full_report (result )


def render_stored_result (result_key ):
    payload =st .session_state .get (result_key )
    if not payload :
        return 
    class ResultProxy :
        def __init__ (self ,data ):
            self .returncode =data ["returncode"]
            self .stdout =data ["stdout"]
            self .stderr =data ["stderr"]
    render_run_result (payload ["cmd"],ResultProxy (payload ))


def extract_batch_payload (stdout_text ):
    text =stdout_text or ""
    if BATCH_JSON_BEGIN not in text or BATCH_JSON_END not in text :
        return None 
    try :
        segment =text .split (BATCH_JSON_BEGIN ,1 )[1 ].split (BATCH_JSON_END ,1 )[0 ].strip ()
        return json .loads (segment )
    except Exception :
        return None 


def extract_protein_payload (stdout_text ):
    text =stdout_text or ""
    if PROTEIN_JSON_BEGIN not in text or PROTEIN_JSON_END not in text :
        return None 
    try :
        segment =text .split (PROTEIN_JSON_BEGIN ,1 )[1 ].split (PROTEIN_JSON_END ,1 )[0 ].strip ()
        return json .loads (segment )
    except Exception :
        return None 


def render_evidence_radar (viz ,lang ):
    if go is None :
        st .info ("Plotly is not installed. Radar plot is unavailable."if lang =="en"else "Plotly is not installed; cannot display evidence radar chart.")
        return 
    radar =(viz or {}).get ("evidence_radar",{})
    scores =radar .get ("scores",{})
    if not scores :
        st .info ("No evidence radar data."if lang =="en"else "No evidence radar chart data.")
        return 
    categories =["literature","domain_go","coexpression","epigenomics","gwas"]
    labels ={
    "literature":"Literature",
    "domain_go":"Domain/GO",
    "coexpression":"Co-expression",
    "epigenomics":"Epigenomics",
    "gwas":"GWAS",
    }
    values =[float (scores .get (item ,0.0 ))for item in categories ]
    categories_closed =[labels [item ]for item in categories ]+[labels [categories [0 ]]]
    values_closed =values +[values [0 ]]
    fig =go .Figure ()
    fig .add_trace (
    go .Scatterpolar (
    r =values_closed ,
    theta =categories_closed ,
    fill ="toself",
    line =dict (color ="#0F766E",width =3 ),
    fillcolor ="rgba(15,118,110,0.20)",
    name ="Evidence",
    )
    )
    fig .update_layout (
    margin =dict (l =20 ,r =20 ,t =20 ,b =20 ),
    polar =dict (radialaxis =dict (range =[0 ,100 ],tickfont =dict (size =10 ))),
    showlegend =False ,
    height =360 ,
    )
    st .plotly_chart (fig ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("evidence_radar"))


def render_promoter_heatmap (viz ,lang ):
    if go is None :
        st .info ("Plotly is not installed. Heatmap is unavailable."if lang =="en"else "Plotly is not installed, cannot display Start sub heatmap.")
        return 
    rows =(viz or {}).get ("promoter_heatmap",[])
    if not rows :
        st .info ("No promoter heatmap data."if lang =="en"else "No Start sub heatmap data.")
        return 
    df =pd .DataFrame (rows )
    value_cols =[c for c in ["open_count","active_count","repressive_count","poised_count","context_count","score"]if c in df .columns ]
    if not value_cols :
        st .info ("No promoter heatmap columns available."if lang =="en"else "No Start sub heatmap column data.")
        return 
    z =df [value_cols ].to_numpy ()
    col_labels ={"open_count":"Open","active_count":"Active","repressive_count":"Repressive",
    "poised_count":"Poised","context_count":"Context","score":"Score"}
    x_labels =[col_labels .get (c ,c )for c in value_cols ]
    customdata =df [["state"]].to_numpy ()if "state"in df .columns else None 
    fig =go .Figure (
    data =go .Heatmap (
    z =z ,
    x =x_labels ,
    y =df ["tissue"],
    colorscale ="YlGnBu",
    colorbar =dict (title ="Value"),
    customdata =customdata ,
    hovertemplate ="Tissue=%{y}<br>Metric=%{x}<br>Value=%{z}<br>State=%{customdata[0]}<extra></extra>",
    )
    )
    fig .update_layout (margin =dict (l =20 ,r =20 ,t =20 ,b =20 ),height =max (320 ,90 +28 *len (df )))
    st .plotly_chart (fig ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("promoter_heatmap"))


def render_promoter_chart (rows ,lang ):
    """Display Start sub epigenetic mark counts and composite scores for each Tissue using grouped horizontal bar charts."""
    if go is None :
        st .info ("Plotly is not installed."if lang =="en"else "Plotly is not installed.")
        return 
    if not rows :
        st .info ("No promoter data."if lang =="en"else "No Start sub data.")
        return 

    tissues =[r ["tissue"]for r in rows ]
    open_counts =[r .get ("open_count",0 )for r in rows ]
    active_counts =[r .get ("active_count",0 )for r in rows ]
    repressive_counts =[r .get ("repressive_count",0 )for r in rows ]
    scores =[r .get ("score",0 )for r in rows ]
    states =[r .get ("state","")for r in rows ]

    STATE_COLORS ={
    "active":"#16A34A",
    "repressed":"#DC2626",
    "dynamic_or_bivalent":"#D97706",
    "poised_or_weakly_active":"#2563EB",
    "unclear":"#6B7280",
    }

    fig =go .Figure ()
    fig .add_trace (go .Bar (
    y =tissues ,x =open_counts ,name ="Open"if lang =="en"else "Openness",
    orientation ="h",marker_color ="#0EA5E9",
    ))
    fig .add_trace (go .Bar (
    y =tissues ,x =active_counts ,name ="Active"if lang =="en"else "Activation marker",
    orientation ="h",marker_color ="#22C55E",
    ))
    fig .add_trace (go .Bar (
    y =tissues ,x =repressive_counts ,name ="Repressive"if lang =="en"else "Inhibition marker",
    orientation ="h",marker_color ="#EF4444",
    ))
    fig .update_layout (
    barmode ="group",
    margin =dict (l =20 ,r =20 ,t =30 ,b =20 ),
    height =max (260 ,60 +50 *len (tissues )),
    xaxis_title ="Mark count"if lang =="en"else "Marker count",
    legend =dict (orientation ="h",yanchor ="bottom",y =1.02 ,xanchor ="left",x =0 ),
    )
    st .plotly_chart (fig ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("promoter_chart"))

    #Detailed table per Tissue: score + status + specific marker list
    table_rows =[]
    for r in rows :
        state =r .get ("state","")
        table_rows .append ({
        "Tissue"if lang =="en"else "Tissue":r ["tissue"],
        "Score"if lang =="en"else "Score":round (r .get ("score",0 ),1 ),
        "State"if lang =="en"else "Status":state ,
        "Open marks"if lang =="en"else "Open markers":", ".join (r .get ("open_marks",[]))or "-",
        "Active marks"if lang =="en"else "Activation marker":", ".join (r .get ("active_marks",[]))or "-",
        "Repressive marks"if lang =="en"else "Inhibition marker":", ".join (r .get ("repressive_marks",[]))or "-",
        })
    st .dataframe (pd .DataFrame (table_rows ),use_container_width =True ,hide_index =True )


def render_trait_effect_plot (viz ,lang ):
    if go is None :
        st .info ("Plotly is not installed. GWAS heatmaps are unavailable."if lang =="en"else "Plotly is not installed, cannot display GWAS heatmaps.")
        return 
    core_rows =(viz or {}).get ("trait_effect_plot",[])
    weak_rows =(viz or {}).get ("suggestive_trait_effect_plot",[])
    rows =core_rows +weak_rows 
    if not rows :
        st .info ("No GWAS trait data."if lang =="en"else "No GWAS trait data.")
        return 
    df =pd .DataFrame (rows )
    if "neg_log10_p"not in df .columns :
        df ["neg_log10_p"]=df ["P_Value"].apply (lambda x :-math .log10 (max (float (x ),1e-300 )))
    df ["abs_beta"]=df ["Effect_Size(BETA)"].apply (lambda x :abs (float (x )))
    if "Trait_Normalized_Evidence"not in df .columns :
        df ["Trait_Normalized_Evidence"]="Unclassified"
    if "Trait_Top_Percent_Max"not in df .columns :
        df ["Trait_Top_Percent_Max"]=""
    trait_order =[
    "Heading_date","Plant_height","Num_panicles","Num_effective_panicles",
    "Yield","Grain_weight","Spikelet_length","Grain_length","Grain_width","Grain_thickness"
    ]
    df =df .drop_duplicates (subset =["Trait"],keep ="first").set_index ("Trait",drop =False )
    completed =[]
    gene_label =str (df ["Gene"].iloc [0 ])if "Gene"in df .columns and not df .empty else "Target"
    for trait in trait_order :
        if trait in df .index :
            item =df .loc [trait ].to_dict ()
        else :
            item ={
            "Trait":trait ,
            "Gene":gene_label ,
            "P_Value":float ("nan"),
            "Effect_Size(BETA)":0.0 ,
            "neg_log10_p":0.0 ,
            "abs_beta":0.0 ,
            "Trait_Normalized_Evidence":"Not significant",
            "Trait_Top_Percent_Max":"",
            }
        is_core =str (item .get ("Trait_Normalized_Evidence",""))in ("Strong","Moderate")
        item ["display_abs_beta"]=float (item .get ("abs_beta")or 0.0 )if is_core else 0.0 
        item ["display_neg_log10_p"]=float (item .get ("neg_log10_p")or 0.0 )if is_core else 0.0 
        item ["display_text_beta"]=f"{item ['display_abs_beta']:.3g}"if is_core else "-"
        item ["display_text_p"]=f"{item ['display_neg_log10_p']:.2f}"if is_core else "-"
        completed .append (item )
    df =pd .DataFrame (completed )
    traits =[str (x ).replace ("_"," ")for x in df ["Trait"]]
    gene_label =str (df ["Gene"].iloc [0 ])if "Gene"in df .columns and not df .empty else "Target"

    beta_text =[df ["display_text_beta"].tolist ()]
    p_text =[df ["display_text_p"].tolist ()]
    hover =df [["Trait","Gene","P_Value","Effect_Size(BETA)","neg_log10_p","abs_beta","Trait_Normalized_Evidence","Trait_Top_Percent_Max"]].values .reshape (1 ,len (df ),8 )

    beta_fig =go .Figure (data =go .Heatmap (
    z =[df ["display_abs_beta"].tolist ()],
    x =traits ,
    y =[gene_label ],
    colorscale ="Blues",
    text =beta_text ,
    texttemplate ="%{text}",
    customdata =hover ,
    hovertemplate ="Trait=%{customdata[0]}<br>Gene=%{customdata[1]}<br>Evidence=%{customdata[6]}<br>Trait top <=%{customdata[7]}%<br>|BETA|=%{customdata[5]:.4g}<br>P=%{customdata[2]:.3e}<extra></extra>",
    colorbar =dict (title ="|BETA|"),
    ))
    beta_fig .update_layout (
    title ="GWAS effect-size magnitude (|BETA|)"if lang =="en"else "GWAS 效应量强度（|BETA|）",
    margin =dict (l =20 ,r =20 ,t =45 ,b =65 ),
    height =240 ,
    xaxis =dict (tickangle =35 ),
    )
    st .plotly_chart (beta_fig ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("gwas_abs_beta_heatmap"))

    p_fig =go .Figure (data =go .Heatmap (
    z =[df ["display_neg_log10_p"].tolist ()],
    x =traits ,
    y =[gene_label ],
    colorscale ="YlOrRd",
    text =p_text ,
    texttemplate ="%{text}",
    customdata =hover ,
    hovertemplate ="Trait=%{customdata[0]}<br>Gene=%{customdata[1]}<br>Evidence=%{customdata[6]}<br>Trait top <=%{customdata[7]}%<br>-log10(P)=%{customdata[4]:.4g}<br>P=%{customdata[2]:.3e}<br>|BETA|=%{customdata[5]:.4g}<extra></extra>",
    colorbar =dict (title ="-log10(P)"),
    ))
    p_fig .update_layout (
    title ="GWAS association significance (-log10(P))"if lang =="en"else "GWAS 关联显著性（-log10(P)）",
    margin =dict (l =20 ,r =20 ,t =45 ,b =65 ),
    height =240 ,
    xaxis =dict (tickangle =35 ),
    )
    st .plotly_chart (p_fig ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("gwas_pvalue_heatmap"))


def render_coexpression_network (viz ,lang ):
    if go is None :
        st .info ("Plotly is not installed. Network plot is unavailable."if lang =="en"else "Plotly is not installed, cannot display co-expression network plot.")
        return 
    network =(viz or {}).get ("coexpression_network",{})
    nodes =network .get ("nodes",[])
    edges =network .get ("edges",[])
    if not nodes :
        st .info ("No co-expression network data."if lang =="en"else "No co-expression network data.")
        return 
    target =next ((n for n in nodes if n .get ("role")=="target"),nodes [0 ])
    partners =[n for n in nodes if n .get ("role")!="target"]
    positions ={target ["id"]:(0.0 ,0.0 )}
    n =max (1 ,len (partners ))
    radius =1.0 
    for idx ,node in enumerate (partners ):
        angle =2 *3.1415926 *idx /n 
        positions [node ["id"]]=(radius *math .cos (angle ),radius *math .sin (angle ))

    fig =go .Figure ()
    for edge in edges :
        if edge ["source"]not in positions or edge ["target"]not in positions :
            continue 
        x0 ,y0 =positions [edge ["source"]]
        x1 ,y1 =positions [edge ["target"]]
        fig .add_trace (
        go .Scatter (
        x =[x0 ,x1 ],
        y =[y0 ,y1 ],
        mode ="lines",
        line =dict (color ="#CBD5E1",width =max (1.0 ,abs (float (edge .get ("weight",0.0 )))*6.0 )),
        hoverinfo ="skip",
        showlegend =False ,
        )
        )

    for node in nodes :
        x ,y =positions [node ["id"]]
        theme =node .get ("theme","development_architecture")
        fig .add_trace (
        go .Scatter (
        x =[x ],
        y =[y ],
        mode ="markers+text",
        text =[node .get ("label",node ["id"])],
        textposition ="top center",
        marker =dict (
        size =24 if node .get ("role")=="target"else 18 ,
        color =THEME_COLORS .get (theme ,"#475569"),
        line =dict (color ="#0F172A",width =1.2 ),
        ),
        customdata =[[node ["id"],THEME_LABELS .get (theme ,theme ),node .get ("correlation",0.0 )]],
        hovertemplate ="ID=%{customdata[0]}<br>Theme=%{customdata[1]}<br>Correlation=%{customdata[2]}<extra></extra>",
        showlegend =False ,
        )
        )
    fig .update_layout (
    margin =dict (l =20 ,r =20 ,t =20 ,b =20 ),
    height =420 ,
    xaxis =dict (visible =False ),
    yaxis =dict (visible =False ),
    plot_bgcolor ="white",
    )
    st .plotly_chart (fig ,use_container_width =True ,config ={"displayModeBar":False },key =_plot_key ("coexpression_network"))


def render_protein_result (result_key ):
    payload =st .session_state .get (result_key )
    if not payload :
        return 

    class ResultProxy :
        def __init__ (self ,data ):
            self .returncode =data ["returncode"]
            self .stdout =data ["stdout"]
            self .stderr =data ["stderr"]

    proxy =ResultProxy (payload )
    render_runtime_panel (payload ["cmd"],proxy )

    protein_payload =extract_protein_payload (proxy .stdout )
    if not protein_payload :
        render_full_report (proxy )
        return 

    lang =get_lang ()
    viz =protein_payload .get ("visualizations",{})

    #Only rendering evidence radar chart (all step cards have been rendered in progressive output).
    render_nav_anchor ("evidence-radar")
    with st .container (border =True ):
        st .subheader ("Evidence Radar"if lang =="en"else "Evidence radar chart")
        render_evidence_radar (viz ,lang )

    if proxy .stderr .strip ():
        with st .container (border =True ):
            st .subheader (t ("error_output"))
            st .code (normalize_display_text (proxy .stderr ,lang ))


def render_batch_result (result_key ):
    payload =st .session_state .get (result_key )
    if not payload :
        return 

    class ResultProxy :
        def __init__ (self ,data ):
            self .returncode =data ["returncode"]
            self .stdout =data ["stdout"]
            self .stderr =data ["stderr"]

    proxy =ResultProxy (payload )
    render_runtime_panel (payload ["cmd"],proxy )

    batch_payload =extract_batch_payload (proxy .stdout )
    if not batch_payload :
        with st .container (border =True ):
            st .subheader ("Batch Report Parse Error")
            st .error ("Structured batch output could not be parsed. Showing raw report below.")
        render_full_report (proxy )
        return 

    qc =batch_payload .get ("qc_report",{})
    scenario_rows =batch_payload .get ("scenario_table",[])
    rank_rows =batch_payload .get ("rank_table",[])
    assignment_rows =batch_payload .get ("assignment_matrix",[])
    selected_scenario_tables =batch_payload .get ("selected_scenario_tables",{})
    llm_summary =batch_payload .get ("llm_summary","")
    cluster_reports =batch_payload .get ("cluster_reports",[])
    benchmark =batch_payload .get ("benchmark",{})
    llm_enabled =llm_enabled_for_ui ()

    render_nav_anchor ("batch-qc")
    with st .container (border =True ):
        st .subheader ("QC Summary")
        qc_view ={
        "input_count":qc .get ("input_count",0 ),
        "unique_count":qc .get ("unique_count",0 ),
        "recognized_count":qc .get ("recognized_count",0 ),
        "unrecognized_count":qc .get ("unrecognized_count",0 ),
        "duplicate_count":qc .get ("duplicate_count",0 ),
        }
        st .dataframe (pd .DataFrame ([qc_view ]),use_container_width =True ,hide_index =True )
        if qc .get ("unrecognized_ids"):
            st .caption ("Unrecognized IDs")
            st .code (", ".join (qc ["unrecognized_ids"][:100 ]))

    render_nav_anchor ("batch-scenario-table")
    with st .container (border =True ):
        st .subheader ("Scenario Table")
        if scenario_rows :
            st .dataframe (pd .DataFrame (scenario_rows ),use_container_width =True ,hide_index =True )
        else :
            st .info ("No scenario rows.")

    render_nav_anchor ("batch-assignment-matrix")
    with st .container (border =True ):
        st .subheader ("Scenario Assignment Matrix")
        if assignment_rows :
            assign_df =pd .DataFrame (assignment_rows )
            st .dataframe (assign_df ,use_container_width =True ,hide_index =True )
        else :
            st .info ("No assignment matrix rows.")

    if selected_scenario_tables :
        render_nav_anchor ("batch-candidate-tables")
        with st .container (border =True ):
            st .subheader ("Selected Scenario Candidate Tables")
            for scenario ,rows in selected_scenario_tables .items ():
                with st .expander (scenario ,expanded =True ):
                    if rows :
                        df =pd .DataFrame (rows )
                        preferred_cols =[
                        "gene_id",
                        "primary_symbol",
                        "assignment_score",
                        "priority_score",
                        "tier",
                        "evidence_types",
                        "traits",
                        "modules",
                        "pmid_count",
                        "risk_tag",
                        "description",
                        ]
                        cols =[c for c in preferred_cols if c in df .columns ]
                        st .dataframe (df [cols ]if cols else df ,use_container_width =True ,hide_index =True )
                    else :
                        st .info ("No candidate genes passed the current assignment threshold.")

    render_nav_anchor ("batch-rank-table")
    with st .container (border =True ):
        st .subheader ("Rank Table")
        if rank_rows :
            rank_df =pd .DataFrame (rank_rows )
            preferred_cols =[
            "rank",
            "gene_id",
            "primary_symbol",
            "total_score",
            "main_scenario",
            "func_score",
            "network_score",
            "gwas_score",
            "literature_boost_score",
            "risk_tag",
            ]
            cols =[c for c in preferred_cols if c in rank_df .columns ]
            st .dataframe (rank_df [cols ]if cols else rank_df ,use_container_width =True ,hide_index =True )
        else :
            st .info ("No rank rows.")

    render_nav_anchor ("batch-integrated-summary")
    with st .container (border =True ):
        st .subheader ("Integrated Summary")
        if not llm_enabled :
            render_api_required_hint (get_lang ())
        else :
            st .markdown (llm_summary if llm_summary else "No summary.")

    if cluster_reports or not llm_enabled :
        with st .container (border =True ):
            st .subheader ("DEG Cluster Reports")
            if not llm_enabled :
                render_api_required_hint (get_lang ())
            else :
                for report in cluster_reports :
                    title =report .get ("cluster_name","Unnamed Cluster")
                    genes =report .get ("genes",[])
                    gene_line =", ".join (genes [:30 ])if genes else "N/A"
                    pmid_count =report .get ("abstracts_pmid_count",0 )
                    with st .expander (f"{title } | genes={report .get ('gene_count',0 )} | pmids={pmid_count }",expanded =False ):
                        st .caption (f"Genes: {gene_line }")
                        st .markdown (report .get ("summary","No cluster summary."))

    if benchmark .get ("available"):
        with st .container (border =True ):
            st .subheader ("Benchmark")
            metric_rows =[]
            for name ,item in benchmark .get ("metrics",{}).items ():
                metric_rows .append ({"set":name ,**item })
            if metric_rows :
                st .dataframe (pd .DataFrame (metric_rows ),use_container_width =True ,hide_index =True )

    if proxy .stderr .strip ():
        with st .container (border =True ):
            st .subheader (t ("error_output"))
            st .code (normalize_display_text (proxy .stderr ,get_lang ()))


def top_runtime_config ():
    render_nav_anchor ("runtime-config")
    st .session_state ["ui_dark_mode"]=False 

    llm_cfg =get_effective_llm_config ()
    llm_enabled =bool (llm_cfg ["api_key"])
    st .session_state ["llm_enabled"]=llm_enabled

    llm_cfg =get_effective_llm_config ()
    llm_enabled =bool (llm_cfg ["api_key"])
    extra_env ={
    "APP_LANG":"en",
    "LLM_BASE_URL":llm_cfg ["base_url"],
    "LLM_API_KEY":llm_cfg ["api_key"],
    "LLM_MODEL_NAME":llm_cfg ["model_name"],
    "BATCH_ENABLE_LLM_SUMMARY":"1"if llm_enabled else "0",
    "RIGOR_LLM_ENABLED":"1"if llm_enabled else "0",
    }
    return extra_env 


@st .cache_resource 
def get_gene_annotator ():
    gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")
    map_file =str (ANNOTATION_DIR /"RAP-MSU.txt")
    annotator =RiceGeneAnnotator (gff_file ,map_file )
    annotator .build_database ()
    return annotator 


@st .cache_resource 
def get_protein_feature_service ():
    return ProteinFeatureService ()


@st .cache_resource 
def get_expression_pattern_service ():
    return ExpressionPatternService ()


@st .cache_data (show_spinner =False )
def get_protein_feature_bundle (rap_id ):
    annotator =get_gene_annotator ()
    gene_info =annotator .query (rap_id )
    if not gene_info :
        return None 
    protein_service =get_protein_feature_service ()
    expression_service =get_expression_pattern_service ()
    return {
    "gene_info":gene_info ,
    "protein_features":protein_service .collect (gene_info ),
    "expression_patterns":expression_service .collect (gene_info ),
    }


def render_protein_feature_panel (gene_id ):
    lang =get_lang ()
    if not gene_id or not gene_id .strip ():
        return 

    bundle =get_protein_feature_bundle (gene_id .strip ())
    if not bundle :
        st .warning ("Gene not found. Protein features unavailable."if lang =="en"else "Gene not found, cannot display protein domain and subcellular localization.")
        return 

    gene_info =bundle ["gene_info"]
    protein =bundle ["protein_features"]
    expression =bundle ["expression_patterns"]

    with st .container (border =True ):
        st .subheader ("Protein Annotations"if lang =="en"else "Protein annotation supplement")
        info_cols =st .columns (4 )
        info_cols [0 ].metric ("RAP ID",gene_info .get ("RAP_ID","N/A"))
        info_cols [1 ].metric ("UniProt",protein .get ("uniprot_id")or "N/A")
        info_cols [2 ].metric ("Protein Length"if lang =="en"else "Protein length",f"{protein ['protein_length']} aa"if protein .get ("protein_length")else "N/A")
        info_cols [3 ].metric ("Pfam Domains"if lang =="en"else "Number of Pfam domains",len (protein .get ("pfam_domains",[])))

        subcell =protein .get ("subcellular_locations",[])
        if subcell :
            st .caption ("Subcellular localization"if lang =="en"else "Subcellular localization")
            st .write (" | ".join (subcell ))
        elif protein .get ("local_subcellular")not in ("","N/A"):
            st .caption ("Subcellular localization (local annotation)"if lang =="en"else "Subcellular localization (local annotation)")
            st .write (normalize_display_text (protein ["local_subcellular"],lang ))
        else :
            st .info ("UniProt localization is unavailable, or network is inaccessible."if lang =="en"else "UniProt did not return subcellular localization annotation, or current network is unavailable.")

        if protein .get ("warnings"):
            for warning in protein ["warnings"]:
                st .warning (normalize_display_text (warning ,lang ))

        interpro_domains =protein .get ("interpro_domains",[])
        if interpro_domains :
            st .caption ("InterPro domain details"if lang =="en"else "InterPro domain details")
            st .dataframe (pd .DataFrame (interpro_domains ),use_container_width =True ,hide_index =True )
        elif protein .get ("local_interpro")not in ("","N/A"):
            st .caption ("InterPro domains (local annotation)"if lang =="en"else "InterPro domains (local annotation)")
            st .write (normalize_display_text (protein ["local_interpro"],lang ))
        else :
            st .info ("No InterPro domain details available."if lang =="en"else "InterPro domain details could not be retrieved.")

        pfam_domains =protein .get ("pfam_domains",[])
        if pfam_domains :
            st .caption ("Pfam core domains"if lang =="en"else "Pfam Core Domain")
            st .dataframe (pd .DataFrame (pfam_domains ),use_container_width =True ,hide_index =True )

        if protein .get ("plot_png"):
            st .caption ("Pfam domain plot"if lang =="en"else "Pfam Domain plot")
            st .image (protein ["plot_png"],use_container_width =True )
        else :
            st .info ("Domain plot is unavailable."if lang =="en"else "Domain plot has not been generated yet.")

    with st .container (border =True ):
        st .subheader ("Transcript Expression Patterns"if lang =="en"else "TranscriptExpression pattern")
        transcript_ids =expression .get ("transcript_ids",[])or gene_info .get ("Transcript_IDs",[])
        expression_model_ids =expression .get ("expression_model_ids",[])
        if lang =="en":
            st .caption (f"Gene {gene_info .get ('RAP_ID','N/A')}: {len (transcript_ids )} annotated transcript(s)")
        else :
            st .caption (f"Gene {gene_info .get ('RAP_ID','N/A')} co-discovered {len (transcript_ids )} eachGenegroupAnnotationTranscript")
        if transcript_ids :
            st .write (" | ".join (transcript_ids ))
        st .caption (
        f"Matched MSU models in expression matrix: {len (expression_model_ids )}"
        if lang =="en"
        else f"Expressionmatrix-matching MSU model count: {len (expression_model_ids )}"
        )
        if expression_model_ids :
            st .write (" | ".join (expression_model_ids ))
        st .caption (
        f"Expression matrix: {expression .get ('matrix_path','N/A')}"
        if lang =="en"
        else f"ExpressionMatrixFile: {expression .get ('matrix_path','N/A')}"
        )

        if expression .get ("warnings"):
            for warning in expression ["warnings"]:
                st .warning (normalize_display_text (warning ,lang ))

        for transcript_result in expression .get ("transcript_results",[]):
            with st .container (border =True ):
                st .markdown (
                f"**Transcript {transcript_result ['transcript_id']}**"
                if lang =="en"
                else f"**Transcript {transcript_result ['transcript_id']}**"
                )
                st .caption (
                f"Expression model: {transcript_result .get ('expression_model_id','N/A')}"
                if lang =="en"
                else f"Expressionmatrix model: {transcript_result .get ('expression_model_id','N/A')}"
                )
                for warning in transcript_result .get ("warnings",[]):
                    st .warning (normalize_display_text (warning ,lang ))

                col1 ,col2 =st .columns (2 )
                with col1 :
                    st .caption ("Organ/time expression pattern"if lang =="en"else "Different organ/time expression patterns")
                    if transcript_result .get ("developmental_plot_png"):
                        st .image (transcript_result ["developmental_plot_png"],use_container_width =True )
                    else :
                        st .info ("No organ/time expression plot."if lang =="en"else "No organism/temporal expression plot generated currently.")
                    if transcript_result .get ("developmental_summary"):
                        st .dataframe (pd .DataFrame (transcript_result ["developmental_summary"]),use_container_width =True ,hide_index =True )

                with col2 :
                    st .caption ("Treatment expression pattern"if lang =="en"else "Different TreatmentExpression patterns")
                    if transcript_result .get ("treatment_plot_png"):
                        st .image (transcript_result ["treatment_plot_png"],use_container_width =True )
                    else :
                        st .info ("No treatment expression plot."if lang =="en"else "No treatment expression plot generated currently.")
                    if transcript_result .get ("treatment_summary"):
                        st .dataframe (pd .DataFrame (transcript_result ["treatment_summary"]),use_container_width =True ,hide_index =True )

def protein_tab (extra_env ):
    lang =get_lang ()
    render_nav_anchor ("protein-input")
    st .subheader ("Protein-coding Gene Function Prediction"if lang =="en"else "Protein coding gene function prediction")
    mode_options =["Gene ID","Coordinate"]if lang =="en"else ["GeneID","Coordinates"]
    mode =st .radio ("Input mode"if lang =="en"else "Input method",mode_options ,horizontal =True ,key ="protein_mode")
    with st .form ("protein_form"):
        if mode ==mode_options [0 ]:
            gene =st .text_input ("Target gene"if lang =="en"else "Target Gene",value ="Os01g0583100")
            st .caption ("MSU/LOC IDs are also supported.")
            submitted =st .form_submit_button ("Run Protein Function"if lang =="en"else "RunningProtein Function Prediction")
            if submitted :
                st .session_state ["protein_last_gene"]=gene .strip ()
                execute_and_store ("protein_result",APP_MAP ["Protein Function Prediction"],["--gene",gene ],extra_env )
        else :
            coord =st .text_input ("Target coordinate"if lang =="en"else "Target coordinates",value ="chr04:18172160-18177058")
            strand =st .selectbox ("Strand"if lang =="en"else "Strand Direction",["+","-"],index =0 )
            submitted =st .form_submit_button ("Run Coordinate Mode"if lang =="en"else "Running coordinate mode")
            if submitted :
                st .session_state ["protein_last_gene"]=""
                execute_and_store ("protein_result",APP_MAP ["Protein Function Prediction"],["--coord",coord ,"--strand",strand ],extra_env )
    render_protein_result ("protein_result")


def _persist_uploaded_csv (uploaded_file ):
    suffix =Path (uploaded_file .name ).suffix or ".csv"
    fd ,temp_path =tempfile .mkstemp (prefix ="batch_genes_",suffix =suffix )
    with os .fdopen (fd ,"wb")as handle :
        handle .write (uploaded_file .getvalue ())
    return temp_path 


def protein_batch_tab (extra_env ):
    lang =get_lang ()
    render_nav_anchor ("batch-input")
    st .subheader ("Batch Protein Analysis"if lang =="en"else "Batch analysis of protein-coding genes")
    st .caption (
    "Two-stage workflow: local scoring + literature boost + LLM synthesis"
    if lang =="en"
    else "Two-stage process: Local scoring + Literature enhancement + LLM comprehensive summary"
    )

    uploaded =st .file_uploader (
    "Upload CSV (single gene_id column)"if lang =="en"else "Upload CSV (single column of gene_id)",
    type =["csv"],
    key ="protein_batch_upload",
    )
    st .caption ("Recommended input size: up to 100 genes."if lang =="en"else "Recommended input size: up to 100 genes.")
    gene_col =st .text_input ("Gene column name"if lang =="en"else "Gene list",value ="gene_id",key ="protein_batch_col")
    selected_scenarios =st .multiselect (
    "Fixed scenarios"if lang =="en"else "Fixed scene selection",
    options =SCENARIO_OPTIONS ,
    default =SCENARIO_OPTIONS ,
    key ="protein_batch_scenarios",
    )

    with st .expander ("Scoring weights"if lang =="en"else "Scoring weights",expanded =False ):
        func_w =st .slider ("func_weight",0.0 ,1.0 ,0.38 ,0.01 ,key ="w_func")
        network_w =st .slider ("network_weight",0.0 ,1.0 ,0.31 ,0.01 ,key ="w_network")
        gwas_w =st .slider ("gwas_weight",0.0 ,1.0 ,0.19 ,0.01 ,key ="w_gwas")
        complete_w =st .slider ("completeness_weight",0.0 ,1.0 ,0.12 ,0.01 ,key ="w_complete")
        top_per =st .slider ("top_per_scenario",1 ,30 ,10 ,1 ,key ="w_top_per")
        max_genes =st .slider ("max_genes",50 ,500 ,300 ,10 ,key ="w_max_genes")

    with st .form ("protein_batch_form"):
        submitted =st .form_submit_button ("Run Batch Analysis"if lang =="en"else "Running batch analysis")
        if submitted :
            if uploaded is None :
                st .error ("Please upload CSV first."if lang =="en"else "Please upload the CSV file first.")
            else :
                temp_csv =_persist_uploaded_csv (uploaded )
                weights ={
                "func_weight":func_w ,
                "network_weight":network_w ,
                "gwas_weight":gwas_w ,
                "completeness_weight":complete_w ,
                "top_per_scenario":top_per ,
                "max_genes":max_genes ,
                "selected_scenarios":selected_scenarios ,
                }
                args =["--input-csv",temp_csv ,"--gene-column",gene_col ,"--weights",json .dumps (weights )]
                execute_and_store ("protein_batch_result",APP_MAP ["Protein Batch Analysis"],args ,extra_env )
    render_batch_result ("protein_batch_result")


def srna_tab (extra_env ):
    lang =get_lang ()
    render_nav_anchor ("srna-input")
    st .subheader ("sRNA Analysis"if lang =="en"else "sRNA Analysis")
    with st .form ("srna_form"):
        inp =st .text_input (
        "miRNA / phasiRNA ID or mature sequence"if lang =="en"else "miRNA/phasiRNA ID or mature sequence",
        value ="osa-mir156a",
        help ="Input ID (e.g., osa-mir156a, Osa-phasi09g0719) or mature sequence (≤25nt)"if lang =="zh"
        else "Enter ID (e.g. osa-mir156a, Osa-phasi09g0719) or mature sequence (≤25nt)"
        )
        top_n =st .slider ("Target gene count"if lang =="en"else "Number of target genes",min_value =1 ,max_value =20 ,value =3 )
        submitted =st .form_submit_button ("Run sRNA Analysis"if lang =="en"else "Running sRNA Analysis")
        if submitted :
            args =["--gene",inp .strip (),"--top-n",str (top_n )]
            execute_and_store ("srna_result",APP_MAP ["sRNA"],args ,extra_env )
    render_stored_result ("srna_result")


def lncrna_tab (extra_env ):
    lang =get_lang ()
    render_nav_anchor ("lncrna-input")
    st .subheader ("lncRNA Analysis"if lang =="en"else "lncRNA Analysis")
    input_mode_label ="Input mode"if lang =="en"else "Input method"
    mode_id ="ID"if lang =="en"else "ID"
    mode_coord ="Coordinate"if lang =="en"else "Genomic coordinates"
    with st .form ("lncrna_form"):
        input_mode =st .radio (input_mode_label ,[mode_id ,mode_coord ],horizontal =True )
        if input_mode ==mode_coord :
            coord_help =(
            "Format: CHR1,12321312-12321999  or  chr1:12321312-12321999"
            if lang =="en"
            else "Format: chr1,12321312-12321999 or chr1:12321312-12321999"
            )
            coord_str =st .text_input (
            "Coordinate"if lang =="en"else "Coordinates",
            value ="chr1,766902-774672",
            help =coord_help ,
            )
            min_overlap =st .slider (
            "Min. overlap similarity (%)"if lang =="en"else "Minimum overlap similarity (%)",
            min_value =10 ,max_value =100 ,value =80 ,step =5 ,
            help =(
            "Overlap length / shorter interval length. Default 80%."
            if lang =="en"
            else "Overlap length / shorter interval length. Default 80%."
            ),
            )
            submitted =st .form_submit_button (
            "Search by Coordinate"if lang =="en"else "Search by coordinates"
            )
            if submitted :
                execute_and_store (
                "lncrna_result",
                APP_MAP ["lncRNA"],
                ["--coord",coord_str .strip (),"--min-overlap",str (min_overlap /100 )],
                extra_env ,
                )
        else :
            lncrna_id =st .text_input ("lncRNA ID",value ="OSA-LNCN0041")
            submitted =st .form_submit_button (
            "Run lncRNA Analysis"if lang =="en"else "Running lncRNA Analysis"
            )
            if submitted :
                execute_and_store ("lncrna_result",APP_MAP ["lncRNA"],["--gene",lncrna_id ],extra_env )
    render_stored_result ("lncrna_result")


def identifier_tab (extra_env ):
    lang =get_lang ()
    render_nav_anchor ("identifier-input")
    st .subheader ("RNA Type Identification"if lang =="en"else "RNA type identification")
    with st .form ("identifier_form"):
        rna_id =st .text_input ("RNA ID (optional)"if lang =="en"else "RNA ID (optional)",value ="")
        sequence =st .text_input ("RNA sequence (optional)"if lang =="en"else "RNA Sequence (optional)",value ="")
        submitted =st .form_submit_button ("Identify Type"if lang =="en"else "Identified type")
        if submitted :
            args =[]
            if rna_id .strip ():
                args .extend (["--id",rna_id .strip ()])
            if sequence .strip ():
                args .extend (["--sequence",sequence .strip ()])
            execute_and_store ("identifier_result",APP_MAP ["RNA Type Identification"],args ,extra_env )
    render_stored_result ("identifier_result")


def help_tab ():
    render_nav_anchor ("help-overview")
    help_path =PROJECT_ROOT /"docs"/"RIGOR_help.md"
    if not help_path .exists ():
        st .warning ("Help document is not available.")
        return 
    content =help_path .read_text (encoding ="utf-8")
    content =content .replace ("## 2. Main Modules",'<div id="rigor-help-modules" class="rigor-anchor"></div>\n\n## 2. Main Modules')
    content =content .replace ("## 3. Web Page Usage",'<div id="rigor-help-usage" class="rigor-anchor"></div>\n\n## 3. Web Page Usage')
    content =content .replace ("## 4. Interpreting Results",'<div id="rigor-help-results" class="rigor-anchor"></div>\n\n## 4. Interpreting Results')
    content =content .replace ("## 5. FAQ",'<div id="rigor-help-faq" class="rigor-anchor"></div>\n\n## 5. FAQ')
    st .markdown (content ,unsafe_allow_html =True )


def main ():
    global NAV_PLACEHOLDER 
    st .session_state ["ui_dark_mode"]=False 
    is_dark =False 

    st .markdown ("""
<style>
/* Premium visual theme */
.stApp {
    background:
      radial-gradient(1100px 500px at 8% -10%, rgba(16,185,129,0.14), transparent 60%),
      radial-gradient(900px 450px at 100% 0%, rgba(14,165,233,0.14), transparent 62%),
      linear-gradient(180deg, #F8FAFC 0%, #F1F5F9 100%);
}
.block-container {
    padding-top: 0.25rem;
    padding-bottom: 2rem;
    max-width: 1350px;
    transition: padding-left 0.22s ease;
}
@media (min-width: 1050px) {
    .block-container {
        padding-left: 1rem;
    }
    .stApp:has(.rigor-nav-toggle:checked) .block-container {
        padding-left: 15.35rem;
    }
}
.rigor-anchor {
    display: block;
    height: 0;
    position: relative;
    top: -1rem;
    visibility: hidden;
}
[data-testid="stSidebar"] {
    border-right: 1px solid rgba(148,163,184,0.32);
}
[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] a {
    color: #0F766E;
    text-decoration: none;
    font-weight: 600;
}
[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] a:hover {
    color: #0891B2;
    text-decoration: underline;
}
.rigor-nav-toggle {
    position: fixed;
    opacity: 0;
    pointer-events: none;
}
.rigor-nav-toggle-label {
    position: fixed;
    left: 0;
    top: 7.2rem;
    z-index: 61;
    display: flex;
    align-items: center;
    justify-content: center;
    width: 32px;
    height: 38px;
    padding: 0;
    border-radius: 0 10px 10px 0;
    background: linear-gradient(135deg, #0F766E 0%, #0891B2 100%);
    color: #FFFFFF;
    font-size: 1rem;
    font-weight: 800;
    cursor: pointer;
    box-shadow: none;
    transition: left 0.22s ease, opacity 0.18s ease;
}
.rigor-nav-toggle-label::after {
    content: "›";
    font-size: 1.45rem;
    line-height: 1;
    transform: translateY(-1px);
}
.rigor-local-nav {
    position: fixed;
    left: 0;
    top: 0;
    width: 232px;
    z-index: 50;
    height: 100vh;
    max-height: 100vh;
    overflow-y: auto;
    padding: 7.1rem 0.82rem 1.25rem 0.82rem;
    border-radius: 0;
    background: #FFFFFF;
    border: none;
    border-right: 1px solid #CBD5E1;
    box-shadow: none;
    backdrop-filter: none;
    opacity: 0;
    transform: translateX(calc(-100% - 1.3rem));
    pointer-events: none;
    transition: transform 0.22s ease, opacity 0.18s ease;
}
.rigor-nav-toggle:checked ~ .rigor-local-nav {
    opacity: 1;
    transform: translateX(0);
    pointer-events: auto;
}
.rigor-nav-toggle:checked ~ .rigor-nav-toggle-label {
    left: 232px;
}
.rigor-nav-toggle:checked ~ .rigor-nav-toggle-label::after {
    content: "‹";
}
.rigor-nav-tabs {
    display: none;
}
.rigor-nav-tab {
    display: block;
    color: #475569;
    background: #F1F5F9;
    border: 1px solid #CBD5E1;
    border-radius: 10px;
    padding: 0.44rem 0.54rem;
    font-size: 0.82rem;
    font-weight: 800;
    line-height: 1.08;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}
.rigor-local-nav-title {
    color: #0F172A;
    font-size: 0.86rem;
    font-weight: 800;
    line-height: 1.2;
    margin: 0.16rem 0.14rem 0.56rem 0.14rem;
    letter-spacing: 0.02em;
}
.rigor-local-nav-link {
    display: block;
    color: #0F766E !important;
    text-decoration: none !important;
    font-size: 0.88rem;
    font-weight: 700;
    line-height: 1.15;
    min-height: 36px;
    padding: 0.48rem 0.56rem;
    border-radius: 11px;
    margin-bottom: 0.2rem;
    background: #F0FDFA;
    border: 1px solid #CCFBF1;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}
.rigor-local-nav-link:hover {
    color: #FFFFFF !important;
    background: linear-gradient(135deg, #0F766E 0%, #0891B2 100%);
}
.rigor-nav-group {
    display: none;
}
.rigor-nav-group-protein {
    display: block;
}
.stApp:has(.stTabs [data-baseweb="tab"][aria-selected="true"]) .rigor-nav-group {
    display: none;
}
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(1)[aria-selected="true"]) .rigor-nav-group-protein,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(2)[aria-selected="true"]) .rigor-nav-group-batch,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(3)[aria-selected="true"]) .rigor-nav-group-srna,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(4)[aria-selected="true"]) .rigor-nav-group-lncrna,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(5)[aria-selected="true"]) .rigor-nav-group-identifier,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(6)[aria-selected="true"]) .rigor-nav-group-help {
    display: block;
}
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(1)[aria-selected="true"]) .rigor-nav-tab-protein,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(2)[aria-selected="true"]) .rigor-nav-tab-batch,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(3)[aria-selected="true"]) .rigor-nav-tab-srna,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(4)[aria-selected="true"]) .rigor-nav-tab-lncrna,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(5)[aria-selected="true"]) .rigor-nav-tab-identifier,
.stApp:has(.stTabs [data-baseweb="tab"]:nth-of-type(6)[aria-selected="true"]) .rigor-nav-tab-help {
    color: #FFFFFF;
    background: linear-gradient(135deg, #0F766E 0%, #0891B2 100%);
    border-color: transparent;
}
@media (max-width: 900px) {
    .rigor-nav-toggle,
    .rigor-nav-toggle-label {
        display: none;
    }
    .rigor-local-nav {
        position: static;
        width: auto;
        opacity: 1;
        transform: none;
        pointer-events: auto;
        max-height: none;
        padding: 0.55rem;
        margin-bottom: 0.75rem;
        border: 1px solid #CBD5E1;
        border-radius: 12px;
    }
    .rigor-nav-tabs {
        grid-template-columns: repeat(2, minmax(0, 1fr));
    }
}
.rigor-top-logo {
    margin-bottom: 0.8rem;
}
.rigor-loader-wrap {
    position: relative;
    width: 100%;
    min-height: 126px;
    border-radius: 16px;
    background: linear-gradient(135deg, rgba(255,255,255,0.92), rgba(241,245,249,0.92));
    border: 1px solid #CBD5E1;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
    gap: 8px;
    overflow: hidden;
    margin: 0.35rem 0 0.55rem 0;
}
.rigor-loader-stage {
    position: relative;
    width: 110px;
    height: 78px;
}
.rigor-loader-orb {
    position: absolute;
    top: 50%;
    left: 50%;
    width: 26px;
    height: 26px;
    border-radius: 50%;
    background: radial-gradient(circle at 30% 30%, #5EEAD4 0%, #0EA5E9 48%, #0F766E 100%);
    box-shadow: 0 0 0 0 rgba(14,165,233,0.55);
    animation: rigor-orb-breathe 1.5s ease-in-out infinite;
    transform: translate(-50%, -50%);
    z-index: 2;
}
.rigor-loader-text {
    margin: 0;
    font-size: 0.95rem;
    color: #0F172A;
    font-weight: 600;
    letter-spacing: 0.2px;
    text-align: center;
    width: 100%;
    z-index: 2;
}
@keyframes rigor-orb-breathe {
    0% { transform: translate(-50%, -50%) scale(0.9); box-shadow: 0 0 0 0 rgba(14,165,233,0.52); }
    70% { transform: translate(-50%, -50%) scale(1.08); box-shadow: 0 0 0 20px rgba(14,165,233,0); }
    100% { transform: translate(-50%, -50%) scale(0.9); box-shadow: 0 0 0 0 rgba(14,165,233,0); }
}
[data-testid="stHeader"] {
    display: none;
}
[data-testid="stToolbar"] {
    display: none;
}
[data-testid="stDecoration"] {
    display: none;
}
#MainMenu {
    visibility: hidden;
}
footer {
    visibility: hidden;
}
.rigor-app-footer {
    margin-top: 2.4rem;
    padding: 1.1rem 0 0.25rem;
    border-top: 1px solid rgba(148,163,184,0.35);
    color: #475569;
    font-size: 0.86rem;
    line-height: 1.55;
    text-align: center;
}
.rigor-app-footer a {
    color: #0F766E;
    font-weight: 650;
    text-decoration: none;
}
.rigor-app-footer a:hover {
    color: #0891B2;
    text-decoration: underline;
}
[data-testid="stVerticalBlock"] > div[data-testid="stHorizontalBlock"] {
    gap: 0.9rem;
}
[data-testid="stForm"] {
    background: rgba(255,255,255,0.88);
    border: 1px solid rgba(148,163,184,0.22);
    border-radius: 16px;
    padding: 1rem 1rem 0.7rem 1rem;
    backdrop-filter: blur(6px);
}
[data-testid="stMetric"] {
    background: #FFFFFF;
    border: 1px solid #E2E8F0;
    border-radius: 12px;
    padding: 0.35rem 0.55rem;
}
[data-testid="stButton"] > button, [data-testid="baseButton-secondary"] {
    border-radius: 12px !important;
}
[data-testid="stFileUploader"] section {
    border-radius: 12px !important;
    border: 1px dashed #94A3B8 !important;
}
/* Tabs */
.stTabs [data-baseweb="tab-list"] {
    gap: 10px;
    display: grid;
    grid-template-columns: repeat(6, minmax(0, 1fr));
    width: 100%;
    background: rgba(255,255,255,0.92);
    border: 1px solid #CBD5E1;
    border-bottom: none;
    border-radius: 16px 16px 0 0;
    padding: 10px;
    margin-bottom: 0;
}
.stTabs [data-baseweb="tab"] {
    width: 100%;
    min-height: 52px;
    font-size: 1.02rem;
    font-weight: 600;
    padding: 10px 20px;
    border-radius: 12px;
    background-color: #E2E8F0;
    border: 1px solid #CBD5E1;
}
.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, #0F766E 0%, #0891B2 100%);
    color: white !important;
    border-color: transparent;
}
.stTabs [data-baseweb="tab-panel"] {
    background: rgba(255,255,255,0.92);
    border: 1px solid #CBD5E1;
    border-top: none;
    border-radius: 0 0 16px 16px;
    padding: 1rem 0.9rem 1.2rem 0.9rem;
    margin-top: 0;
}
</style>
""",unsafe_allow_html =True )

    if is_dark :
        st .markdown ("""
<style>
.stApp {
    background:
      radial-gradient(1100px 520px at 10% -8%, rgba(16,185,129,0.18), transparent 58%),
      radial-gradient(980px 460px at 98% 0%, rgba(14,165,233,0.20), transparent 60%),
      linear-gradient(180deg, #0B1220 0%, #0F172A 100%) !important;
}
[data-testid="stAppViewContainer"], [data-testid="stAppViewContainer"] * {
    color: #E2E8F0;
}
[data-testid="stForm"] {
    background: rgba(15,23,42,0.70) !important;
    border: 1px solid #334155 !important;
}
[data-testid="stMetric"] {
    background: rgba(15,23,42,0.76) !important;
    border: 1px solid #334155 !important;
}
.stTabs [data-baseweb="tab-list"] {
    background: rgba(15,23,42,0.86) !important;
    border: 1px solid #334155 !important;
}
.stTabs [data-baseweb="tab"] {
    background-color: #1E293B !important;
    border: 1px solid #334155 !important;
    color: #CBD5E1 !important;
}
.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, #0F766E 0%, #0EA5E9 100%) !important;
    color: #F8FAFC !important;
}
.stTabs [data-baseweb="tab-panel"] {
    background: rgba(15,23,42,0.78) !important;
    border: 1px solid #334155 !important;
    border-top: none !important;
}
.rigor-loader-wrap {
    background: linear-gradient(135deg, rgba(15,23,42,0.92), rgba(30,41,59,0.90)) !important;
    border: 1px solid #334155 !important;
}
.rigor-loader-text {
    color: #E2E8F0 !important;
}
.rigor-local-nav {
    background: #0F172A !important;
    border: none !important;
    border-right: 1px solid #334155 !important;
    box-shadow: none !important;
}
.rigor-nav-toggle-label {
    background: linear-gradient(135deg, #0F766E 0%, #0EA5E9 100%) !important;
    color: #F8FAFC !important;
}
.rigor-nav-tab {
    color: #CBD5E1 !important;
    background: #1E293B !important;
    border-color: #475569 !important;
}
.rigor-local-nav-title {
    color: #E2E8F0 !important;
}
.rigor-local-nav-link {
    color: #5EEAD4 !important;
    background: #134E4A !important;
    border-color: #0F766E !important;
}
/* Keep input text and action buttons non-inverted in dark mode */
[data-testid="stTextInput"] input,
[data-testid="stTextArea"] textarea,
[data-testid="stNumberInput"] input {
    color: #0F172A !important;
    background: #FFFFFF !important;
    caret-color: #0F172A !important;
}
[data-testid="stButton"] > button,
[data-testid="stFormSubmitButton"] > button,
[data-baseweb="button"],
[data-testid="baseButton-secondary"],
[data-testid="baseButton-primary"] {
    color: #0F172A !important;
    background: #F8FAFC !important;
    border: 1px solid #CBD5E1 !important;
    opacity: 1 !important;
    -webkit-text-fill-color: #0F172A !important;
}
[data-testid="stButton"] > button *,
[data-testid="stFormSubmitButton"] > button *,
[data-baseweb="button"] *,
[data-testid="baseButton-secondary"] *,
[data-testid="baseButton-primary"] * {
    color: #0F172A !important;
    fill: #0F172A !important;
    -webkit-text-fill-color: #0F172A !important;
}
[data-testid="stButton"] > button:hover,
[data-testid="stFormSubmitButton"] > button:hover,
[data-baseweb="button"]:hover,
[data-testid="baseButton-secondary"]:hover,
[data-testid="baseButton-primary"]:hover {
    color: #0F172A !important;
    background: #EEF2F7 !important;
}
[data-testid="stButton"] > button:disabled,
[data-testid="stFormSubmitButton"] > button:disabled,
[data-baseweb="button"]:disabled {
    color: #0F172A !important;
    opacity: 1 !important;
    -webkit-text-fill-color: #0F172A !important;
}
/* File uploader text: force black and replace size hint behavior */
[data-testid="stFileUploader"] section,
[data-testid="stFileUploader"] section * {
    color: #0F172A !important;
    -webkit-text-fill-color: #0F172A !important;
}
[data-testid="stFileUploader"] small {
    display: none !important;
}
</style>
""",unsafe_allow_html =True )

    render_nav_anchor ("top")
    st .markdown ('<div class="rigor-top-logo"></div>',unsafe_allow_html =True )
    render_logo_row ()
    extra_env =top_runtime_config ()
    NAV_PLACEHOLDER =st .empty ()
    render_global_fixed_nav (NAV_PLACEHOLDER )

    tab1 ,tab2 ,tab3 ,tab4 ,tab5 ,tab6 =st .tabs (
    [t ("tab_protein"),t ("tab_batch"),t ("tab_srna"),t ("tab_lncrna"),t ("tab_identifier"),t ("tab_help")]
    )

    with tab1 :
        protein_tab (extra_env )
    with tab2 :
        protein_batch_tab (extra_env )
    with tab3 :
        srna_tab (extra_env )
    with tab4 :
        lncrna_tab (extra_env )
    with tab5 :
        identifier_tab (extra_env )
    with tab6 :
        help_tab ()

    render_app_footer ()


if __name__ =="__main__":
    main ()
