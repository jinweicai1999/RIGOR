from __future__ import annotations 

from io import BytesIO 
import os 
from pathlib import Path 
from typing import Any 

import matplotlib 
matplotlib .use ("Agg")
import matplotlib .pyplot as plt 
import pandas as pd 

from src .project_paths import EXPRESSION_DIR 


SAMPLE_DESCRIPTION_MAP ={
"SRR1618546":"Anthers",
"SRR1618547":"Pistils",
"SRR1618548":"Seeds-5DAP",
"SRR1618549":"Shoots-14 DAG",
"SRR8490277":"Seedlings-Control-Rep1",
"SRR8490276":"Seedlings-Control-Rep2",
"SRR8490271":"Seedlings-Control-Rep3",
"SRR8490275":"Seedlings-Drought-Rep1",
"SRR8490274":"Seedlings-Drought-Rep2",
"SRR8490281":"Seedlings-Drought-Rep3",
"SRR7059297":"Mature seed-Control-Rep1",
"SRR7059298":"Mature seed-Control-Rep2",
"SRR7059299":"Mature seed-Control-Rep3",
"SRR7059300":"Mature seed-Hypoxia-Rep1",
"SRR7059295":"Mature seed-Hypoxia-Rep2",
"SRR7059296":"Mature seed-Hypoxia-Rep3",
"SRR8051550":"Infloresence-Rep1",
"SRR8051547":"Infloresence-Rep2",
"SRR8051546":"Infloresence-Rep3",
"SRR8051557":"Leaf-Rep1",
"SRR8051556":"Leaf-Rep2",
"SRR8051555":"Leaf-Rep3",
"SRR8051560":"Root-Rep1",
"SRR8051559":"Root-Rep2",
"SRR8051558":"Root-Rep3",
"SRR7974062":"Seed-Rep1",
"SRR7974059":"Seed-Rep2",
"SRR8051561":"Seed-Rep3",
"SRR8051549":"Seedling-Rep1",
"SRR8051548":"Seedling-Rep2",
"SRR8051551":"Seedling-Rep3",
"SRR8051554":"Stem-Rep1",
"SRR8051553":"Stem-Rep2",
"SRR8051552":"Stem-Rep3",
"SRR9002107":"Embryo-Rep1",
"SRR9002108":"Embryo-Rep2",
"SRR9106180":"Germinating Seed-0h Imbibation-Rep1",
"SRR9106181":"Germinating Seed-0h Imbibation-Rep2",
"SRR9106182":"Germinating Seed-0h Imbibation-Rep3",
"SRR9117223":"Germinating Seed-0h Imbibation-Rep1",
"SRR9117222":"Germinating Seed-0h Imbibation-Rep2",
"SRR9117221":"Germinating Seed-0h Imbibation-Rep3",
"SRR9117218":"Germinating Seed-8h Imbibition-Rep1",
"SRR9117220":"Germinating Seed-8h Imbibition-Rep2",
"SRR9117219":"Germinating Seed-8h Imbibition-Rep3",
"SRR10991573":"Seedling-Control-Rep1",
"SRR10991574":"Seedling-Control-Rep2",
"SRR10991575":"Seedling-Control-Rep3",
"SRR10991576":"Seedling-Heat Stress-Rep1",
"SRR10991577":"Seedling-Heat Stress-Rep2",
"SRR10991578":"Seedling-Heat Stress-Rep3",
"SRR11784188":"Anther",
"SRR11784189":"Glume",
"SRR13808035":"Leaf-Booting stage-Rep2",
"SRR13808036":"Leaf-Booting stage-Rep1",
"SRR13808042":"Spikelet-Booting stage-Rep2",
"SRR13808043":"Spikelet-Booting stage-Rep1",
"SRR13808034":"Leaf-V2",
"SRR13808040":"Leaf Sheath-V2-Rep2",
"SRR13808041":"Leaf Sheath-V2-Rep1",
"SRR13808037":"Root-V2",
"SRR13808039":"Seedling-V2",
"SRR12545732":"Leaf",
"SRR12545731":"Panicle",
"SRR12545730":"Root",
"SRR16526866":"Leaf-flag-rep1",
"SRR16526865":"Leaf-flag-rep2",
"SRR16526864":"Leaf-flag-rep3",
"SRR16526863":"Panicles, young-Rep1",
"SRR16526872":"Panicles, young-Rep2",
"SRR16526871":"Panicles, young-Rep3",
"SRR17151219":"Endosperm",
"SRR18827799":"Leaf-seedling, Cold Stress-24 hr-Rep1",
"SRR18827798":"Leaf-seedling, Cold Stress-24 hr-Rep2",
"SRR18827800":"Leaf-seedling, Cold Stress-24 hr-Rep3",
"SRR18827805":"Leaf-seedling, Control-24 hr-Rep1",
"SRR18827806":"Leaf-seedling, Control-24 hr-Rep2",
"SRR18827797":"Leaf-seedling, Control-24 hr-Rep3",
"SRR19614426":"Seed-22 DAP-Rep3",
"SRR19614428":"Seed-22 DAP-Rep1",
"SRR19614427":"Seed-22 DAP-Rep2",
"SRR25884998":"Seedling-Control-Rep1",
"SRR25884997":"Seedling-Control-Rep2",
"SRR25884988":"Seedling-Control-Rep3",
"SRR25884981":"Seedling-Cold-Rep1",
"SRR25884996":"Seedling-Cold-Rep2",
"SRR25884995":"Seedling-Cold-Rep3",
"SRR26772170":"Seedling-Salt stress-12 hr-Rep3",
"SRR26772168":"Seedling-Salt stress-24 hr-Rep2",
"SRR26772159":"Seedling-Salt stress-0 hr-Rep1",
"SRR26772174":"Seedling-Salt stress-0 hr-Rep2",
"SRR26772173":"Seedling-Salt stress-0 hr-Rep3",
"SRR26772172":"Seedling-Salt stress-12 hr-Rep1",
"SRR26772171":"Seedling-Salt stress-12 hr-Rep2",
"SRR26772169":"Seedling-Salt stress-24 hr-Rep1",
"SRR26772167":"Seedling-Salt stress-24 hr-Rep3",
}

TREATMENT_KEYWORDS =(
"control",
"stress",
"drought",
"hypoxia",
"cold",
"salt",
"heat",
"imbib",
"0 hr",
"8h",
"24 hr",
)


class ExpressionPatternService :
    def __init__ (self ,matrix_path :str |None =None ):
        env_path =os .getenv ("TRANSCRIPT_EXPRESSION_XLSX")
        self .matrix_path =matrix_path or env_path or self ._discover_matrix_path ()
        self .matrix :pd .DataFrame |None =None 

    def collect (self ,gene_info :dict [str ,Any ])->dict [str ,Any ]:
        transcript_ids =[tid .strip ()for tid in gene_info .get ("Transcript_IDs",[])if str (tid ).strip ()]
        if not transcript_ids and gene_info .get ("Transcript_ID"):
            transcript_ids =[gene_info ["Transcript_ID"].strip ()]
        expression_model_ids =[mid .strip ()for mid in gene_info .get ("MSU_IDs",[])if str (mid ).strip ()]
        result ={
        "gene_id":gene_info .get ("RAP_ID",""),
        "transcript_id":transcript_ids [0 ]if transcript_ids else "",
        "transcript_ids":transcript_ids ,
        "transcript_count":len (transcript_ids ),
        "expression_model_ids":expression_model_ids ,
        "expression_model_count":len (expression_model_ids ),
        "matrix_path":self .matrix_path ,
        "available":False ,
        "transcript_results":[],
        "warnings":[],
        "summary_text":"",
        }

        if not transcript_ids :
            result ["warnings"].append ("This gene's Transcript ID was not found in the genome annotation.")
        if not expression_model_ids :
            result ["warnings"].append ("The corresponding MSU Expression model ID for this gene was not found; unable to search in the expression matrix.")
            return result 

        path =Path (self .matrix_path )if self .matrix_path else None 
        if not path or not path .exists ():
            result ["warnings"].append (f"Not FoundExpressionMatrixFile: {self .matrix_path or 'No .xlsx file automatically detected'}")
            return result 

        try :
            matrix =self ._load_matrix ()
        except Exception as exc :
            result ["warnings"].append (f"readExpressionMatrixFailed: {exc }")
            return result 

        if len (transcript_ids )!=len (expression_model_ids ):
            result ["warnings"].append (
            f"GenegroupAnnotationTranscriptnumber ({len (transcript_ids )}) withExpressionmatrix model number ({len (expression_model_ids )}) inconsistent，will beExpressionmatrix model mapping one by one。"
            )

        for index ,expression_model_id in enumerate (expression_model_ids ):
            transcript_label =transcript_ids [index ]if index <len (transcript_ids )else expression_model_id 
            transcript_result =self ._collect_one_transcript (matrix ,transcript_label ,expression_model_id )
            if transcript_result ["available"]:
                result ["transcript_results"].append (transcript_result )
            else :
                result ["warnings"].extend (transcript_result ["warnings"])

        result ["available"]=bool (result ["transcript_results"])
        if not result ["transcript_results"]and not result ["warnings"]:
            result ["warnings"].append ("No transcripts were found in the expression matrix that can be used for plotting.")
        result ["summary_text"]=self ._build_summary_text (result )
        return result 

    def _load_matrix (self )->pd .DataFrame :
        if self .matrix is None :
            self .matrix =pd .read_excel (self .matrix_path )
            first_col =self .matrix .columns [0 ]
            self .matrix =self .matrix .rename (columns ={first_col :"transcript_id"})
            self .matrix ["transcript_id"]=self .matrix ["transcript_id"].astype (str ).str .strip ()
            self .matrix =self .matrix .set_index ("transcript_id")
        return self .matrix 

    def _collect_one_transcript (self ,matrix :pd .DataFrame ,transcript_id :str ,expression_model_id :str )->dict [str ,Any ]:
        result ={
        "transcript_id":transcript_id ,
        "expression_model_id":expression_model_id ,
        "available":False ,
        "developmental_plot_png":None ,
        "treatment_plot_png":None ,
        "developmental_summary":[],
        "treatment_summary":[],
        "warnings":[],
        }
        if expression_model_id not in matrix .index :
            result ["warnings"].append (f"Expressionin the matrixNot Foundmodel: {expression_model_id }（correspondAnnotationTranscript {transcript_id }）")
            return result 

        expr =matrix .loc [expression_model_id ]
        sample_df =self ._build_sample_frame (expr )
        if sample_df .empty :
            result ["warnings"].append (f"{expression_model_id } notHasandSampleAnnotationtable showing overlap of SRA Sample。")
            return result 

        dev_df =self ._aggregate (sample_df [sample_df ["sample_type"]=="developmental"])
        treatment_df =self ._aggregate (sample_df [sample_df ["sample_type"]=="treatment"])
        result ["available"]=True 
        result ["developmental_summary"]=dev_df .to_dict ("records")if not dev_df .empty else []
        result ["treatment_summary"]=treatment_df .to_dict ("records")if not treatment_df .empty else []
        if not dev_df .empty :
            result ["developmental_plot_png"]=self ._plot_grouped_bars (
            dev_df ,
            title =f"Expression Pattern by Organ/Stage: {transcript_id } | {expression_model_id }",
            color ="#4e79a7",
            )
        else :
            result ["warnings"].append (f"{expression_model_id } organ not identified/time-relatedSample。")
        if not treatment_df .empty :
            result ["treatment_plot_png"]=self ._plot_grouped_bars (
            treatment_df ,
            title =f"Expression Pattern by Treatment: {transcript_id } | {expression_model_id }",
            color ="#e15759",
            )
        else :
            result ["warnings"].append (f"{expression_model_id } not identifiedTreatmentclassSample。")
        return result 

    def _build_sample_frame (self ,expr :pd .Series )->pd .DataFrame :
        rows =[]
        for sample_id ,value in expr .items ():
            sample_key =str (sample_id ).strip ()
            description =SAMPLE_DESCRIPTION_MAP .get (sample_key )
            if not description :
                continue 
            rows .append (
            {
            "sample_id":sample_key ,
            "description":description ,
            "group_label":self ._normalize_description (description ),
            "sample_type":self ._sample_type (description ),
            "expression":self ._safe_float (value ),
            }
            )
        return pd .DataFrame (rows )

    def _aggregate (self ,df :pd .DataFrame )->pd .DataFrame :
        if df .empty :
            return pd .DataFrame (columns =["group_label","mean_tpm","replicate_count"])
        grouped =(
        df .groupby ("group_label",as_index =False )
        .agg (mean_tpm =("expression","mean"),replicate_count =("expression","size"))
        .sort_values ("mean_tpm",ascending =False )
        )
        grouped ["mean_tpm"]=grouped ["mean_tpm"].round (3 )
        return grouped 

    def _plot_grouped_bars (self ,df :pd .DataFrame ,title :str ,color :str )->bytes :
        fig_height =max (4.5 ,len (df )*0.45 )
        fig ,ax =plt .subplots (figsize =(11.5 ,fig_height ))
        bars =ax .barh (df ["group_label"],df ["mean_tpm"],color =color ,edgecolor ="black",alpha =0.85 )
        ax .set_title (title ,fontsize =14 ,fontweight ="bold",pad =12 )
        ax .set_xlabel ("Mean TPM",fontsize =11 ,fontweight ="bold")
        ax .set_ylabel ("Condition",fontsize =11 ,fontweight ="bold")
        ax .grid (axis ="x",linestyle ="--",alpha =0.35 )
        ax .invert_yaxis ()

        max_value =float (df ["mean_tpm"].max ())if not df .empty else 0 
        offset =max (max_value *0.015 ,0.05 )
        for bar ,value in zip (bars ,df ["mean_tpm"]):
            ax .text (value +offset ,bar .get_y ()+bar .get_height ()/2 ,f"{value :.2f}",va ="center",fontsize =9 )

        plt .tight_layout ()
        output =BytesIO ()
        fig .savefig (output ,format ="png",dpi =220 ,bbox_inches ="tight")
        plt .close (fig )
        output .seek (0 )
        return output .getvalue ()

    @staticmethod 
    def _normalize_description (description :str )->str :
        normalized =description .strip ()
        for token in ("-Rep1","-Rep2","-Rep3","-rep1","-rep2","-rep3"):
            if normalized .endswith (token ):
                normalized =normalized [:-len (token )]
        return normalized .strip ()

    @staticmethod 
    def _sample_type (description :str )->str :
        text =description .lower ()
        if any (keyword in text for keyword in TREATMENT_KEYWORDS ):
            return "treatment"
        return "developmental"

    @staticmethod 
    def _safe_float (value :Any )->float :
        number =pd .to_numeric (value ,errors ="coerce")
        if pd .isna (number ):
            return 0.0 
        return float (number )

    @staticmethod 
    def _discover_matrix_path ()->str :
        exact_candidates =sorted (EXPRESSION_DIR .glob ("osa1_r7.all_models.gene_exp_matrix*.xlsx"))
        if exact_candidates :
            return str (exact_candidates [0 ])
        generic_candidates =sorted (EXPRESSION_DIR .glob ("*.xlsx"))
        if generic_candidates :
            return str (generic_candidates [0 ])
        return str (EXPRESSION_DIR /"osa1_r7.all_models.gene_exp_matrix (2).xlsx")

    @staticmethod 
    def _build_summary_text (result :dict [str ,Any ])->str :
        lines =["[Expression Pattern Evidence]"]
        lines .append (f"- GenegroupAnnotationTranscriptnumber: {result .get ('transcript_count',0 )}")
        lines .append (f"- Expressionmatrix model number: {result .get ('expression_model_count',0 )}")
        for item in result .get ("transcript_results",[]):
            dev =item .get ("developmental_summary",[])
            treat =item .get ("treatment_summary",[])
            top_dev =dev [0 ]if dev else None 
            top_treat =treat [0 ]if treat else None 
            lines .append (f"- Transcript {item ['transcript_id']} / model {item .get ('expression_model_id','N/A')}:")
            if top_dev :
                lines .append (
                f"  * organ/highest periodExpression: {top_dev ['group_label']} (mean TPM={top_dev ['mean_tpm']})"
                )
            if top_treat :
                lines .append (
                f"  * TreatmentConditionhighestExpression: {top_treat ['group_label']} (mean TPM={top_treat ['mean_tpm']})"
                )
        for warning in result .get ("warnings",[]):
            lines .append (f"- warning: {warning }")
        if not result .get ("transcript_results"):
            lines .append ("- No usable expression pattern plot has been generated currently.")
        return "\n".join (lines )
