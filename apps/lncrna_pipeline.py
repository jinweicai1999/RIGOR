#!/usr/bin/env python3
"""
lncRNA trait-analysis entrypoint.

Novel lncRNA:
- local co-expression matrix
- cis/trans neighborhood inference
- supplementary cis/trans and direct GWAS evidence
- literature-backed functional inference
- proxy GWAS scan on protein partners

Known lncRNA:
- location lookup only
"""

import argparse 
import json 
import sys 
from pathlib import Path 

from rich .console import Console 
from rich .panel import Panel 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .agents .gwas_agent import GwasAgent 
from src .agents .literature_agent import LiteratureAgent 
from src .agents .lncrna_agent import LncRnaAgent 
from src .project_paths import ANNOTATION_DIR ,NCRNA_DIR 
from src .utils .gene_annotator import RiceGeneAnnotator 
from src .utils .known_lncrna_annotator import KnownLncRnaAnnotator 
from src .utils .lncrna_annotator import LncRnaAnnotator 
from src .utils .lncrna_supplement_service import LncRnaSupplementService 

STEP_BEGIN ="PROTEIN_STEP_JSON_BEGIN"
STEP_END ="PROTEIN_STEP_JSON_END"


def emit_step (step ,title ,data ):
    payload ={"step":step ,"title":title ,"status":"done","data":data }
    print (STEP_BEGIN )
    print (json .dumps (payload ,ensure_ascii =False ))
    print (STEP_END )
    sys .stdout .flush ()


class LncRnaTraitAgent :
    def __init__ (self ):
        self .console =Console ()

        self .gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")
        self .map_file =str (ANNOTATION_DIR /"RAP-MSU.txt")
        self .annotator =RiceGeneAnnotator (self .gff_file ,self .map_file )
        self .annotator .build_database ()

        self .lnc_annotator =LncRnaAnnotator (str (NCRNA_DIR /"novel_lncRNA.gtf"))
        self .lnc_annotator .build_database ()

        self .known_lnc_annotator =KnownLncRnaAnnotator ()
        self .known_lnc_annotator .build_database ()

        self .lnc_agent =LncRnaAgent ()
        self .lit_agent =LiteratureAgent ()
        self .gwas_agent =GwasAgent ()
        self .supplement_service =LncRnaSupplementService ()

    def run (self ,lncrna_id ):
        if "lncn"in lncrna_id .lower ():
            self .run_novel_lncrna (lncrna_id )
            return 
        if self .known_lnc_annotator .has_lncrna (lncrna_id ):
            self .run_known_lncrna (lncrna_id )
            return 
        self .console .print (f"[red]❌ Nomethod identify lncRNA: {lncrna_id }[/red]")
        self .console .print ("[yellow]Tip: Novel lncRNA uses OSA-LNCN*, known lncRNA uses OSA-LNC*[/yellow]")
        sys .exit (1 )

    def run_by_coordinate (self ,chrom ,start ,end ,min_overlap =0.8 ):
        self .console .print (
        Panel (
        f"[bold cyan]lncRNA coordinatesSearchStart[/bold cyan]\n"
        f"query interval: {chrom }:{start }-{end }\n"
        f"minimum similarity threshold: {min_overlap *100 :.0f}%",
        border_style ="cyan",
        )
        )

        hits =self .lnc_annotator .query_by_coordinate (chrom ,start ,end ,min_overlap =min_overlap )
        if not hits :
            self .console .print (
            f"[yellow]⚠ Not Foundwith {chrom }:{start }-{end } similarity ≥ {min_overlap *100 :.0f}% 's novel lncRNA。[/yellow]"
            )
            sys .exit (0 )

        best =hits [0 ]
        self .console .print (
        f"[green]✓ best hit: {best ['id']} {best ['chromosome']}:{best ['start']}-{best ['end']} "
        f"({best ['strand']}), similarity {best ['similarity']*100 :.1f}%[/green]"
        )
        self .run_novel_lncrna (best ["id"],coord_hit =best )

    def run_novel_lncrna (self ,lncrna_id ,coord_hit =None ):
        self .console .print (
        Panel (
        f"[bold cyan]lncRNA IntelligentAnalysisdeviceStart[/bold cyan]\n"
        f"target: {lncrna_id }\n"
        f"type: Novel lncRNA\n"
        f"pattern: totalExpressionNetwork + supplementary table integrationAnalysis",
        border_style ="cyan",
        )
        )

        self .console .print ("\n[Phase 1] Analysis of the local co-expression network and cis/trans relationship...")
        network_genes ,target_ids =self ._collect_local_network_genes (lncrna_id )
        if not network_genes :
            self .console .print ("[red]❌ Not Found: Effective co-expression protein gene[/red]")
            sys .exit (1 )

        lnc_location =self .lnc_annotator .query (lncrna_id )
        nearby_genes =self .lnc_annotator .find_nearby_genes (lncrna_id ,self .annotator ,distance =10000 )if lnc_location else []
        nearby_gene_ids =[item ["Gene_ID"]for item in nearby_genes ]
        cis_count =sum (1 for item in network_genes if item ["RAP_ID"]in nearby_gene_ids )
        trans_count =len (network_genes )-cis_count 

        local_coexp_list =[]
        for gene in network_genes :
            rap =gene .get ("RAP_ID","")
            local_coexp_list .append (
            {
            "rap_id":rap ,
            "symbol":gene .get ("Primary_Symbol",""),
            "correlation":round (float (gene .get ("Correlation",0 )),4 ),
            "regulation":"Cis"if rap in nearby_gene_ids else "Trans",
            "description":gene .get ("Description",""),
            }
            )

        supplement_coexp =self .supplement_service .get_coexpression_evidence (lncrna_id )
        supplement_gwas =self .supplement_service .get_gwas_evidence (lncrna_id )
        supplement_report =self ._build_supplement_report (lncrna_id ,supplement_coexp ,supplement_gwas )

        emit_step (
        "lncrna_basic",
        "Basic Information on lncRNA",
        {
        "lncrna_id":lncrna_id ,
        "chromosome":lnc_location ["Chromosome"]if lnc_location else "",
        "start":lnc_location ["Start"]if lnc_location else 0 ,
        "end":lnc_location ["End"]if lnc_location else 0 ,
        "strand":lnc_location ["Strand"]if lnc_location else ".",
        "length":lnc_location ["Length"]if lnc_location else 0 ,
        "coord_hit":coord_hit ,
        "coexp_genes":local_coexp_list ,
        "cis_count":cis_count ,
        "trans_count":trans_count ,
        "supplement_mapping_ids":supplement_coexp .get ("summary",{}).get ("mapped_supplement_ids",[]),
        "supplement_pair_count":supplement_coexp .get ("summary",{}).get ("pair_count",0 ),
        "supplement_direct_trait_count":supplement_gwas .get ("summary",{}).get ("direct_trait_count",0 ),
        },
        )

        if supplement_report :
            emit_step (
            "lncrna_supplement",
            "Supplement: Additional Evidence",
            {
            "report":supplement_report ,
            "coexpression_pairs":supplement_coexp .get ("pairs",[]),
            "trait_summary":supplement_gwas .get ("trait_summary",[]),
            "pve_hits":supplement_gwas .get ("pve_hits",[]),
            },
            )

        self .console .print ("\n[Phase 2] Combine local and supplementary tables partner, FetchLiterature and perform FunctionInference...")
        literature_targets =self ._build_literature_targets (network_genes ,supplement_coexp )
        abstracts =self .lit_agent .fetch_abstracts_route_b (literature_targets [:8 ])
        mechanism_report =""
        refs =[]
        if abstracts .strip ():
            extra_context =self ._build_literature_extra_context (
            lncrna_id =lncrna_id ,
            local_network_genes =network_genes ,
            nearby_gene_ids =nearby_gene_ids ,
            supplement_coexp =supplement_coexp ,
            supplement_gwas =supplement_gwas ,
            )
            mechanism_report =self .lit_agent .summarize_ncrna (
            lncrna_id ,
            literature_targets [:8 ],
            abstracts ,
            rna_type ="lncRNA",
            extra_context =extra_context ,
            )
            refs =self ._parse_references (abstracts )

        emit_step (
        "lncrna_coexpression",
        "ExpressionFunctionInference",
        {
        "report":mechanism_report ,
        "references":refs ,
        },
        )

        self .console .print ("\n[Phase 3] Evaluation proxy GWAS and direct lncRNA GWAS evidence...")
        gwas_scan_ids =self ._build_gwas_scan_targets (target_ids ,supplement_coexp )
        gwas_df =self .gwas_agent .scan_network_gwas (gwas_scan_ids ,self .annotator )if gwas_scan_ids else None 
        proxy_gwas_report =""
        if gwas_df is not None and not gwas_df .empty :
            proxy_gwas_report =self .gwas_agent .summarize_ncrna (lncrna_id ,gwas_scan_ids ,gwas_df )

        integrated_gwas_report =self ._build_integrated_gwas_report (
        lncrna_id =lncrna_id ,
        supplement_gwas =supplement_gwas ,
        proxy_gwas_report =proxy_gwas_report ,
        proxy_gwas_df =gwas_df ,
        gwas_scan_ids =gwas_scan_ids ,
        )

        emit_step (
        "lncrna_gwas",
        "GWAS breeding effect evaluation",
        {
        "report":integrated_gwas_report ,
        "has_signal":bool (
        (gwas_df is not None and not gwas_df .empty )or supplement_gwas .get ("available",False )
        ),
        },
        )

        self .console .print ("\n[bold green]✅ lncRNA AnalysisCompleted。[/bold green]")

    def run_known_lncrna (self ,lncrna_id ):
        self .console .print (
        Panel (
        f"[bold cyan]lncRNA IntelligentAnalysisdeviceStart[/bold cyan]\n"
        f"target: {lncrna_id }\n"
        f"type: Known lncRNA\n"
        f"pattern: position mappingAnalysis",
        border_style ="cyan",
        )
        )

        lnc_info =self .known_lnc_annotator .query (lncrna_id )
        if not lnc_info :
            self .console .print (f"[red]❌ Not Found {lncrna_id } 'sInformation[/red]")
            sys .exit (1 )

        loc_id =lnc_info ["LOC_ID"]
        loc_info =self .annotator .query (loc_id )
        if not loc_info or "Location"not in loc_info :
            self .console .print (f"[yellow]⚠ Nomethod to obtain {loc_id } 'sGenegroup position[/yellow]")
            sys .exit (1 )

        emit_step (
        "lncrna_basic",
        "Basic Information on lncRNA",
        {
        "lncrna_id":lncrna_id ,
        "chromosome":"",
        "start":0 ,
        "end":0 ,
        "strand":".",
        "length":0 ,
        "coord_hit":None ,
        "coexp_genes":[],
        "cis_count":0 ,
        "trans_count":0 ,
        },
        )
        emit_step (
        "lncrna_coexpression",
        "ExpressionFunctionInference",
        {
        "report":(
        f"Hasknown lncRNA `{lncrna_id }` Currently only supports positional resolution。\n\n"
        f"- RefSeq/LOC map: `{loc_id }`\n"
        f"- Annotationposition: `{loc_info ['Location']}`\n"
        f"- limitation: currently lacks dedicatedExpressionMatrix，NoImaging method novel lncRNA like performing full coExpressionInference。"
        ),
        "references":[],
        },
        )
        emit_step (
        "lncrna_gwas",
        "GWAS breeding effect evaluation",
        {
        "report":"The lncRNA branch currently lacks direct GWAS or ExpressionNetwork supplementary evidence.",
        "has_signal":False ,
        },
        )

    def _collect_local_network_genes (self ,lncrna_id ):
        network_genes =[]
        target_ids =[]
        search_n =5 

        while len (network_genes )<5 and search_n <=50 :
            top_genes ,error =self .lnc_agent .get_top_coexpressed_genes (lncrna_id ,top_n =search_n )
            if error :
                self .console .print (f"[yellow]⚠ {error }[/yellow]")
                return [],[]

            for _ ,row in top_genes .iterrows ():
                gene_id =row ["Gene_ID"]
                if gene_id in target_ids :
                    continue 
                gene_info =self .annotator .query (gene_id )
                if not gene_info :
                    continue 
                gene_info ["Correlation"]=float (row ["Correlation"])
                network_genes .append (gene_info )
                target_ids .append (gene_id )
                if len (network_genes )>=5 :
                    break 
            search_n +=10 

        return network_genes ,target_ids 

    def _build_literature_targets (self ,local_network_genes ,supplement_coexp ):
        literature_targets =[]
        seen =set ()

        for gene in local_network_genes :
            rap_id =gene .get ("RAP_ID","")
            if rap_id and rap_id not in seen :
                literature_targets .append (gene )
                seen .add (rap_id )

        for gene_id in supplement_coexp .get ("partner_gene_ids",[])[:8 ]:
            clean_id =str (gene_id ).strip ()
            if not clean_id or clean_id in seen :
                continue 
            gene_info =self .annotator .query (clean_id )
            if not gene_info :
                continue 
            gene_info ["Correlation"]=0.0 
            literature_targets .append (gene_info )
            seen .add (clean_id )

        return literature_targets 

    def _build_gwas_scan_targets (self ,local_target_ids ,supplement_coexp ):
        scan_ids =[]
        seen =set ()
        for gene_id in list (local_target_ids )+supplement_coexp .get ("partner_gene_ids",[])[:10 ]:
            clean_id =str (gene_id ).strip ()
            if not clean_id or clean_id in seen :
                continue 
            seen .add (clean_id )
            scan_ids .append (clean_id )
        return scan_ids 

    def _build_literature_extra_context (
    self ,
    lncrna_id ,
    local_network_genes ,
    nearby_gene_ids ,
    supplement_coexp ,
    supplement_gwas ,
    ):
        lines =[]
        local_bits =[]
        for gene in local_network_genes [:5 ]:
            rap_id =gene .get ("RAP_ID","")
            corr =float (gene .get ("Correlation",0 ))
            local_bits .append (
            f"{rap_id } (r={corr :.3f}, {'cis candidates'if rap_id in nearby_gene_ids else 'trans candidates'})"
            )
        if local_bits :
            lines .append ("Local Expression matrix top total ExpressionGene:"+"; ".join (local_bits ))

        supp_summary =supplement_coexp .get ("summary",{})
        if supplement_coexp .get ("available"):
            lines .append (
            f"Supplement Sheet5 Hit {supp_summary .get ('pair_count',0 )} line lncRNA-mRNA pairing，"
            f"involve {supp_summary .get ('unique_partner_count',0 )} a unique mRNA；"
            f"cis={supp_summary .get ('cis_count',0 )}, trans={supp_summary .get ('trans_count',0 )}, "
            f"positive correlation={supp_summary .get ('positive_pair_count',0 )}, negatively correlated={supp_summary .get ('negative_pair_count',0 )}。"
            )
            pair_bits =[]
            for item in supplement_coexp .get ("pairs",[])[:6 ]:
                pair_bits .append (
                f"{item ['partner_gene']} ({item ['regulation']}, r={item ['correlation']}, tissue={item ['tissue']}, state={item ['state_class']})"
                )
            if pair_bits :
                lines .append ("Representative pairings for supplements:"+"; ".join (pair_bits ))

        if supplement_gwas .get ("available"):
            gwas_summary =supplement_gwas .get ("summary",{})
            lines .append (
            f"Supplement Sheet9/10 for {lncrna_id } providing direct GWAS Evidence："
            f"{gwas_summary .get ('direct_snp_hit_count',0 )} each SNP Hit，"
            f"{gwas_summary .get ('direct_trait_count',0 )} trait，"
            f"best P={gwas_summary .get ('best_p_value')}, maximum PVE={gwas_summary .get ('best_pve')}。"
            )
        return "\n".join (lines )if lines else "No additional contextual evidence."

    def _build_supplement_report (self ,lncrna_id ,supplement_coexp ,supplement_gwas ):
        lines =[]

        if supplement_coexp .get ("available"):
            summary =supplement_coexp .get ("summary",{})
            lines .append (f"## {lncrna_id } 's Supplement totalExpressionadditional evidence")
            lines .append (
            f"- map to {len (summary .get ('mapped_supplement_ids',[]))} supplementary table lncRNA entry"
            )
            lines .append (
            f"- Hit {summary .get ('pair_count',0 )} line lncRNA-mRNA pairing，involve "
            f"{summary .get ('unique_partner_count',0 )} a unique mRNA"
            )
            lines .append (
            f"- cis={summary .get ('cis_count',0 )}，trans={summary .get ('trans_count',0 )}，"
            f"positive correlation={summary .get ('positive_pair_count',0 )}，negatively correlated={summary .get ('negative_pair_count',0 )}"
            )
            if summary .get ("tissues"):
                lines .append (f"- PrimaryTissue: {', '.join (summary ['tissues'])}")
            if summary .get ("state_classes"):
                lines .append (f"- Primarystate: {', '.join (summary ['state_classes'])}")
            lines .append ("")
            lines .append ("### Representative lncRNA-mRNA pairs")
            for item in supplement_coexp .get ("pairs",[])[:8 ]:
                lines .append (
                f"- {item ['supplement_id']} -> {item ['partner_gene']} | {item ['regulation']} | "
                f"r={item ['correlation']} | tissue={item ['tissue']} | state={item ['state_class']}"
                )

        if supplement_gwas .get ("available"):
            summary =supplement_gwas .get ("summary",{})
            if lines :
                lines .append ("")
            lines .append (f"## {lncrna_id } 's Supplement directly GWAS Evidence")
            lines .append (
            f"- directly SNP hit count: {summary .get ('direct_snp_hit_count',0 )}"
            )
            lines .append (
            f"- coverage trait number: {summary .get ('direct_trait_count',0 )}"
            )
            if summary .get ("best_p_value")is not None :
                lines .append (f"- strongest GWAS P value: {summary ['best_p_value']:.3e}")
            if summary .get ("best_pve")is not None :
                lines .append (f"- maximum PVE: {summary ['best_pve']:.2f}%")
            if supplement_gwas .get ("trait_summary"):
                lines .append ("")
                lines .append ("### Primary traits")
                for item in supplement_gwas ["trait_summary"][:8 ]:
                    lines .append (
                    f"- {item ['trait']}: best P={item ['min_p_value']:.3e}, SNP number={item ['snp_count']}"
                    )
            if supplement_gwas .get ("pve_hits"):
                lines .append ("")
                lines .append ("### High PVE Records")
                for item in supplement_gwas ["pve_hits"][:5 ]:
                    lines .append (
                    f"- {item ['trait']}: {item ['snp']} (PVE={item ['pve']}%)"
                    )

        return "\n".join (lines ).strip ()

    def _build_integrated_gwas_report (
    self ,
    lncrna_id ,
    supplement_gwas ,
    proxy_gwas_report ,
    proxy_gwas_df ,
    gwas_scan_ids ,
    ):
        sections =[f"# {lncrna_id } GWAS comprehensiveEvaluation"]

        if supplement_gwas .get ("available"):
            sections .append ("## 1. Direct GWAS evidence for lncRNA")
            summary =supplement_gwas .get ("summary",{})
            sections .append (
            f"- Supplement Hit {summary .get ('direct_snp_hit_count',0 )} each lncRNA direct association SNP，"
            f"cover {summary .get ('direct_trait_count',0 )} trait。"
            )
            if summary .get ("best_p_value")is not None :
                sections .append (f"- strongestSignal P value: {summary ['best_p_value']:.3e}")
            if summary .get ("best_pve")is not None :
                sections .append (f"- maximum PVE: {summary ['best_pve']:.2f}%")
            for item in supplement_gwas .get ("trait_summary",[])[:8 ]:
                sections .append (
                f"- {item ['trait']}: best P={item ['min_p_value']:.3e}, SNP number={item ['snp_count']}"
                )
        else :
            sections .append ("## 1. Direct GWAS evidence for lncRNA")
            sections .append ("- No direct SNP evidence mapping to this novel lncRNA was found in the current Supplement GWAS table.")

        sections .append ("")
        sections .append ("## 2. Proxy GWAS evidence for shared expression partners")
        if proxy_gwas_df is not None and not proxy_gwas_df .empty :
            sections .append (
            f"- use {len (gwas_scan_ids )} each partner Geneconduct proxyScan，capture {len (proxy_gwas_df )} significant GWAS record。"
            )
            sections .append ("")
            sections .append (proxy_gwas_report .strip ()or "- Proxy GWAS signal detected, but the model summary is empty.")
        else :
            sections .append (
            f"- HasScan {len (gwas_scan_ids )} each partner Gene，but no significant was found proxy GWAS Signal，or local GWAS FileInadequate。"
            )

        sections .append ("")
        sections .append ("## 3. Comprehensive Explanation")
        if supplement_gwas .get ("available")and proxy_gwas_df is not None and not proxy_gwas_df .empty :
            sections .append (
            "- This lncRNA has both direct lncRNA site association and partner mRNA proxy association evidence, "
            "making it a candidate closer to real trait-regulation loci than modules inferred only from co-expression partners."
            )
        elif supplement_gwas .get ("available"):
            sections .append (
            "- The current stronger evidence comes from direct GWAS hits of lncRNA itself, indicating that this module should not rely solely on proxy explanations from partner mRNA."
            )
        elif proxy_gwas_df is not None and not proxy_gwas_df .empty :
            sections .append (
            "- Current GWAS explanations still primarily rely on partner mRNA proxy signals, which are suitable as functional clues but should not replace direct lncRNA causal evidence."
            )
        else :
            sections .append (
            "- Currently, there is a lack of direct or proxy GWAS strong evidence, and the trait significance of lncRNA still needs to rely on expression, co-expression, and literature-level support."
            )

        return "\n".join (sections ).strip ()

    def _parse_references (self ,abstracts ):
        refs =[]
        import re 

        for block in re .split (r"\n(?=- PMID:)",abstracts .strip ()):
            pmid_m =re .search (r"PMID:\s*(\d+)",block )
            title_m =re .search (r"Title:\s*(.+)",block )
            journal_m =re .search (r"Journal:\s*(.+)",block )
            if not pmid_m :
                continue 
            raw_journal =journal_m .group (1 ).strip ()if journal_m else ""
            clean_journal =re .sub (r"\s*\(DOI:[^)]*\)\s*$","",raw_journal ).strip ()
            refs .append (
            {
            "pmid":pmid_m .group (1 ),
            "title":title_m .group (1 ).strip ()if title_m else "",
            "journal":clean_journal ,
            }
            )
        return refs 


if __name__ =="__main__":
    parser =argparse .ArgumentParser (description ="Rice lncRNA Trait Gene Analysis System")
    group =parser .add_mutually_exclusive_group (required =True )
    group .add_argument ("--gene",help ="lncRNA ID (OSA-LNCN* or OSA-LNC*)")
    group .add_argument ("--coord",help ="Genomic coordinates, format: CHR1,12321312-12321999")
    parser .add_argument ("--min-overlap",type =float ,default =0.8 ,help ="Minimum overlap ratio for coordinate matching")
    args =parser .parse_args ()

    agent =LncRnaTraitAgent ()

    if args .coord :
        coord_str =args .coord .strip ().replace (":",",").replace (" ",",")
        try :
            parts =[part .strip ()for part in coord_str .split (",")if part .strip ()]
            if len (parts )==2 :
                chrom =parts [0 ]
                start_s ,end_s =parts [1 ].split ("-")
            elif len (parts )==3 :
                chrom ,start_s ,end_s =parts 
            else :
                raise ValueError ("bad coordinate format")
            start =int (start_s .replace (",",""))
            end =int (end_s .replace (",",""))
        except Exception :
            print (f"❌ Nomethod to resolve coordinates: {args .coord }")
            print ("Supported format: CHR1,12321312-12321999 | chr1:12321312-12321999")
            sys .exit (1 )
        agent .run_by_coordinate (chrom ,start ,end ,min_overlap =args .min_overlap )
    else :
        agent .run (args .gene )
