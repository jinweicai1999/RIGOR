import argparse 
import json 
import os 
import sys 
import subprocess 
from pathlib import Path 
from rich .console import Console 
from rich .panel import Panel 
from rich .markdown import Markdown 

STEP_BEGIN ="PROTEIN_STEP_JSON_BEGIN"
STEP_END ="PROTEIN_STEP_JSON_END"


def _emit_step (step ,title_zh ,title_en ,data ):
    lang ="en"if os .getenv ("APP_LANG","zh").lower ().startswith ("en")else "zh"
    payload ={
    "step":step ,
    "title":title_en if lang =="en"else title_zh ,
    "status":"done",
    "data":data ,
    }
    print (STEP_BEGIN )
    print (json .dumps (payload ,ensure_ascii =False ))
    print (STEP_END )
    sys .stdout .flush ()

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .agents .ncrna_agent import NcRnaAgent 
from src .utils .gene_annotator import RiceGeneAnnotator 
from src .agents .literature_agent import LiteratureAgent 
from src .agents .gwas_agent import GwasAgent 
from src .agents .mirna_sequence_agent import MiRnaSequenceAgent 
from src .utils .rna_type_detector import RnaTypeDetector 
from src .utils .cross_omics_engine import CrossOmicsEngine 
from src .project_paths import ANNOTATION_DIR ,PROCESSED_DIR ,NCRNA_DIR 

class MiRnaTraitAgent :
    def __init__ (self ):
        self .console =Console ()

        #Phase 0: Ensure the Database exists (must be done before initializing NcRnaAgent)
        self ._ensure_clean_database ()

        #Initialize GeneAnnotationSystem
        self .gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")
        self .map_file =str (ANNOTATION_DIR /"RAP-MSU.txt")
        self .annotator =RiceGeneAnnotator (self .gff_file ,self .map_file )
        self .annotator .build_database ()

        #Initialize the four major engines
        self .ncrna_agent =NcRnaAgent ()
        self .lit_agent =LiteratureAgent ()
        self .gwas_agent =GwasAgent ()
        self .seq_agent =MiRnaSequenceAgent (self .annotator )
        self .detector =RnaTypeDetector ()
        self .cross_engine =CrossOmicsEngine ()

    @staticmethod 
    def _to_float (value ,default =999.0 ):
        try :
            if value in (None ,"","NA","nan"):
                return default 
            return float (value )
        except (TypeError ,ValueError ):
            return default 

    @classmethod 
    def _database_rank_score (cls ,target ):
        degradome_score =cls ._to_float (target .get ("Degradome_score"),default =999.0 )
        psr_score =cls ._to_float (target .get ("psRNAtarget_score"),default =999.0 )
        lib_count =cls ._to_float (target .get ("Lib_count"),default =0.0 )

        if degradome_score >=999 and psr_score >=999 :
            return 0.0 

        evidence =max (0.0 ,60.0 -degradome_score *10.0 )
        pairing =max (0.0 ,30.0 -psr_score *5.0 )
        support =min (lib_count ,10.0 )
        return round (evidence +pairing +support ,2 )

    @classmethod 
    def _format_target_evidence (cls ,target ):
        degradome =str (target .get ("Degradome","N/A")).strip ()
        psr =str (target .get ("psRNAtarget","N/A")).strip ()

        if degradome .endswith (".0"):
            degradome =degradome [:-2 ]
        if psr .endswith (".0"):
            psr =psr [:-2 ]

        rank_score =cls ._database_rank_score (target )
        return f"Score:{rank_score }, Deg:{degradome }, psR:{psr }"

    def _ensure_clean_database (self ):
        """Automatically detect and generate cleaned database"""
        if not (PROCESSED_DIR /"clean_mirna_targets.csv").exists ():
            self .console .print ("[yellow]⚠️ Database not initialized detected, currently auto-cleaning...[/yellow]")
            try :
                subprocess .run ([sys .executable ,str (PROJECT_ROOT /"scripts"/"build_mirna_target_db.py")],check =True )
                self .console .print ("[green]✅ Database initialization completed![/green]")
            except subprocess .CalledProcessError as e :
                self .console .print (f"[bold red]❌ DatabaseinitializeFailed: {e }[/bold red]")
                sys .exit (1 )

    def run_phasirna_mode (self ,phasirna_id :str ,top_n :int =5 ):
        """phasiRNA Analysis workflow (ID-driven, prioritizing negative correlation in expression, falling back to plant scoring complementarity)"""
        sequence =self .detector .get_sequence_from_id (phasirna_id )
        if not sequence :
            self .console .print (f"[red]Not Found {phasirna_id } 'sSequence，pleaseCheck ID format[/red]")
            return 

        self .console .print (Panel (
        f"[bold cyan]sRNA IntelligentAnalysisdevice[/bold cyan]\ntarget: {phasirna_id }\n"
        f"Sequence: {sequence } ({len (sequence )}nt)\npattern: phasiRNA — Expressionnegatively correlated + SequenceComplementary",
        border_style ="cyan"
        ))

        #Phase 1: Candidate target genes (prioritizing negative correlation)
        self .console .print ("\n[bold cyan]Part 1: Candidate Target Gene Identification[/bold cyan]")
        self .console .print ("[cyan]Prioritize using Expression matrix negative correlation Relationship to screen candidate targets...[/cyan]")

        partners_df =self .cross_engine .find_cross_omics_partners (phasirna_id ,top_n =top_n )
        targets =[]
        evidence_mode ="expression_anticorrelation"

        if partners_df is not None and not partners_df .empty :
            for _ ,row in partners_df .head (top_n ).iterrows ():
                gene_info =self .annotator .query (row ["Gene_ID"])
                if not gene_info :
                    continue 
                symbol =gene_info .get ("Primary_Symbol")or (
                gene_info .get ("Synonyms",[""])[0 ]if gene_info .get ("Synonyms")else "No"
                )
                targets .append ({
                "RAP_ID":gene_info ["RAP_ID"],
                "Primary_Symbol":symbol ,
                "PMIDs":gene_info .get ("PMIDs",[]),
                "Correlation":float (row ["Correlation"]),
                "Evidence_Type":"expression_anticorrelation",
                })

        if not targets :
            self .console .print ("[yellow]⚠ Expression matrix No stable candidate, downgraded to plant standard sequence complementary scan.[/yellow]")
            evidence_mode ="plant_complementarity"
            pred_df =self .seq_agent .predict_targets (sequence ,top_n =top_n ,rna_label ="phasiRNA")
            for _ ,row in pred_df .iterrows ():
                targets .append ({
                "RAP_ID":row ["Gene_ID"],
                "Primary_Symbol":row ["Symbol"],
                "PMIDs":row ["PMIDs"],
                "Penalty":row ["Penalty"],
                "Evidence_Type":"plant_complementarity",
                })

        if not targets :
            self .console .print ("[red]Not Found candidate target Gene[/red]")
            return 

        target_ids =[]
        target_rows =[]
        for i ,t in enumerate (targets ,1 ):
            symbol =t .get ("Primary_Symbol","N/A")
            has_paper ="📚"if t .get ("PMIDs")else "🌌"
            if evidence_mode =="expression_anticorrelation":
                corr =t .get ("Correlation",0.0 )
                self .console .print (f"  {i }. {t ['RAP_ID']} ({symbol }) | AntiCorr={corr :.3f} {has_paper }")
                target_rows .append ({"rap_id":t ["RAP_ID"],"symbol":symbol ,
                "anti_corr":round (corr ,3 ),"has_pmid":bool (t .get ("PMIDs"))})
            else :
                pen =t .get ("Penalty",0 )
                self .console .print (f"  {i }. {t ['RAP_ID']} ({symbol }) | Penalty={pen :.1f} {has_paper }")
                target_rows .append ({"rap_id":t ["RAP_ID"],"symbol":symbol ,
                "penalty":pen ,"has_pmid":bool (t .get ("PMIDs"))})
            target_ids .append (t ["RAP_ID"])

        _emit_step ("phasirna_targets","Candidate target gene identification","Target Gene Identification",{
        "phasirna_id":phasirna_id ,"sequence":sequence ,
        "evidence_mode":evidence_mode ,"target_count":len (targets ),"targets":target_rows ,
        })

        # Phase 2: Literature
        self .console .print ("\n[bold cyan]Part 2: Target Gene Literature Fetch[/bold cyan]")
        abstracts =self .lit_agent .fetch_abstracts_route_b (targets )

        #Phase 3: Reverse Inference
        self .console .print ("\n[bold cyan]Part Three: Reverse Inference of phasiRNA Function[/bold cyan]")
        if abstracts and abstracts .strip ():
            report =self .lit_agent .summarize_ncrna (phasirna_id ,targets ,abstracts ,rna_type ="phasiRNA")
            self .console .print (Markdown (report ))
        else :
            report =""
            self .console .print ("[yellow]Check NoLiterature, skip FunctionInference[/yellow]")
        _emit_step ("phasirna_function","Function reverse inference","Functional Reverse Inference",{"report":report })

        # Phase 4: GWAS
        self .console .print ("\n[bold cyan]Part Four: GWAS Breeding Value Evaluation[/bold cyan]")
        gwas_df =self .gwas_agent .scan_network_gwas (target_ids ,self .annotator )
        if gwas_df .empty :
            self .console .print ("[yellow]No significant GWAS signal detected[/yellow]")
            gwas_report ,gwas_records ="",[]
        else :
            gwas_report =self .gwas_agent .summarize_ncrna (phasirna_id ,target_ids ,gwas_df )
            self .console .print (Markdown (gwas_report ))
            gwas_records =gwas_df .to_dict ("records")
        _emit_step ("phasirna_gwas","GWAS breeding value evaluation","GWAS Breeding Assessment",
        {"report":gwas_report ,"gwas_signals":gwas_records })

        self .console .print (Panel (f"[bold green]{phasirna_id } AnalysisCompleted[/bold green]",border_style ="green"))

    def run_database_mode (self ,mirna_id ,top_n =3 ):
        """Database target gene prediction mode"""
        self .console .print (Panel (
        f"[bold cyan]miRNA IntelligentAnalysisdevice[/bold cyan]\ntarget: {mirna_id }\npattern: Databasequery",
        border_style ="cyan"
        ))

        #Phase 1: Target Gene Query
        self .console .print (f"\n[bold cyan]first part：targetGenequery[/bold cyan]")
        targets =self .ncrna_agent .find_top_targets (mirna_id ,self .annotator ,top_n =top_n )

        if not targets :
            self .console .print (f"[bold red]Not Found {mirna_id } target ofGene[/bold red]")
            sys .exit (1 )

        annotated =[t for t in targets if t .get ('Primary_Symbol')and t ['Primary_Symbol']!=t ['RAP_ID']]
        self .console .print (f"found {len (targets )} target siteGene（{len (annotated )} eachHasAnnotation）")
        target_ids =[]
        target_rows =[]
        for t in targets :
            symbol =t .get ('Primary_Symbol')or (t .get ('Synonyms',[''])[0 ]if t .get ('Synonyms')else 'No')
            has_paper ="📚"if t .get ('PMIDs')else "🌌"
            evidence =self ._format_target_evidence (t )
            self .console .print (f"  • {t ['RAP_ID']} ({symbol }) [{evidence }] {has_paper }")
            target_ids .append (t ['RAP_ID'])
            target_rows .append ({
            "rap_id":t ['RAP_ID'],
            "symbol":symbol ,
            "evidence":evidence ,
            "has_pmid":bool (t .get ('PMIDs')),
            })

        self .console .print ("ℹ️ Score is a comprehensive ranking score (Degradome + psRNATarget + library support), with lower Deg indicating stronger cleavage evidence.")
        _emit_step ("mirna_targets","Target Gene Query","Target Gene Query",{
        "mirna_id":mirna_id ,
        "mode":"database",
        "target_count":len (targets ),
        "annotated_count":len (annotated ),
        "targets":target_rows ,
        })

        # Phase 2: LiteratureFetch
        self .console .print (f"\n[bold cyan]part two：targetGeneLiteratureSearch[/bold cyan]")
        abstracts =self .lit_agent .fetch_abstracts_route_b (targets )

        #Phase 3: Reverse Inference
        self .console .print (f"\n[bold cyan]third part：miRNA FunctionreverseInference[/bold cyan]")
        if abstracts .strip ():
            mechanism_report =self .lit_agent .summarize_ncrna (mirna_id ,targets ,abstracts ,rna_type ="miRNA")
            self .console .print (Markdown (mechanism_report ))
        else :
            mechanism_report =""
            self .console .print ("[yellow]Check NoLiterature, skip FunctionInference[/yellow]")
        _emit_step ("mirna_function","Function reverse inference","Functional Reverse Inference",{
        "report":mechanism_report ,
        })

        # Phase 4: GWAS Evaluation
        self .console .print (f"\n[bold cyan]fourth part：GWAS breeding effectEvaluation[/bold cyan]")
        self .console .print ("Currently scanning GWAS signals for the target gene...")
        gwas_df =self .gwas_agent .scan_network_gwas (target_ids ,self .annotator )

        if not gwas_df .empty :
            self .console .print (f"capture [bold green]{len (gwas_df )}[/bold green] a significant GWAS Signal，CurrentlyEvaluation miRNA passedExpressionbreeding effect of...")
            yield_report =self .gwas_agent .summarize_ncrna (mirna_id ,target_ids ,gwas_df )
            self .console .print (Markdown (yield_report ))
            gwas_records =gwas_df .to_dict ("records")
        else :
            yield_report =""
            gwas_records =[]
            self .console .print ("[yellow]No significant GWAS signal found[/yellow]")
        _emit_step ("mirna_gwas","GWAS breeding effect evaluation","GWAS Breeding Assessment",{
        "report":yield_report ,
        "gwas_signals":gwas_records ,
        })

        self .console .print (Panel (f"[bold green]{mirna_id } AnalysisCompleted[/bold green]",border_style ="green"))

    def run_sequence_mode (self ,mirna_seq ,top_n =10 ,known_id =None ):
        """Sequence prediction mode (plant standard scoring)"""
        label =f"Hasknown {known_id }"if known_id else "New Sequence"
        self .console .print (Panel (
        f"[bold cyan]miRNA SequencepredictSystem[/bold cyan]\nSequence: {mirna_seq }\nlength: {len (mirna_seq )}nt\npattern: Sequencescoring prediction（{label }）",
        border_style ="cyan"
        ))

        #Phase 1: Sequence prediction
        self .console .print (f"\n[bold cyan]first part：targetGenepredict[/bold cyan]")
        try :
            pred_df =self .seq_agent .predict_targets (mirna_seq ,top_n =top_n )
        except ValueError as e :
            self .console .print (f"[bold red]SequenceValidationFailed: {e }[/bold red]")
            sys .exit (1 )

        if pred_df .empty :
            self .console .print ("[bold red]Not Found candidate target Gene[/bold red]")
            sys .exit (1 )

        self .console .print (f"predicted to {len (pred_df )} candidate targetsGene")
        targets =[]
        target_ids =[]
        target_rows =[]
        for _ ,row in pred_df .iterrows ():
            symbol =row ['Symbol']
            penalty =row ['Penalty']
            has_paper ="📚"if row ['PMIDs']else "🌌"
            self .console .print (f"  • {row ['Gene_ID']} ({symbol }) [Penalty:{penalty :.1f}] {has_paper }")

            targets .append ({
            'RAP_ID':row ['Gene_ID'],
            'Primary_Symbol':symbol ,
            'PMIDs':row ['PMIDs']
            })
            target_ids .append (row ['Gene_ID'])
            target_rows .append ({
            "rap_id":row ['Gene_ID'],
            "symbol":symbol ,
            "penalty":penalty ,
            "has_pmid":bool (row ['PMIDs']),
            })

        self .console .print ("ℹ️ Penalty is the plant standard penalty (≤4.0 retained, lower indicates a more perfect pairing, 0 = completely complementary).")
        _emit_step ("mirna_targets","Target Gene Prediction","Target Gene Prediction",{
        "mirna_id":known_id or "novel-sequence",
        "mode":"sequence",
        "sequence":mirna_seq ,
        "target_count":len (targets ),
        "annotated_count":len (targets ),
        "targets":target_rows ,
        })

        # Phase 2: LiteratureFetch
        self .console .print (f"\n[bold cyan]part two：targetGeneLiteratureSearch[/bold cyan]")
        abstracts =self .lit_agent .fetch_abstracts_route_b (targets )

        #Phase 3: Reverse Inference
        self .console .print (f"\n[bold cyan]third part：miRNA FunctionreverseInference[/bold cyan]")
        rna_label =known_id or "novel-sequence"
        if abstracts .strip ():
            mechanism_report =self .lit_agent .summarize_ncrna (rna_label ,targets ,abstracts )
            self .console .print (Markdown (mechanism_report ))
        else :
            mechanism_report =""
            self .console .print ("[yellow]Check NoLiterature, skip FunctionInference[/yellow]")
        _emit_step ("mirna_function","Function reverse inference","Functional Reverse Inference",{
        "report":mechanism_report ,
        })

        # Phase 4: GWAS Evaluation
        self .console .print (f"\n[bold cyan]fourth part：GWAS breeding effectEvaluation[/bold cyan]")
        self .console .print ("Currently scanning GWAS signals for the target gene...")
        gwas_df =self .gwas_agent .scan_network_gwas (target_ids ,self .annotator )

        if not gwas_df .empty :
            self .console .print (f"capture [bold green]{len (gwas_df )}[/bold green] a significant GWAS Signal，CurrentlyEvaluation miRNA passedExpressionbreeding effect of...")
            yield_report =self .gwas_agent .summarize_ncrna (rna_label ,target_ids ,gwas_df )
            self .console .print (Markdown (yield_report ))
            gwas_records =gwas_df .to_dict ("records")
        else :
            yield_report =""
            gwas_records =[]
            self .console .print ("[yellow]No significant GWAS signal found[/yellow]")
        _emit_step ("mirna_gwas","GWAS breeding effect evaluation","GWAS Breeding Assessment",{
        "report":yield_report ,
        "gwas_signals":gwas_records ,
        })

        self .console .print (Panel ("[bold green]miRNA Sequence Prediction Completed[/bold green]",border_style ="green"))

def _is_sequence (s ):
    """Determine whether Input is a Sequence or an ID"""
    return all (c in 'AUCGTaucgt'for c in s .strip ())


if __name__ =="__main__":
    parser =argparse .ArgumentParser (description ="Rice miRNA Trait Gene Analysis System")
    parser .add_argument ("--gene",required =True ,
    help ="miRNA ID (e.g., osa-mir156a) or mature sequence (e.g., UGACAGAAGAGAGUGAGCAC)")
    parser .add_argument ("--top-n",type =int ,default =3 ,help ="Return target gene number (default 3)")
    args =parser .parse_args ()

    agent =MiRnaTraitAgent ()
    inp =args .gene .strip ()

    #── Branch 1: Input is a Sequence ────────────────────────────────────────────────
    if _is_sequence (inp ):
        seq =inp .upper ().replace ('T','U')
        seq_len =len (seq )

        #Length ≤ 25 nt: Mature form, proceed directly to prediction pipeline
        if seq_len <=25 :
        #Search library with annotation (exact match)
            db_type =agent .detector .detect_by_database (seq )
            known_id =None 
            if db_type =="miRNA":
            #Reverse lookup ID (search FASTA for corresponding ID)
                for fasta in ('mirbase_osa_mature.fa','novel_mature_miRNA.fa'):
                    sd =agent .detector ._load_fasta (fasta )
                    for sid ,sseq in sd .items ():
                        if sseq ==seq :
                            known_id =sid 
                            break 
                    if known_id :
                        break 
                agent .console .print (
                f"[green]✓ HitDatabaseHasknown miRNA：{known_id }（{seq_len }nt），additionalHasknownAnnotationcontinue prediction after[/green]"
                )
            elif db_type =="phasiRNA":
                agent .console .print (
                f"[green]✓ HitDatabaseHasknown phasiRNA（{seq_len }nt），additionalHasknownAnnotationcontinue prediction after[/green]"
                )
            else :
                agent .console .print (
                f"[cyan]ℹ SequencemissedDatabase，as newSequencemake predictions（{seq_len }nt）[/cyan]"
                )
            agent .run_sequence_mode (seq ,args .top_n ,known_id =known_id )

            #Length > 25 nt: may be a precursor, check the corresponding mature form in the database
        else :
            agent .console .print (
            f"[yellow]⚠ InputSequencelength of {seq_len }nt，over mature miRNA range（≤25nt），may be precursorSequence。[/yellow]"
            )
            agent .console .print ("[cyan]The corresponding mature miRNA in CurrentlySearchDatabase...[/cyan]")

            #Precursor exact match search in database to find corresponding mature form
            pre_dict =agent .detector ._load_fasta ('novel_precursor.fa')
            mature_id =None 
            for pid ,pseq in pre_dict .items ():
                if pseq ==seq :
                #Look for a mature counterpart with the same name (try matching by removing suffixes from the ID).
                    base_id =pid .rsplit ('_',1 )[0 ]if '_'in pid else pid 
                    mature_seq =agent .detector .get_sequence_from_id (base_id )
                    if mature_seq and len (mature_seq )<=25 :
                        mature_id =base_id 
                        seq =mature_seq .upper ().replace ('T','U')
                        break 

            if mature_id :
                agent .console .print (
                f"[green]✓ find corresponding mature miRNA：{mature_id }（{len (seq )}nt），will use mature formsSequenceConductAnalysis。[/green]"
                )
                agent .run_sequence_mode (seq ,args .top_n ,known_id =mature_id )
            else :
                agent .console .print (
                "[red]✗ No mature miRNA record found in the database corresponding to Not Found. Please provide the mature sequence (≤25nt) and try again.[/red]"
                )
                sys .exit (0 )

                #── Branch 2: Input is ID ─────────────────────────────────────────────────
    else :
        rna_type =agent .detector .detect_by_id (inp )
        if rna_type =="phasiRNA":
            agent .run_phasirna_mode (inp ,args .top_n )
        else :
            agent .run_database_mode (inp .lower (),args .top_n )
