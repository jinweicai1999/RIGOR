import argparse 
import sys 
from pathlib import Path 
from rich .console import Console 
from rich .panel import Panel 
from rich .markdown import Markdown 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

    #Import the exclusive Four Heavenly Kings engines from Protein GeneAnalysis
from src .utils .gene_annotator import RiceGeneAnnotator 
from src .agents .macro_network_agent import MacroNetworkAgent 
from src .agents .literature_agent import LiteratureAgent 
from src .agents .gwas_agent import GwasAgent 
from src .project_paths import ANNOTATION_DIR 

class ProteinTraitAgent :
    def __init__ (self ):
        self .console =Console ()

        #--- Initialization Phase 1: Household Registration System ---
        self .gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")
        self .map_file =str (ANNOTATION_DIR /"RAP-MSU.txt")
        self .annotator =RiceGeneAnnotator (self .gff_file ,self .map_file )
        self .annotator .build_database ()

        #--- Initialization Phase 2: Macro Network God's Eye View ---
        self .macro_agent =MacroNetworkAgent ()

        #--- Initialization Phase 3 & 4: Literature Hub and GWAS Yield Evaluation ---
        self .lit_agent =LiteratureAgent ()
        self .gwas_agent =GwasAgent ()

    def is_star_gene (self ,gene_info ):
        """Core decision logic: Whether the gene has officially certified PMID"""
        pmids =gene_info .get ('PMIDs',[])
        return len (pmids )>0 

    def run (self ,target_gene ):
        self .console .print (Panel (f"[bold cyan]🌾 proteinCodingGeneIntelligentAnalysisdevice (Protein-Portal) Start[/bold cyan]\ntargetGene: {target_gene }",border_style ="cyan"))

        # ==========================================
        #Phase 1: Household Registration and Routing
        # ==========================================
        self .console .print ("\n[bold green][Phase 1] Currently performing local gene identity verification...[/bold green]")
        gene_info =self .annotator .query (target_gene )

        if not gene_info :
            self .console .print (f"[bold red]❌ fatalError：Databasenot found in '{target_gene }'，pleaseConfirmInputisHaseffective rice proteinCodingGene。[/bold red]")
            sys .exit (1 )

        self .console .print (f"✅ identityVerifySuccessful: [bold]{gene_info ['RAP_ID']}[/bold] ({gene_info ['Location']})")

        self .console .print (f"  [cyan]• official brand (Primary Symbol):[/cyan] {gene_info .get ('Primary_Symbol','N/A')}")
        self .console .print (f"  [cyan]• Description (Note):[/cyan] {gene_info .get ('Description','N/A')}")
        syns =", ".join (gene_info .get ('Synonyms',[]))
        self .console .print (f"  [cyan]• alias (Synonyms):[/cyan] {syns if syns else 'No'}")
        self .console .print (f"  [cyan]• InterPro Domain:[/cyan] {gene_info .get ('InterPro','N/A')}")

        go_terms =gene_info .get ('GO',[])
        if go_terms :
            self .console .print ("  [cyan]• GO Annotation:[/cyan]")
            for go in go_terms :
                self .console .print (f"    - {go }")
        else :
            self .console .print ("  [cyan]• GO Annotation:[/cyan] No")

        pmids =gene_info .get ('PMIDs',[])
        self .console .print (f"  [cyan]• includedLiterature (PMID):[/cyan] {', '.join (pmids )if pmids else 'No'}")

        #Core routing and classification
        is_route_a =self .is_star_gene (gene_info )
        target_ids =[]
        abstracts =""

        if is_route_a :
            self .console .print (f"\n[bold magenta]🟢 routing determination：【pathway A - Hasknown starTarget】[/bold magenta]")
            self .console .print ("System determines that the gene has official literature (PMID), will directly extract core ontology literature.")
            target_ids =[gene_info ['RAP_ID']]
            abstracts =self .lit_agent .fetch_abstracts_route_a (pmids ,gene_info =gene_info )

        else :
            self .console .print (f"\n[bold yellow]🟡 routing determination：【pathway B - unknown dark matterTarget】[/bold yellow]")
            self .console .print ("System determines that the gene is a function unknown 'dark matter', currently starting Phase 2 macro network god's eye view expansion...")

            target_module ,total_same_mod ,same_mod_found ,network_genes =self .macro_agent .get_top_hub_partners (
            target_gene ,gene_info ,self .annotator ,top_n =5 
            )
            if not network_genes :
            #🌟 Delete the original sys.exit(1) and replace it with a gentle downgrade message!
                self .console .print ("[bold yellow]⚠️ Macro network does not include this gene (possibly due to extremely low expression or spatiotemporal specificity). Automatically downgraded to pure domain blind guess mode.[/bold yellow]")
                #Keep target_ids and abstracts empty, do not exit!
            else :
                friendly_names =[]
                #... The following traversal assembly logic remains unchanged ...

            friendly_names =[]
            target_ids =[]
            for d in network_genes :
                rap =d ['RAP_ID']
                syn_str =d .get ('Primary_Symbol')or (d .get ('Synonyms')[0 ]if d .get ('Synonyms')else "No")
                has_paper ="📚HasLiterature"if d .get ('PMIDs')else "🌌 Dark Matter"
                mod_tag ="🏠 Same School"if d .get ('Is_Same_Mod')else "🌍 External Ally"
                friendly_names .append (f"{rap } ({syn_str }) [{has_paper }|{mod_tag }]")
                target_ids .append (rap )

            if same_mod_found >=5 :
                self .console .print (f"🎯 totalHas [bold green]{total_same_mod }[/bold green] same school comrades，Successfulextract interaction weightsRankingprevious 5 and withHas Symbol core hub:\n   {friendly_names }")
            else :
                self .console .print (f"🎯 totalHas [bold green]{total_same_mod }[/bold green] same school comrades，found only in {same_mod_found } name with Symbol core。Hasautomatically borrow from outside，Successfulassemble {len (network_genes )} key node:\n   {friendly_names }")

            abstracts =self .lit_agent .fetch_abstracts_route_b (network_genes ,target_gene_info =gene_info )

            # ==========================================
            #Phase 3: PubMed Literature Intelligence and Mechanism Network Analysis
            # ==========================================
        self .console .print (f"\n[bold blue]🚀 Currentlywill target set {target_ids } transfer to large model center...[/bold blue]")
        self .console .print ("\n[bold cyan]--- [Phase 3] PubMed Literature Intelligence and Mechanism Network Analysis ---[/bold cyan]")

        if abstracts .strip ():
            self .console .print ("🧠 LiteratureFetch complete! Currently undergoing reading and deep inference by the large model hub...")
            mechanism_report =self .lit_agent .summarize (gene_info ,is_route_a ,abstracts )
            self .console .print ("\n"+"="*80 )
            self .console .print (Markdown (mechanism_report ))
            self .console .print ("="*80 +"\n")
        else :
            self .console .print ("[yellow]⚠️ No literature found, skipping mechanism inference step.[/yellow]")

            # ==========================================
            #Phase 4: Network Community Characteristics and GWAS Field Effect Evaluation
            # ==========================================
        self .console .print ("\n[bold cyan]--- [Phase 4] Network Community Characteristics and GWAS Field Effect Evaluation ---[/bold cyan]")

        # 1. Extract macro community characteristics
        macro_evidence =self .gwas_agent .get_macro_evidence (target_gene ,target_ids ,self .annotator )
        self .console .print ("[cyan]🔍 Network community phenotypic feature extraction completed.[/cyan]")
        self .console .print (macro_evidence )

        # 2. Start local radar scan for quantitative micro evidence
        self .console .print ("\n[cyan]📡 Currently starting local GWAS radar scan on target and its allies' quantitative effect values...[/cyan]")
        all_scan_targets =list (set ([target_gene ]+target_ids ))
        gwas_df =self .gwas_agent .scan_network_gwas (all_scan_targets ,self .annotator )

        if not gwas_df .empty :
            self .console .print (f"🎯 Successfulcapture [bold green]{len (gwas_df )}[/bold green] a significant GWAS quantitative associationSignal！")

            # 3. Transfer evidence to large model
            # 3. Transfer evidence to large model
        self .console .print ("\n🧠 Data collection complete! Currently evaluating yield cost and breeding application by large model...")
        yield_report =self .gwas_agent .summarize (target_gene ,target_ids ,gwas_df ,macro_evidence )

        # 🌟 The following four lines are the print statements I accidentally "lost"—let me fix them right away!
        self .console .print ("\n"+"="*80 )
        self .console .print (Markdown (yield_report ))
        self .console .print ("="*80 +"\n")
        self .console .print ("\n[bold green]✅ Protein-coding gene full process completed! (Protein-Portal Final Comprehensive Report Generated)[/bold green]")
if __name__ =="__main__":
    parser =argparse .ArgumentParser ()
    parser .add_argument ("--gene",type =str ,required =True ,help ="Input rice protein-coding gene for analysis (e.g., Os04g0372800)")
    args =parser .parse_args ()

    agent =ProteinTraitAgent ()
    agent .run (args .gene )
