import os 
import urllib .parse 
from collections import Counter 
from rich .console import Console 
from rich .table import Table 
from src .project_paths import ANNOTATION_DIR 

console =Console ()

def explore_gff_attributes (gff_path ):
    console .print (f"🔍 CurrentlyScantreasure trove: {gff_path } ...\n")

    known_genes =[]
    unknown_genes =[]

    known_keys_counter =Counter ()
    unknown_keys_counter =Counter ()

    with open (gff_path ,'r')as f :
        for line in f :
            if line .startswith ('#'):continue 
            parts =line .strip ().split ('\t')

            #We only focus on mRNA level because rich information is concentrated here
            if len (parts )<9 or parts [2 ]!='mRNA':continue 

            attr_string =parts [8 ]
            attrs ={}
            for item in attr_string .split (';'):
                if '='in item :
                    k ,v =item .split ('=',1 )
                    attrs [k ]=urllib .parse .unquote (v )#Decode web character %2C

                    #Core decision logic: Whether it contains PMID or DOI
            is_known ='Literature_PMID'in attrs 

            if is_known :
                known_genes .append (attrs )
                for k in attrs .keys ():
                    known_keys_counter [k ]+=1 
            else :
                unknown_genes .append (attrs )
                for k in attrs .keys ():
                    unknown_keys_counter [k ]+=1 

                    #================= Print Statistics Report =================
    console .print (f"📊 [bold green]ScanResultoverview[/bold green]")
    console .print (f"✅ discover 'HasknownGene' (withLiterature PMID): {len (known_genes )} each")
    console .print (f"❓ discover 'unknown/uncharacterizedGene' (NoLiterature): {len (unknown_genes )} each\n")

    #Print the most frequently carried treasure fields of Known Genes
    table =Table (title ="💎 Known Genes (with PMIDs) field frequency report",show_header =True ,header_style ="bold magenta")
    table .add_column ("Attribute Field",style ="cyan")
    table .add_column ("Occurrence Count",justify ="right")
    table .add_column ("Coverage Rate",justify ="right")

    for k ,v in known_keys_counter .most_common (15 ):
        coverage =(v /len (known_genes ))*100 
        table .add_row (k ,str (v ),f"{coverage :.1f}%")
    console .print (table )

    #Display one exceptionally rich Known Gene at random
    if known_genes :
        showcase =known_genes [0 ]#Select the first one to display
        console .print (f"\n🌟 [bold yellow]luxury archive display (Locus: {showcase .get ('Locus_id')})[/bold yellow]")
        for k ,v in showcase .items ():
            if k not in ['ID','Name']:#Omit some No. chat IDs
                console .print (f"  [bold cyan]{k }:[/bold cyan] {v }")

if __name__ =="__main__":
#Use the transcripts.gff (or transcripts_exon.gff) you discovered
    gff_file =str (ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts.gff")

    if not os .path .exists (gff_file ):
        console .print (f"[red]❌ cannot findFile: {gff_file }[/red]")
    else :
        explore_gff_attributes (gff_file )
