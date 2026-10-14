#!/usr/bin/env python3
"""
Unknown Sequence Identification Tool

Helps users identify the type of unknown RNA sequence and recommends appropriate analysis scripts.
"""

import argparse 
import sys 
from pathlib import Path 
from rich .console import Console 
from rich .panel import Panel 
from rich .table import Table 

PROJECT_ROOT =Path (__file__ ).resolve ().parents [1 ]
if str (PROJECT_ROOT )not in sys .path :
    sys .path .insert (0 ,str (PROJECT_ROOT ))

from src .utils .rna_type_detector import RnaTypeDetector 


def main ():
    parser =argparse .ArgumentParser (
    description ="Identify unknown RNA sequence type",
    formatter_class =argparse .RawDescriptionHelpFormatter ,
    epilog ="""
Example usage:
  python apps/sequence_identifier.py --sequence ACCAGTTGACGTGGCGCTGA
  python apps/sequence_identifier.py --id Osa-miRN01_RiceNCexp
  python apps/sequence_identifier.py --id Osa-miRN01 --sequence ACCAGTTGACGTGGCGCTGA
        """
    )

    parser .add_argument ("--sequence",help ="RNA Sequence")
    parser .add_argument ("--id",help ="RNA ID (optional)")

    args =parser .parse_args ()

    if not args .sequence and not args .id :
        parser .error ("Either --sequence or --id parameter must be provided")

    console =Console ()
    detector =RnaTypeDetector ()

    #If only ID is provided, attempt to load sequence from database
    sequence =args .sequence 
    if not sequence and args .id :
        sequence =detector .get_sequence_from_id (args .id )
        if sequence :
            console .print (f"[green]✓ fromDatabaseloadSequence: {sequence }[/green]\n")

            #Identify type
    rna_type ,method =detector .detect (rna_id =args .id ,sequence =sequence )

    #Display result
    console .print (Panel .fit (
    "[bold cyan]RNA Sequence Identification Result[/bold cyan]",
    border_style ="cyan"
    ))

    table =Table (show_header =True ,header_style ="bold magenta")
    table .add_column ("Project",style ="cyan")
    table .add_column ("Information",style ="green")

    if args .id :
        table .add_row ("Input ID",args .id )
    if sequence :
        table .add_row ("Sequence",sequence )
        table .add_row ("Length",f"{len (sequence )} nt")

    if rna_type :
        table .add_row ("Identified type",rna_type )
        table .add_row ("Identification method",method )

        #Recommended script
        script_map ={
        'miRNA':'apps/mirna_pipeline.py',
        'miRNA_precursor':'apps/mirna_pipeline.py',
        'phasiRNA':'apps/mirna_pipeline.py',
        'lncRNA':'apps/lncrna_pipeline.py'
        }

        if rna_type in script_map :
            script =script_map [rna_type ]
            table .add_row ("Recommended script",f"[bold yellow]{script }[/bold yellow]")

            #Generate example command
            if args .id :
                example_cmd =f"python {script } --gene {args .id }"
            elif sequence and rna_type in ['miRNA','phasiRNA']:
                example_cmd =f"python {script } --mode sequence --sequence {sequence }"
            else :
                example_cmd =f"python {script } --gene <ID>"

            table .add_row ("Example command",example_cmd )
    else :
        table .add_row ("Identified type","[red]Unable to identify[/red]")

    console .print (table )
    console .print ()


if __name__ =="__main__":
    main ()
