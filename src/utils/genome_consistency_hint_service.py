from __future__ import annotations 

import gzip 
from functools import lru_cache 
from pathlib import Path 
from typing import Any 

from src .project_paths import ANNOTATION_DIR 


class GenomeConsistencyHintService :
    def __init__ (
    self ,
    transcripts_exon_gff :str |None =None ,
    external_rna_fasta :str ="",
    assembly_fasta :str ="",
    ):
        self .transcripts_exon_gff =transcripts_exon_gff or str (
        ANNOTATION_DIR /"IRGSP-1.0_representative"/"transcripts_exon.gff"
        )
        self .external_rna_fasta =external_rna_fasta 
        self .assembly_fasta =assembly_fasta 

    def collect (self ,gene_info :dict [str ,Any ])->dict [str ,Any ]:
        transcript_ids =[item for item in gene_info .get ("Transcript_IDs",[])if item ]
        msu_ids =[item for item in gene_info .get ("MSU_IDs",[])if item ]
        symbol =(gene_info .get ("Primary_Symbol")or "").strip ()
        rap_id =gene_info .get ("RAP_ID","")

        result ={
        "gene_id":rap_id ,
        "transcript_ids":transcript_ids ,
        "msu_ids":msu_ids ,
        "assembly_available":False ,
        "external_rna_match_count":0 ,
        "external_rna_matches":[],
        "warnings":[],
        "summary_text":"",
        }

        assembly_path =Path (self .assembly_fasta )
        if assembly_path .exists ():
            try :
                transcript_sequences =self ._extract_transcript_sequences (
                transcript_ids ,
                gene_info ["Chromosome"],
                gene_info ["Strand"],
                )
                if transcript_sequences :
                    result ["assembly_available"]=True 
                else :
                    result ["warnings"].append ("Failed to extract the Transcript sequence for this gene from OSA1 assembly.")
            except Exception as exc :
                result ["warnings"].append (f"read OSA1 assembleFailed: {exc }")
        else :
            result ["warnings"].append (f"Not Found OSA1 assembleFile: {assembly_path }")

        external_path =Path (self .external_rna_fasta )
        if external_path .exists ():
            try :
                matches =self ._search_external_rna_headers (rap_id ,msu_ids ,transcript_ids ,symbol )
                result ["external_rna_matches"]=matches [:10 ]
                result ["external_rna_match_count"]=len (matches )
            except Exception as exc :
                result ["warnings"].append (f"ScanExternal RNA cacheFailed: {exc }")
        else :
            result ["warnings"].append (f"Not FoundExternal RNA File: {external_path }")

        lines =["[Cross-Sequence Resource Consistency Notice]"]
        lines .append (f"- GenegroupAnnotationTranscriptnumber: {len (transcript_ids )}")
        lines .append (f"- Expressionmatrix model number: {len (msu_ids )}")
        if result ["assembly_available"]:
            lines .append ("- OSA1 assembly: Successfully extracted the corresponding Transcript/Gene region sequence.")
        else :
            lines .append ("- OSA1 assembly: Sequence extraction not yet successful.")

        if result ["external_rna_match_count"]>0 :
            lines .append (f"- GCF RNA cache: automatically matched to {result ['external_rna_match_count']} candidate entriesTranscriptheadInformation。")
            for header in result ["external_rna_matches"][:3 ]:
                lines .append (f"  * {header }")
            lines .append ("- Auto-notice: External RNA candidate detected; please manually verify if it belongs to the same gene.")
        else :
            lines .append ("- GCF RNA library: Items not found can be directly used via RAP/MSU/symbol auto-mapping.")
            lines .append ("- Auto-notice: Currently unable to automatically determine full consistency across the three sequence resources; for rigorous conclusions, manual comparison is required.")

        if result ["warnings"]:
            for warning in result ["warnings"]:
                lines .append (f"- warning: {warning }")

        result ["summary_text"]="\n".join (lines )
        return result 

    def _extract_transcript_sequences (
    self ,
    transcript_ids :list [str ],
    chromosome :str ,
    strand :str ,
    )->dict [str ,str ]:
        exons_by_transcript =self ._load_exons (transcript_ids )
        if not exons_by_transcript :
            return {}
        fasta =self ._load_assembly_fasta (self .assembly_fasta )
        chrom_key =self ._normalize_chromosome (chromosome )
        chrom_seq =fasta .get (chrom_key )
        if not chrom_seq :
            return {}

        sequences ={}
        for transcript_id in transcript_ids :
            exons =exons_by_transcript .get (transcript_id ,[])
            if not exons :
                continue 
            ordered =sorted (exons ,key =lambda item :item [0 ],reverse =(strand =="-"))
            seq ="".join (chrom_seq [start -1 :end ]for start ,end in ordered )
            if strand =="-":
                seq =self ._reverse_complement (seq )
            sequences [transcript_id ]=seq 
        return sequences 

    def _search_external_rna_headers (
    self ,
    rap_id :str ,
    msu_ids :list [str ],
    transcript_ids :list [str ],
    symbol :str ,
    )->list [str ]:
        patterns =[rap_id ,*msu_ids ,*transcript_ids ]
        if symbol :
            patterns .append (symbol )
        patterns =[item .lower ()for item in patterns if item ]
        matches =[]
        with gzip .open (self .external_rna_fasta ,"rt")as handle :
            for line in handle :
                if not line .startswith (">"):
                    continue 
                lowered =line .lower ()
                if any (pattern in lowered for pattern in patterns ):
                    matches .append (line .strip ())
        return matches 

    @lru_cache (maxsize =1 )
    def _load_assembly_fasta (self ,path :str )->dict [str ,str ]:
        sequences :dict [str ,list [str ]]={}
        current =None 
        with gzip .open (path ,"rt")as handle :
            for line in handle :
                if line .startswith (">"):
                    current =line [1 :].strip ().split ()[0 ]
                    sequences [current ]=[]
                elif current is not None :
                    sequences [current ].append (line .strip ().upper ())
        return {key :"".join (value )for key ,value in sequences .items ()}

    def _load_exons (self ,transcript_ids :list [str ])->dict [str ,list [tuple [int ,int ]]]:
        exon_map ={tid :[]for tid in transcript_ids }
        with open (self .transcripts_exon_gff ,"r")as handle :
            for line in handle :
                if line .startswith ("#")or "\texon\t"not in line :
                    continue 
                parts =line .rstrip ("\n").split ("\t")
                attrs =parts [8 ]
                for transcript_id in transcript_ids :
                    if f"Parent={transcript_id }"in attrs :
                        exon_map [transcript_id ].append ((int (parts [3 ]),int (parts [4 ])))
        return exon_map 

    @staticmethod 
    def _normalize_chromosome (chromosome :str )->str :
        chrom =chromosome .strip ()
        if chrom .lower ().startswith ("chr")and len (chrom )>=5 :
            suffix =chrom [3 :].lstrip ("0")
            return f"Chr{suffix or '0'}"
        return chrom 

    @staticmethod 
    def _reverse_complement (seq :str )->str :
        table =str .maketrans ("ACGTNacgtn","TGCANtgcan")
        return seq .translate (table )[::-1 ]
