from __future__ import annotations 

import os 
import re 
import math 
from typing import Dict ,List ,Optional 
from src .project_paths import CHROMATIN_DIR ,HISTONE_DIR 


class EpigenomicAnalyzer :
    def __init__ (
    self ,
    chromatin_dir :str =str (CHROMATIN_DIR ),
    histone_dir :str =str (HISTONE_DIR ),
    promoter_upstream :int =2000 ,
    promoter_downstream :int =500 ,
    ):
        self .chromatin_dir =chromatin_dir 
        self .histone_dir =histone_dir 
        self .promoter_upstream =promoter_upstream 
        self .promoter_downstream =promoter_downstream 

    def analyze (
    self ,
    gene_info :Optional [Dict ]=None ,
    coord :Optional [str ]=None ,
    strand :Optional [str ]=None ,
    )->Dict :
        promoter_region =self ._resolve_promoter_region (gene_info =gene_info ,coord =coord ,strand =strand )

        chromatin_hits =self ._scan_directory (self .chromatin_dir ,promoter_region ,"chromatin_accessibility")
        histone_hits =self ._scan_directory (self .histone_dir ,promoter_region ,"histone_modification")
        tissue_scores =self ._build_tissue_scores (chromatin_hits ,histone_hits )

        return {
        "promoter_region":promoter_region ,
        "chromatin_accessibility":chromatin_hits ,
        "histone_modification":histone_hits ,
        "tissue_scores":tissue_scores ,
        "summary_text":self ._build_summary_text (promoter_region ,chromatin_hits ,histone_hits ,tissue_scores ),
        }

    def _resolve_promoter_region (
    self ,
    gene_info :Optional [Dict ]=None ,
    coord :Optional [str ]=None ,
    strand :Optional [str ]=None ,
    )->Dict :
        if gene_info :
            chrom =self ._normalize_chromosome (gene_info ["Chromosome"])
            start =int (gene_info ["Start"])
            end =int (gene_info ["End"])
            gene_strand =gene_info .get ("Strand","+")
            gene_id =gene_info .get ("RAP_ID","unknown_gene")
        elif coord :
            chrom ,start ,end =self ._parse_coord (coord )
            gene_strand =strand or "+"
            gene_id =coord 
        else :
            raise ValueError ("Either gene_info or coord must be provided")

        if gene_strand not in {"+","-"}:
            gene_strand ="+"

        tss =start if gene_strand =="+"else end 
        if gene_strand =="+":
            promoter_start =max (1 ,tss -self .promoter_upstream )
            promoter_end =tss +self .promoter_downstream 
        else :
            promoter_start =max (1 ,tss -self .promoter_downstream )
            promoter_end =tss +self .promoter_upstream 

        return {
        "gene_id":gene_id ,
        "chromosome":chrom ,
        "gene_start":start ,
        "gene_end":end ,
        "strand":gene_strand ,
        "tss":tss ,
        "promoter_start":promoter_start ,
        "promoter_end":promoter_end ,
        "promoter_label":f"{chrom }:{promoter_start }-{promoter_end }",
        }

    def _parse_coord (self ,coord :str ):
        match =re .match (r"^(chr[0-9A-Za-z]+):(\d+)-(\d+)$",coord .strip ())
        if not match :
            raise ValueError (f"Nomethod to parse coordinate format: {coord }")
        chrom =self ._normalize_chromosome (match .group (1 ))
        start =int (match .group (2 ))
        end =int (match .group (3 ))
        if start >end :
            start ,end =end ,start 
        return chrom ,start ,end 

    def _normalize_chromosome (self ,chrom :str )->str :
        chrom =chrom .strip ()
        match =re .match (r"^chr0*([0-9]+)$",chrom ,flags =re .IGNORECASE )
        if match :
            return f"chr{int (match .group (1 ))}"
        return chrom 

    def _scan_directory (self ,directory :str ,promoter_region :Dict ,evidence_type :str )->List [Dict ]:
        if not os .path .isdir (directory ):
            return []

        hits =[]
        for file_name in sorted (os .listdir (directory )):
            if not file_name .endswith (".bed"):
                continue 
            file_path =os .path .join (directory ,file_name )
            overlap_info =self ._scan_bed_file (file_path ,promoter_region )
            if overlap_info :
                parsed_meta =self ._parse_bed_metadata (file_name ,overlap_info )
                parsed_meta ["evidence_type"]=evidence_type 
                hits .append (parsed_meta )
        return hits 

    def _scan_bed_file (self ,file_path :str ,promoter_region :Dict )->Optional [Dict ]:
        promoter_start0 =promoter_region ["promoter_start"]-1 
        promoter_end0 =promoter_region ["promoter_end"]
        promoter_chr =promoter_region ["chromosome"]
        overlaps =[]

        with open (file_path ,"r")as handle :
            for line in handle :
                if not line .strip ()or line .startswith ("#"):
                    continue 
                parts =line .rstrip ("\n").split ("\t")
                if len (parts )<3 :
                    continue 

                peak_chr =self ._normalize_chromosome (parts [0 ])
                if peak_chr !=promoter_chr :
                    continue 

                try :
                    peak_start =int (parts [1 ])
                    peak_end =int (parts [2 ])
                except ValueError :
                    continue 

                if peak_start >=promoter_end0 or peak_end <=promoter_start0 :
                    continue 

                overlap_bp =min (peak_end ,promoter_end0 )-max (peak_start ,promoter_start0 )
                overlaps .append (
                {
                "peak_chr":peak_chr ,
                "peak_start":peak_start ,
                "peak_end":peak_end ,
                "overlap_bp":overlap_bp ,
                "score":parts [4 ]if len (parts )>4 else "NA",
                "signal_value":parts [6 ]if len (parts )>6 else "NA",
                "peak_label":parts [8 ]if len (parts )>8 else os .path .basename (file_path ),
                "tissue_label":parts [9 ]if len (parts )>9 else "NA",
                }
                )

        if not overlaps :
            return None 

        overlaps .sort (key =lambda item :(item ["overlap_bp"],self ._safe_float (item ["signal_value"])),reverse =True )
        return {
        "file_path":file_path ,
        "hit_count":len (overlaps ),
        "best_peak":overlaps [0 ],
        "all_peaks":overlaps [:5 ],
        }

    def _parse_bed_metadata (self ,file_name :str ,overlap_info :Dict )->Dict :
        base_name =os .path .splitext (file_name )[0 ]
        tokens =base_name .split ("_")
        assay_or_mark =tokens [1 ]if len (tokens )>1 and tokens [0 ].startswith ("SRR")else tokens [0 ]
        tissue =tokens [-1 ].replace ("-"," ")if tokens else "unknown"
        mark_class =self ._classify_mark (assay_or_mark )

        return {
        "file_name":file_name ,
        "assay_or_mark":assay_or_mark ,
        "tissue":self ._normalize_tissue (tissue ),
        "mark_class":mark_class ,
        "hit_count":overlap_info ["hit_count"],
        "best_peak":overlap_info ["best_peak"],
        "all_peaks":overlap_info ["all_peaks"],
        }

    def _classify_mark (self ,assay_or_mark :str )->str :
        mark =assay_or_mark .lower ()
        if any (term in mark for term in ["atac","faire","mnase","h3k27ac","h3k4me3","h3k4me2","h3k9ac","h3k4ac","h4k12ac","h3k23ac","h4k16ac"]):
            return "active_or_open"
        if any (term in mark for term in ["h3k27me3","h3k9me2","h3k9me1"]):
            return "repressive"
        if "h3k36me3"in mark or "h3k4me1"in mark :
            return "gene_body_or_poised"
        return "context_dependent"

    def _build_tissue_scores (self ,chromatin_hits :List [Dict ],histone_hits :List [Dict ])->List [Dict ]:
        tissue_map :Dict [str ,Dict ]={}
        for hit in chromatin_hits +histone_hits :
            tissue =hit ["tissue"]
            tissue_entry =tissue_map .setdefault (
            tissue ,
            {
            "tissue":tissue ,
            "score":0.0 ,
            "open_evidence":[],
            "active_marks":[],
            "repressive_marks":[],
            "poised_marks":[],
            "context_marks":[],
            },
            )
            weight =self ._score_hit (hit )
            tissue_entry ["score"]+=weight 

            assay =hit ["assay_or_mark"]
            if assay .lower ()in {"atac-seq","faire","mnase-seq"}:
                tissue_entry ["open_evidence"].append (assay )
            elif hit ["mark_class"]=="active_or_open":
                tissue_entry ["active_marks"].append (assay )
            elif hit ["mark_class"]=="repressive":
                tissue_entry ["repressive_marks"].append (assay )
            elif hit ["mark_class"]=="gene_body_or_poised":
                tissue_entry ["poised_marks"].append (assay )
            else :
                tissue_entry ["context_marks"].append (assay )

        tissue_scores =[]
        for tissue ,tissue_entry in tissue_map .items ():
            tissue_entry ["state"]=self ._classify_tissue_state (tissue_entry )
            tissue_entry ["score"]=round (tissue_entry ["score"],2 )
            tissue_scores .append (tissue_entry )

        tissue_scores .sort (key =lambda item :item ["score"],reverse =True )
        return tissue_scores 

    def _score_hit (self ,hit :Dict )->float :
        best_peak =hit ["best_peak"]
        overlap =max (1 ,best_peak ["overlap_bp"])
        overlap_bonus =min (1.5 ,math .log10 (overlap +1 ))

        assay =hit ["assay_or_mark"].lower ()
        if assay in {"atac-seq","faire","mnase-seq"}:
            base =2.5 
        elif hit ["mark_class"]=="active_or_open":
            base =2.0 
        elif hit ["mark_class"]=="repressive":
            base =-2.0 
        elif hit ["mark_class"]=="gene_body_or_poised":
            base =0.75 
        else :
            base =0.25 

        return base +(overlap_bonus if base >=0 else -overlap_bonus *0.5 )

    def _classify_tissue_state (self ,tissue_entry :Dict )->str :
        score =tissue_entry ["score"]
        has_open =bool (tissue_entry ["open_evidence"])
        has_active =bool (tissue_entry ["active_marks"])
        has_repressive =bool (tissue_entry ["repressive_marks"])

        if (has_open or has_active )and has_repressive :
            return "dynamic_or_bivalent"
        if score >=5 :
            return "active"
        if score <=-2 :
            return "repressed"
        if has_open or has_active or tissue_entry ["poised_marks"]:
            return "poised_or_weakly_active"
        return "unclear"

    def _translate_mark (self ,mark :str )->str :
        """Translate technical terms into plain language names"""
        translations ={
        "H3K4me3":"Start site activation marker",
        "H3K27ac":"Active enhancer marker",
        "H3K9ac":"Transcription activation marker",
        "H3K4ac":"Transcription activation marker",
        "H3K23ac":"Transcription activation marker",
        "H4K16ac":"Transcription activation marker",
        "H4K12ac":"Transcription activation marker",
        "H3K27me3":"Gene silencing marker",
        "H3K9me2":"Heterochromatin inhibition marker",
        "H3K9me1":"Transcription inhibition marker",
        "H3K36me3":"Transcription elongation marker",
        "H3K4me1":"Enhancer marker",
        "ATAC-Seq":"Chromatin openness",
        "FAIRE":"Chromatin openness",
        "MNase-Seq":"Nucleosome positioning"
        }
        for key ,value in translations .items ():
            if key .lower ()in mark .lower ():
                return f"{mark }({value })"
        return mark 

    def _build_summary_text (self ,promoter_region :Dict ,chromatin_hits :List [Dict ],histone_hits :List [Dict ],tissue_scores :List [Dict ])->str :
        lines =[
        "【Start site epigenetic regulation evidence】",
        f"- Startsub-defense zone: {promoter_region ['promoter_label']} (TSS={promoter_region ['tss']}, strand={promoter_region ['strand']})",
        f"- chromatin openness hit count: {len (chromatin_hits )}",
        f"- histone modification hit count: {len (histone_hits )}",
        ]

        if chromatin_hits :
            lines .append ("- Openness evidence:")
            for hit in chromatin_hits [:8 ]:
                peak =hit ["best_peak"]
                readable_assay =self ._translate_mark (hit ['assay_or_mark'])
                lines .append (
                f"  * {hit ['file_name']} | detection type={readable_assay } | Tissue={hit ['tissue']} | "
                f"hit count={hit ['hit_count']} | best peak={peak ['peak_chr']}:{peak ['peak_start']}-{peak ['peak_end']} | "
                f"overlap={peak ['overlap_bp']} bp | Signalstrength={peak ['signal_value']}"
                )
        else :
            lines .append ("- Openness evidence: No start site peaks found in current chromatin openness data.")

        if histone_hits :
            lines .append ("- Histone modification evidence:")
            for hit in histone_hits [:12 ]:
                peak =hit ["best_peak"]
                readable_mark =self ._translate_mark (hit ['assay_or_mark'])
                class_translation ={
                "active_or_open":"Activation/Open",
                "repressive":"Inhibition",
                "gene_body_or_poised":"Gene body/Pending activation",
                "context_dependent":"Context-dependent"
                }
                readable_class =class_translation .get (hit ['mark_class'],hit ['mark_class'])
                lines .append (
                f"  * {hit ['file_name']} | modification type={readable_mark } | Function={readable_class } | Tissue={hit ['tissue']} | "
                f"hit count={hit ['hit_count']} | best peak={peak ['peak_chr']}:{peak ['peak_start']}-{peak ['peak_end']} | "
                f"overlap={peak ['overlap_bp']} bp | Signalstrength={peak ['signal_value']}"
                )
        else :
            lines .append ("- Histone modification evidence: No start site peaks found in current histone modification data.")

        if tissue_scores :
            lines .append ("- Start site activity scores aggregated by tissue:")
            for tissue_info in tissue_scores [:8 ]:
                state_translation ={
                "active":"Active state",
                "repressed":"Inhibited state",
                "poised_or_weakly_active":"Pending activation/Weakly active",
                "dynamic_or_bivalent":"Dynamic Regulation/Dual State",
                "unclear":"Unclear"
                }
                readable_state =state_translation .get (tissue_info ['state'],tissue_info ['state'])

                #Translation marker list
                open_marks =[self ._translate_mark (m )for m in sorted (set (tissue_info ['open_evidence']))]if tissue_info ['open_evidence']else []
                active_marks =[self ._translate_mark (m )for m in sorted (set (tissue_info ['active_marks']))]if tissue_info ['active_marks']else []
                repressive_marks =[self ._translate_mark (m )for m in sorted (set (tissue_info ['repressive_marks']))]if tissue_info ['repressive_marks']else []

                lines .append (
                f"  * Tissue={tissue_info ['tissue']} | score={tissue_info ['score']} | state={readable_state } | "
                f"openness={','.join (open_marks )or 'No'} | "
                f"activation marker={','.join (active_marks )or 'No'} | "
                f"Inhibitionmarker={','.join (repressive_marks )or 'No'}"
                )
        else :
            lines .append ("- Start sub-activity scores aggregated by Tissue: No hits, no scores available.")

        interpretation =self ._build_interpretation (chromatin_hits ,histone_hits )
        lines .append ("- Preliminary explanation:")
        lines .append (f"  * {interpretation }")
        return "\n".join (lines )

    def _build_interpretation (self ,chromatin_hits :List [Dict ],histone_hits :List [Dict ])->str :
        active_marks =[hit for hit in histone_hits if hit ["mark_class"]=="active_or_open"]
        repressive_marks =[hit for hit in histone_hits if hit ["mark_class"]=="repressive"]
        poised_marks =[hit for hit in histone_hits if hit ["mark_class"]=="gene_body_or_poised"]

        if chromatin_hits and active_marks and not repressive_marks :
            return "The promoter region shows both open chromatin and active histone modifications, suggesting this site has transcriptional activation potential in the corresponding tissue."
        if repressive_marks and not chromatin_hits and not active_marks :
            return "The promoter region primarily exhibits inhibitory histone modifications, indicating this site may be repressed or silenced in the corresponding tissue."
        if chromatin_hits and repressive_marks :
            return "The promoter region simultaneously shows both open and inhibitory signals, suggesting it may be subject to tissue-specific or dynamic epigenetic regulation."
        if poised_marks and not chromatin_hits and not active_marks :
            return "The area near the promoter is dominated by poised/gene-body-related modifications, suggesting this site may be in a ready-to-activate or transcription-associated state."
        if chromatin_hits or histone_hits :
            return "An epigenomic peak is detected in the promoter region, but the directional interpretation across different data types requires integration with TissueExpression and co-expression results."
        return "No promoter peak was observed in the current data, which may indicate that this site is inactive in the included tissue, or that the peak set does not match this gene."

    def _safe_float (self ,value :str )->float :
        try :
            return float (value )
        except (TypeError ,ValueError ):
            return float ("-inf")

    def _normalize_tissue (self ,tissue :str )->str :
        value =tissue .strip ().lower ().replace ("-"," ")
        aliases ={
        "panicles":"panicles",
        "roots":"roots",
        "leaves":"leaves",
        "mature leaves":"mature_leaves",
        "seedlings":"seedlings",
        "sam":"sam",
        "young leaves":"seedlings",
        }
        return aliases .get (value ,value .replace (" ","_"))
