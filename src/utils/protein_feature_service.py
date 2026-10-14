from __future__ import annotations 

from io import BytesIO 
from typing import Any 

import matplotlib 
matplotlib .use ("Agg")
import matplotlib .pyplot as plt 
import matplotlib .patches as patches 
import requests 


class ProteinFeatureService :
    def __init__ (self ,timeout :int =15 ):
        self .timeout =timeout 
        self .session =requests .Session ()
        self .base_interpro_url ="https://www.ebi.ac.uk/interpro/api"
        self ._json_cache :dict [str ,dict [str ,Any ]]={}

    def collect (self ,gene_info :dict [str ,Any ])->dict [str ,Any ]:
        uniprot_id =(gene_info .get ("UniProt_ID")or "").strip ().upper ()
        result ={
        "rap_id":gene_info .get ("RAP_ID",""),
        "uniprot_id":uniprot_id ,
        "local_interpro":gene_info .get ("InterPro","N/A"),
        "local_subcellular":gene_info .get ("Subcellular","N/A"),
        "subcellular_locations":[],
        "interpro_domains":[],
        "pfam_domains":[],
        "protein_length":None ,
        "plot_png":None ,
        "warnings":[],
        }

        if not uniprot_id :
            result ["warnings"].append ("This gene lacks a UniProt ID, so protein domain and subcellular localization cannot be queried online.")
            return result 

        length =self ._safe_call (self .get_length ,uniprot_id ,result ["warnings"])
        subcell =self ._safe_call (self .fetch_subcellular_locations ,uniprot_id ,result ["warnings"],default =[])
        interpro =self ._safe_call (self .fetch_interpro_domains ,uniprot_id ,result ["warnings"],default =[])
        pfam =self ._safe_call (self .fetch_pfam_domains ,uniprot_id ,result ["warnings"],default =[])

        result ["protein_length"]=length 
        result ["subcellular_locations"]=subcell 
        result ["interpro_domains"]=interpro 
        result ["pfam_domains"]=pfam 

        if length and pfam :
            result ["plot_png"]=self ._safe_call (
            self .build_domain_plot ,
            uniprot_id ,
            length ,
            pfam ,
            result ["warnings"],
            default =None ,
            )
        else :
            result ["warnings"].append ("Insufficient Pfam or protein length information; domain diagram not generated.")

        return result 

    def get_length (self ,uniprot_id :str )->int |None :
        data =self ._get_json (f"https://rest.uniprot.org/uniprotkb/{uniprot_id }.json")
        return data .get ("sequence",{}).get ("length")

    def fetch_subcellular_locations (self ,uniprot_id :str )->list [str ]:
        data =self ._get_json (f"https://rest.uniprot.org/uniprotkb/{uniprot_id }.json")
        locations =[]
        for comment in data .get ("comments",[]):
            if comment .get ("commentType")!="SUBCELLULAR LOCATION":
                continue 
            for subcell in comment .get ("subcellularLocations",[]):
                loc_name =subcell .get ("location",{}).get ("value")
                if loc_name and loc_name not in locations :
                    locations .append (loc_name )
        return locations 

    def fetch_interpro_domains (self ,uniprot_id :str )->list [dict [str ,Any ]]:
        return self ._fetch_interpro_entries (
        f"{self .base_interpro_url }/entry/interpro/protein/uniprot/{uniprot_id }/",
        uniprot_id ,
        )

    def fetch_pfam_domains (self ,uniprot_id :str )->list [dict [str ,Any ]]:
        return self ._fetch_interpro_entries (
        f"{self .base_interpro_url }/entry/pfam/protein/uniprot/{uniprot_id }/",
        uniprot_id ,
        )

    def build_domain_plot (
    self ,
    uniprot_id :str ,
    length :int ,
    domains :list [dict [str ,Any ]],
    )->bytes |None :
        if not length or not domains :
            return None 

        fig_width =12 
        fig_height =max (3.2 ,len (domains )*0.65 )
        fig ,ax =plt .subplots (figsize =(fig_width ,fig_height ))

        backbone_y =0.5 
        backbone_height =0.08 
        rect_backbone =patches .Rectangle (
        (1 ,backbone_y -backbone_height /2 ),
        length ,
        backbone_height ,
        facecolor ="#f0f0f0",
        edgecolor ="#cccccc",
        linewidth =1 ,
        alpha =0.9 ,
        zorder =1 ,
        )
        ax .add_patch (rect_backbone )

        colors =["#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f","#edc948"]
        domain_height =0.22 
        label_y_step =0.12 

        for idx ,domain in enumerate (domains ):
            start =domain ["Start"]
            end =domain ["End"]
            width =end -start +1 
            center =start +width /2 
            label_y =backbone_y +0.16 +(idx %2 )*label_y_step 

            rect_domain =patches .Rectangle (
            (start ,backbone_y -domain_height /2 ),
            width ,
            domain_height ,
            facecolor =colors [idx %len (colors )],
            edgecolor ="black",
            linewidth =1.2 ,
            zorder =3 ,
            )
            ax .add_patch (rect_domain )

            label =f"{domain ['Name']}\n({domain ['Accession']})\n[{start }-{end } aa]"
            ax .text (
            center ,
            label_y ,
            label ,
            ha ="center",
            va ="center",
            fontsize =9 ,
            fontweight ="bold",
            bbox =dict (facecolor ="white",edgecolor ="none",alpha =0.7 ,pad =1 ),
            zorder =4 ,
            )
            ax .plot (
            [center ,center ],
            [backbone_y ,label_y ],
            color ="grey",
            linestyle ="--",
            linewidth =0.5 ,
            alpha =0.5 ,
            zorder =2 ,
            )

        ax .set_xlim (-(length *0.02 ),length *1.02 )
        ax .set_ylim (0.2 ,backbone_y +(label_y_step *2 )+0.3 )
        ax .set_xlabel ("Amino Acid Position (aa)",fontsize =11 ,fontweight ="bold",labelpad =10 )
        ax .tick_params (axis ="x",labelsize =10 )
        ax .set_yticks ([])
        ax .spines ["left"].set_visible (False )
        ax .spines ["right"].set_visible (False )
        ax .spines ["top"].set_visible (False )
        plt .title (
        f"Protein Architecture: {uniprot_id } (Based on Pfam Domains)",
        fontsize =14 ,
        fontweight ="bold",
        pad =25 ,
        )
        plt .tight_layout ()

        output =BytesIO ()
        fig .savefig (output ,format ="png",dpi =220 ,bbox_inches ="tight")
        plt .close (fig )
        output .seek (0 )
        return output .getvalue ()

    def _fetch_interpro_entries (self ,url :str ,uniprot_id :str )->list [dict [str ,Any ]]:
        results =[]
        next_url =url 
        headers ={"Accept":"application/json"}

        while next_url :
            data =self ._get_json (next_url ,headers =headers )
            for result in data .get ("results",[]):
                domain_acc =result .get ("metadata",{}).get ("accession","Unknown")
                domain_name =result .get ("metadata",{}).get ("name","Unknown")
                for entry_protein in result .get ("proteins",[]):
                    if entry_protein .get ("accession","").upper ()!=uniprot_id .upper ():
                        continue 
                    for location in entry_protein .get ("entry_protein_locations",[]):
                        for fragment in location .get ("fragments",[]):
                            start =fragment .get ("start")
                            end =fragment .get ("end")
                            if start is None or end is None :
                                continue 
                            results .append (
                            {
                            "Accession":domain_acc ,
                            "Name":domain_name ,
                            "Start":start ,
                            "End":end ,
                            }
                            )
            next_url =data .get ("next")

        return sorted (results ,key =lambda item :(item ["Start"],item ["End"],item ["Accession"]))

    def _get_json (self ,url :str ,headers :dict [str ,str ]|None =None )->dict [str ,Any ]:
        cache_key =url 
        if cache_key in self ._json_cache :
            return self ._json_cache [cache_key ]
        response =self .session .get (url ,headers =headers ,timeout =self .timeout )
        response .raise_for_status ()
        data =response .json ()
        self ._json_cache [cache_key ]=data 
        return data 

    @staticmethod 
    def _safe_call (func ,*args ,default =None ):
        warnings =args [-1 ]
        call_args =args [:-1 ]
        try :
            return func (*call_args )
        except Exception as exc :
            warnings .append (f"{func .__name__ } queryFailed: {exc }")
            return default 
