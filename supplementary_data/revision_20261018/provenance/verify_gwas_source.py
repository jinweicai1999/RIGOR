"""Bounded source-identity audit, not a new RIGOR benchmark or full GWAS rerun.

Usage: python verify_gwas_source.py --local-gwas-dir /path/to/data/processed
Public inputs: RiceVarMap_phenos.csv, RiceVarMap_genotypes_1.csv alongside script.
No API credentials or private data are sent to any service.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess
import numpy as np

def beta_fraction(a,b,x):
    qab=a+b;qap=a+1;qam=a-1;c=1.;d=1.-qab*x/qap
    tiny=1e-300
    if abs(d)<tiny:d=tiny
    d=1/d;h=d
    for m in range(1,1001):
        aa=m*(b-m)*x/((qam+2*m)*(a+2*m));d=1+aa*d
        if abs(d)<tiny:d=tiny
        c=1+aa/c
        if abs(c)<tiny:c=tiny
        d=1/d;h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+2*m)*(qap+2*m));d=1+aa*d
        if abs(d)<tiny:d=tiny
        c=1+aa/c
        if abs(c)<tiny:c=tiny
        d=1/d;delta=d*c;h*=delta
        if abs(delta-1)<2e-14:return h
    raise RuntimeError('Incomplete beta did not converge')

def ibeta(a,b,x):
    if x==0:return 0.
    if x==1:return 1.
    bt=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
    if x<(a+1)/(a+b+2):return bt*beta_fraction(a,b,x)/a
    return 1-bt*beta_fraction(b,a,1-x)/b

def t_pvalue(t,df):return ibeta(df/2,.5,df/(df+t*t))

def run(local_dir,folder):
    phrows=list(csv.DictReader((folder/'RiceVarMap_phenos.csv').open()))
    phenos={r['id_name']:r for r in phrows}
    rows=list(csv.reader((folder/'RiceVarMap_genotypes_1.csv').open()))
    variants=['vg0100001151','vg0100001184','vg0311647322','vg0507081467']
    # The fifth queried marker has deletion/multi-allelic portal calls; it is not
    # silently recoded. Four biallelic markers are used for the identity test.
    traits=list(phrows[0])[2:]
    files=[local_dir/('gwas_'+t+'.assoc.linear') for t in traits]
    args=['rg','-m','4','-N','--with-filename']
    for v in variants:args += ['-e',v]
    args += [str(p) for p in files]
    found=subprocess.run(args,capture_output=True,text=True,check=True).stdout
    local={}
    for line in found.splitlines():
        f,s=line.split(':',1);fields=s.split();trait=Path(f).name[len('gwas_'):-len('.assoc.linear')]
        chr_,v,bp,allele,test,n,b,t,p=fields
        local[(trait,v)]=dict(allele=allele,n=int(n),beta=b,t=t,p=p)
    assert len(local)==40
    comparisons=[]
    for trait in traits:
        for variant in variants:
            j=rows[0].index(variant);record=local[(trait,variant)];allele=record['allele']
            major,minor=rows[5][j],rows[6][j]
            assert allele in {major,minor}
            other=major if allele==minor else minor
            points=[]
            for row in rows[7:]:
                if not row[0].endswith('Imputed'):continue
                accession=row[0].split(' ',1)[0];value=float(phenos[accession][trait]);call=row[j]
                if value==-9 or not math.isfinite(value):continue
                if call==allele:dosage=2
                elif call==other:dosage=0
                elif call=='N':continue
                else:raise ValueError(('Unspecified genotype',variant,call))
                points.append((dosage,value))
            a=np.asarray(points,dtype=float);x,y=a.T;n=len(y);x=x-x.mean();y=y-y.mean()
            beta=float((x*y).sum()/(x*x).sum())
            se=float(np.sqrt(((y-beta*x)**2).sum()/(n-2)/(x*x).sum()))
            t=beta/se;p=t_pvalue(t,n-2)
            matches=[n==record['n'],float(format(beta,'.4g'))==float(record['beta']),float(format(t,'.4g'))==float(record['t']),float(format(p,'.4g'))==float(record['p'])]
            comparisons.append(dict(trait=trait,variant=variant,effect_allele=allele,
                local_n=record['n'],reconstructed_n=n,local_beta=record['beta'],reconstructed_beta=format(beta,'.16g'),
                local_stat=record['t'],reconstructed_stat=format(t,'.16g'),local_p=record['p'],reconstructed_p=format(p,'.16g'),
                n_match=matches[0],beta_match_4_significant_digits=matches[1],stat_match_4_significant_digits=matches[2],p_match_4_significant_digits=matches[3],all_match=all(matches)))
    with (folder/'GWAS_Source_Reconstruction_40_checks_1019.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    result=dict(verification_date='2026-10-19',purpose='Source/population identification using a bounded numerical fingerprint; not a full GWAS reproduction',
        phenotype_source_url='https://ricevarmap.ncpgr.cn/media/Phenos/phenos.csv',genotype_source_url='https://ricevarmap.ncpgr.cn/v3/vars_genotype/',
        phenotype_source_sha256=hashlib.sha256((folder/'RiceVarMap_phenos.csv').read_bytes()).hexdigest(),
        genotype_table_sha256=hashlib.sha256((folder/'RiceVarMap_genotypes_1.csv').read_bytes()).hexdigest(),
        panel_size=529,study='Xie et al. (2015), DOI 10.1073/pnas.1515919112',
        population='All 529 accessions retained in the published 533-accession panel after exclusion of four accessions; not the combined 1479-accession panel or the 3K cohort',
        variants=variants,traits=traits,total_checks=len(comparisons),all_fields_matching_checks=sum(r['all_match'] for r in comparisons),
        fields_checked=['NMISS','BETA','STAT','P'],precision='Exact NMISS; BETA, STAT and P matched to the 4-significant-digit precision retained in the local PLINK output',
        reconstruction='OLS with an intercept and additive 0/2 dosage for homozygous effect-allele counts; nonmissing imputed genotypes; phenotype -9 excluded; two-sided t test with n-2 degrees of freedom; no covariates',
        phenotype_nonmissing_counts={t:sum(float(r[t])!=-9 for r in phrows) for t in traits},
        limitations=['Bounded fingerprint does not independently validate every SNP or recover the historical PLINK command/version.','Portal precomputed LR/LMM P-values are not the local unadjusted PLINK ADD output; a difference alone does not refute source-population identity.','One of five queried markers (vg0503927501) has deletion/multi-allelic display calls and was not recoded for reconstruction.'])
    (folder/'GWAS_Source_Reconstruction_1019.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps(result,ensure_ascii=False,indent=2))
    assert result['all_fields_matching_checks']==40,[r for r in comparisons if not r['all_match']]

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--local-gwas-dir',type=Path,required=True)
    run(p.parse_args().local_gwas_dir,Path(__file__).resolve().parent)
