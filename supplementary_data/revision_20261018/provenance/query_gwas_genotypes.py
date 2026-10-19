from pathlib import Path
import csv,hashlib,http.cookiejar,json,re,urllib.parse,urllib.request
from lxml import html

out=Path('/private/tmp/rigor_source_check1018')
url='https://ricevarmap.ncpgr.cn/v3/vars_genotype/'
jar=http.cookiejar.CookieJar()
opener=urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
page=opener.open(urllib.request.Request(url,headers={'User-Agent':'RIGOR source verification'}),timeout=40).read()
root=html.fromstring(page)
form=root.xpath('//form[@id="geno_form"]')[0]
token=form.xpath('.//input[@name="csrfmiddlewaretoken"]/@value')[0]
phenos=list(csv.DictReader((out/'RiceVarMap_phenos.csv').open()))
ids={x['id_name'] for x in phenos}
selected=[]
for option in form.xpath('.//select[@name="cultivar"]/option'):
 label=' '.join(option.itertext()).strip()
 m=re.search(r'\(([^()]*)\)$',label)
 if m and m.group(1) in ids:selected.append((m.group(1),option.get('value')))
assert len(selected)==len(ids),(len(selected),len(ids))
variants=['vg0100001151','vg0100001184','vg0311647322','vg0503927501','vg0507081467']
fields=[('csrfmiddlewaretoken',token),('var_ids',', '.join(variants)),('type','all')]
fields += [('cultivar',v) for _,v in selected]
req=urllib.request.Request(url,data=urllib.parse.urlencode(fields).encode(),headers={'User-Agent':'RIGOR source verification','Referer':url})
response=opener.open(req,timeout=55).read()
(out/'RiceVarMap_529_five_SNP_genotypes.html').write_bytes(response)
root=html.fromstring(response)
result={'url':url,'phenotype_n':len(phenos),'genotype_accessions_requested':len(selected),'variants':variants,'response_bytes':len(response),'response_sha256':hashlib.sha256(response).hexdigest(),'tables':[]}
for table in root.xpath('//table'):
 rows=[[' '.join(x.itertext()).strip() for x in row.xpath('./td|./th')] for row in table.xpath('.//tr')]
 result['tables'].append({'id':table.get('id'),'rows':len(rows),'first_rows':rows[:3]})
 if len(rows)>10:
  with (out/('RiceVarMap_genotypes_'+str(len(result['tables']))+'.csv')).open('w',newline='') as f:csv.writer(f).writerows(rows)
(out/'Genotype_query_record.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(result,ensure_ascii=False,indent=2))
