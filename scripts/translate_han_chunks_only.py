#!/usr/bin/env python3
import re
from pathlib import Path
from dotenv import dotenv_values
from openai import OpenAI

ROOT=Path(__file__).resolve().parents[1]
HAN_CHUNK_RE=re.compile(r'[\u4e00-\u9fff]+')
SKIP_DIRS={'.venv','.git','data','figures'}
TEXT_EXT={'.py','.md','.txt','.toml','.yaml','.yml','.json','.csv'}


def iter_files():
    for p in ROOT.rglob('*'):
        if not p.is_file():
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if p.suffix.lower() not in TEXT_EXT:
            continue
        yield p


def collect_chunks(files):
    chunks=set()
    for p in files:
        txt=p.read_text(encoding='utf-8',errors='ignore')
        chunks.update(HAN_CHUNK_RE.findall(txt))
    return sorted(chunks, key=len, reverse=True)


def translate_chunks(chunks):
    cfg=dotenv_values(ROOT/'.env')
    client=OpenAI(api_key=cfg.get('LLM_API_KEY'), base_url=cfg.get('LLM_BASE_URL'), timeout=60.0)
    out={}
    bs=40
    total=(len(chunks)+bs-1)//bs
    for i in range(0,len(chunks),bs):
        sub=chunks[i:i+bs]
        lines='\n'.join([f"{j}\t{s}" for j,s in enumerate(sub)])
        prompt=(
            "Translate each Chinese term/phrase to concise English. "
            "Return exactly same number of lines as: <index>\\t<translation>. No extra text.\n\n"
            + lines
        )
        try:
            rsp=client.chat.completions.create(model='LongCat-Flash-Lite',messages=[{'role':'user','content':prompt}],temperature=0)
            txt=(rsp.choices[0].message.content or '').strip().splitlines()
            local={}
            for ln in txt:
                if '\t' not in ln:
                    continue
                k,v=ln.split('\t',1)
                if k.strip().isdigit():
                    idx=int(k.strip())
                    if 0<=idx<len(sub):
                        local[sub[idx]]=v.strip() or sub[idx]
            for s in sub:
                out[s]=local.get(s,s)
        except Exception:
            for s in sub:
                out[s]=s
        print(f'batch {i//bs+1}/{total}', flush=True)
    return out


def apply(files, mapping):
    changed=0
    for p in files:
        txt=p.read_text(encoding='utf-8',errors='ignore')
        new=txt
        for zh,en in mapping.items():
            if zh in new:
                new=new.replace(zh,en)
        if new!=txt:
            p.write_text(new,encoding='utf-8')
            changed+=1
    return changed


def main():
    files=list(iter_files())
    chunks=collect_chunks(files)
    print('han chunks',len(chunks), flush=True)
    if not chunks:
        return
    mapping=translate_chunks(chunks)
    c=apply(files,mapping)
    print('changed files',c, flush=True)

if __name__=='__main__':
    main()
