#!/usr/bin/env python3
import io
import re
import tokenize
from pathlib import Path
from collections import OrderedDict
from dotenv import dotenv_values
from openai import OpenAI

ROOT=Path(__file__).resolve().parents[1]
HAN_RE=re.compile(r'[\u4e00-\u9fff]')

TARGET_EXT={'.py','.md','.txt'}
EXCLUDE_DIRS={'.venv','.git','data','figures'}


def has_han(s):
    return bool(HAN_RE.search(s or ''))


def split_string_token(tok:str):
    m=re.match(r"(?is)^([rubf]*)('''|\"\"\"|'|\")(.*)(\2)$",tok)
    return (m.group(1),m.group(2),m.group(3),m.group(4)) if m else None


def iter_files():
    for p in ROOT.rglob('*'):
        if not p.is_file():
            continue
        if any(part in EXCLUDE_DIRS for part in p.parts):
            continue
        if p.suffix.lower() not in TARGET_EXT:
            continue
        yield p


def collect_segments(path):
    txt=path.read_text(encoding='utf-8',errors='ignore')
    segs=[]
    if path.suffix.lower()=='.py':
        toks=list(tokenize.generate_tokens(io.StringIO(txt).readline))
        for ttype,tstr,*_ in toks:
            if ttype==tokenize.COMMENT and has_han(tstr):
                segs.append(('comment',tstr[1:]))
            elif ttype==tokenize.STRING and has_han(tstr):
                sp=split_string_token(tstr)
                if sp:
                    segs.append(('string',sp[2]))
    else:
        for line in txt.splitlines():
            if has_han(line):
                segs.append(('line',line))
    return segs


def translate_batch(client, items):
    prompt_lines=[]
    for i,s in enumerate(items):
        prompt_lines.append(f"{i}\t{s}")
    payload='\n'.join(prompt_lines)
    prompt=(
        "Translate Chinese text into natural English. Preserve placeholders like {x}, markdown symbols, and code tokens.\n"
        "Return exactly the same number of lines in the format: <index>\\t<translation>. No extra text.\n\n"
        f"{payload}"
    )
    try:
        rsp=client.chat.completions.create(model='LongCat-Flash-Lite',messages=[{'role':'user','content':prompt}],temperature=0)
        out=(rsp.choices[0].message.content or '').strip().splitlines()
    except Exception:
        return {s:s for s in items}

    mapping={s:s for s in items}
    for line in out:
        if '\t' not in line:
            continue
        idx,trans=line.split('\t',1)
        idx=idx.strip()
        if not idx.isdigit():
            continue
        i=int(idx)
        if 0<=i<len(items):
            mapping[items[i]]=trans.strip() or items[i]
    return mapping


def apply_file(path, tr):
    txt=path.read_text(encoding='utf-8',errors='ignore')
    changed=False
    if path.suffix.lower()=='.py':
        toks=list(tokenize.generate_tokens(io.StringIO(txt).readline))
        new=[]
        for ttype,tstr,*rest in toks:
            n=tstr
            if ttype==tokenize.COMMENT and has_han(tstr):
                body=tstr[1:]
                n='#'+tr.get(body,body)
            elif ttype==tokenize.STRING and has_han(tstr):
                sp=split_string_token(tstr)
                if sp:
                    p,q1,b,q2=sp
                    n=f"{p}{q1}{tr.get(b,b)}{q2}"
            if n!=tstr:
                changed=True
            new.append((ttype,n))
        if changed:
            path.write_text(tokenize.untokenize(new),encoding='utf-8')
        return changed

    lines=txt.splitlines(keepends=True)
    out=[]
    for ln in lines:
        body=ln[:-1] if ln.endswith('\n') else ln
        if has_han(body):
            nb=tr.get(body,body)
            if nb!=body:
                changed=True
            out.append(nb + ('\n' if ln.endswith('\n') else ''))
        else:
            out.append(ln)
    if changed:
        path.write_text(''.join(out),encoding='utf-8')
    return changed


def main():
    cfg=dotenv_values(ROOT/'.env')
    client=OpenAI(api_key=cfg.get('LLM_API_KEY'), base_url=cfg.get('LLM_BASE_URL'), timeout=60.0)

    segs=OrderedDict()
    files=list(iter_files())
    for p in files:
        for _,s in collect_segments(p):
            if has_han(s):
                segs.setdefault(s,None)

    uniq=list(segs.keys())
    print('remaining segments',len(uniq), flush=True)
    if not uniq:
        return

    tr={}
    bs=20
    total=(len(uniq)+bs-1)//bs
    for i in range(0,len(uniq),bs):
        chunk=uniq[i:i+bs]
        m=translate_batch(client,chunk)
        tr.update(m)
        print(f'batch {i//bs+1}/{total}', flush=True)

    changed=0
    for p in files:
        try:
            if apply_file(p,tr):
                changed+=1
        except Exception:
            pass
    print('changed files',changed, flush=True)

if __name__=='__main__':
    main()
