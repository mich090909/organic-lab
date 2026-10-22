"""Restore audited extension from the pinned source and recorded corrections.
Runtime changes must hash to the previously tested extension before publication.
"""
import re,json,collections,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
from chemistry import *
ROOT=Path(__file__).resolve().parents[1]
def details(m,idx,numbering=None,sides=None):
    mh,order=priorities(m,idx);nums={v:str(i+1) for i,v in enumerate(numbering or [])};out=[]
    for rank,j in enumerate(order,1):
        a=mh.GetAtomWithIdx(j);sym=a.GetSymbol();shell=[]
        for b in a.GetBonds():
            nb=b.GetOtherAtom(a);times=int(b.GetBondTypeAsDouble())-(1 if nb.GetIdx()==idx else 0);shell += [nb.GetSymbol()]*max(0,times)
        zs=Chem.GetPeriodicTable();shell.sort(key=lambda s:zs.GetAtomicNumber(s),reverse=True)
        group=('verso C'+nums[j] if j in nums else {'O':'O (gruppo ossigenato)','N':'N (gruppo amminico)'}.get(sym,sym))
        if shell:group+=' ['+','.join(shell)+']'
        row={'priority':rank,'group':group,'atomIndex':j,'firstShell':shell}
        if sides:row['side']=sides[j]
        out.append(row)
    ans=label(m)[idx];d={'answer':ans,'priorities':out,'engine':'RDKit 2026.09.1 modern CIPLabeler','view':'Fischer' if sides else 'Scheletrica'}
    if sides:d['group4']=out[-1]['side'];d['geometryAnswer']=geometric_rs([v['side'] for v in out]);d['orientation']='orizzontali davanti; verticali dietro'
    return d
def fix_fischer(mol):
    ref,nm,match,sts=name_stereo_check(mol);n=len(parse_name(mol['name'])[2]);h=116+74*(n-2)
    pieces=[f'<svg class="fischer-figure" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 360 {h}" width="360" height="{h}" role="img" aria-label="Proiezione di Fischer"><rect width="360" height="100%" fill="white"/>']
    def text(x,y,t):return f'<text x="{x}" y="{y}" text-anchor="middle" font-family="Arial,sans-serif" font-size="20" fill="#173347">{t}</text>'
    def line(x1,y1,x2,y2):return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#233c4a" stroke-width="2.7"/>'
    pieces += [line(180,51,180,h-45),text(180,35,'COOH' if mol['name'].startswith('acido') else 'CH₃'),text(180,h-19,'CH₃')]
    for i in range(1,n-1):
        y=29+74*i;outside=[a for a in nm.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx()>=n]
        if len(outside)>1:raise ValueError('unsupported Fischer branch')
        group={'C':'CH₃','O':'OH','N':'NH₂','F':'F','Cl':'Cl','Br':'Br'}[outside[0].GetSymbol()] if outside else 'H'
        pieces += [line(94,y,266,y),text(67,y+6,'H'),text(293,y+6,group)]
        if group=='H':pieces += ['<rect x="168" y="'+str(y-16)+'" width="24" height="32" fill="white"/>',text(180,y+6,'C')]
        else:pieces += [f'<text x="12" y="{y+5}" fill="#617084" font-size="14">C{i+1}</text>']
    svg=''.join(pieces)+'</svg>';fm,oldrs=read_fischer(svg);target={str(p):v for p,v in sts if v in 'RS'};root=ET.fromstring(svg);ns='{http://www.w3.org/2000/svg}'
    rows=sorted(set(float(t.get('y')) for t in root.iter(ns+'text') if t.get('x') in ['67','293']))
    for key,d in oldrs.items():
        if d['answer']!=target[key]:
            y=rows[int(key)-2];left=next(t for t in root.iter(ns+'text') if t.get('x')=='67' and float(t.get('y'))==y);right=next(t for t in root.iter(ns+'text') if t.get('x')=='293' and float(t.get('y'))==y);left.text,right.text=right.text,left.text
    ET.register_namespace('',ns[1:-1]);svg=ET.tostring(root,encoding='unicode');fm,rs=read_fischer(svg);assert canonical(fm)==canonical(ref);meta={}
    for key,d in rs.items():
        idx=d['atom'];mh,order=priorities(fm,idx);sd=dict(zip(order,d['sides']));meta[key]=details(fm,idx,list(range(len(rows)+2)),sd)
    return svg,meta
def repair():
    import subprocess
    ext=json.loads(subprocess.check_output(['git','show','640c01c42d4666cc01ecaa2d660796c0af80639b:organic-v25-extension.js'],cwd=ROOT,text=True).split('=',1)[1].strip().rstrip(';'))
    findings=json.loads((ROOT/'audit/changes.json').read_text());namefix={c['id']:c['after'] for c in findings if c['type'] in ['numbering','name punctuation']};whyfix={c['id']:c['after'] for c in findings if c['type']=='explanation precision'}
    for m in ext['molecules']:
        old=m['name'];m['name']=namefix.get(m['id'],m['name']);g=ext['guide'][str(m['id'])]
        if m['name']!=old:g['stereo']=[{'position':p,'value':v} for p,v in parse_name(m['name'])[1]]
        if m.get('fischerSvg'):m['fischerSvg'],m['cipFischer']=fix_fischer(m)
        g['prefix']=prefix_for_name(m['name']);locs=re.findall(r'\d+',re.sub(r'\([^)]*\)-','',m['name']));g['locants']=','.join(map(str,sorted(map(int,locs)))) or 'nessuno'
        ref=Chem.MolFromSmiles(m['smiles']);pot=Chem.FindPotentialStereo(ref);m['stereoCenters']=sum(x.type==Chem.StereoType.Atom_Tetrahedral for x in pot);m['unspecifiedStereo']=any(x.specified==Chem.StereoSpecified.Unspecified for x in pot);m['stereoDescriptors']=[p+v for p,v in parse_name(m['name'])[1]]
        if m['id'] in whyfix:m['why']=whyfix[m['id']]
        m['auditStatus']='checked'
    byid={m['id']:m for m in ext['molecules']}
    for q in ext['extraRS']:
        m=byid[q['moleculeId']];q['name']=m['name'];q['centers']=[{'position':int(p),'answer':v} for p,v in parse_name(m['name'])[1] if v in 'RS'];q['answer']=','.join(str(c['position'])+c['answer'] for c in q['centers']);ref,nm,match,sts=name_stereo_check(m)
        q['cipDetails']={p:details(ref,match[int(p)-1],[match[i] for i in range(len(parse_name(m['name'])[2]))]) for p,v in sts if v in 'RS'}
        q['why']='Confronta i sostituenti secondo CIP fino al primo punto di differenza; i legami multipli si espandono con atomi duplicati. Ricava R/S dai legami effettivamente disegnati.';q['auditStatus']='checked';q['stereoCenters']=m['stereoCenters']
    ext['legacyFixes']={};ext['legacyGuide']={};ext['existingRS']={};ext['release']='2.5.1'
    out=('window.ORG_LAB_V25='+json.dumps(ext,ensure_ascii=False,separators=(',',':'))+';\n').encode();sha=hashlib.sha1(b'blob '+str(len(out)).encode()+b'\0'+out).hexdigest();print('restored extension blob',sha,'bytes',len(out));assert sha=='51cd906e8c91db72450be48f54502272ece8abf1','Restoration differs from tested artifact';(ROOT/'organic-v25-extension.js').write_bytes(out)
if __name__=='__main__':repair()
