"""Independent bounded graph/name/SVG checks for Organic Lab.
Unsupported syntax raises; saved answers are never geometry inputs.
"""
import re,collections,xml.etree.ElementTree as ET
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler,rdMolDescriptors
ROOTS={'met':1,'et':2,'prop':3,'but':4,'pent':5,'es':6,'ept':7,'ott':8,'non':9,'dec':10,'undec':11,'dodec':12,'tridec':13}
BT={1:Chem.BondType.SINGLE,2:Chem.BondType.DOUBLE,3:Chem.BondType.TRIPLE}
PARENT_RE=r'(ciclo)?('+'|'.join(sorted(ROOTS,key=len,reverse=True))+r')(an(?:o|ale|one|olo|ammina|ammide|nitrile|-)|ano(?:ico|ato)|(?:-\d+(?:,\d+)*-)?(?:ene|ino)|ene|ino|anale|anammide|anonitrile)'
def canonical(m,stereo=True):return Chem.MolToSmiles(m,isomericSmiles=stereo)
def label(m):
    for a in m.GetAtoms():
        if a.HasProp('_CIPCode'):a.ClearProp('_CIPCode')
    rdCIPLabeler.AssignCIPLabels(m)
    return {a.GetIdx():a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode')}
def prefix_for_name(name):
    base=re.sub(r'\([^)]*\)-','',name).replace('acido ','').split(' di ')[0]
    for e in ['benzene','fenolo','benzenammina','benzoico','benzoato','benzonitrile']:
        if base.endswith(e):return base[:-len(e)].rstrip('-') or 'nessuno'
    ms=list(re.finditer(PARENT_RE,base))
    if not ms:raise ValueError('unsupported parent '+base)
    return base[:ms[-1].start()].rstrip('-') or 'nessuno'
def parse_name(name):
    name=name.replace('acido ','');stereo=[];mt=re.match(r'^\(([^)]+)\)-',name)
    if mt:stereo=re.findall(r'(\d*)([RSEZ])',mt[1]);name=name[mt.end():]
    ester=None
    if ' di ' in name:name,ester=name.split(' di ')
    arom=False;ring=False;suffix='';prefix='';n=0
    ends={'benzene':('aromatic',6,''),'fenolo':('aromatic',6,'olo'),'benzenammina':('aromatic',6,'ammina'),'benzoico':('aromatic',6,'oico'),'benzoato':('aromatic',6,'oato'),'benzonitrile':('aromatic',6,'nitrile')}
    for ending,(kind,num,suf) in ends.items():
        if name.endswith(ending):prefix=name[:-len(ending)];n=num;arom=True;ring=True;suffix=suf;break
    else:
        ms=list(re.finditer(PARENT_RE,name))
        if not ms:raise ValueError('unsupported name '+name)
        mt=ms[-1];prefix=name[:mt.start()];ring=bool(mt[1]);n=ROOTS[mt[2]];suffix=name[mt.start()+len(mt[1] or '')+len(mt[2]):]
    rw=Chem.RWMol()
    def atom(sym):return rw.AddAtom(Chem.Atom(sym))
    def attach(pos,sym,order=1):
        j=atom(sym);rw.AddBond(pos,j,BT[order]);return j
    for _ in range(n):atom('C')
    for i in range(n-1):rw.AddBond(i,i+1,BT[1])
    if ring:rw.AddBond(n-1,0,BT[1])
    if arom:
        for a in rw.GetAtoms():a.SetIsAromatic(True);a.SetNoImplicit(False)
        for b in rw.GetBonds():b.SetBondType(Chem.BondType.AROMATIC)
    specs=[]
    if suffix in ['ano','']:pass
    elif suffix in ['ene','ino']:specs=[(0,'bond',2 if suffix=='ene' else 3)]
    elif re.fullmatch(r'-\d+-(ene|ino)',suffix):
        p,t=re.fullmatch(r'-(\d+)-(ene|ino)',suffix).groups();specs=[(int(p)-1,'bond',2 if t=='ene' else 3)]
    elif suffix in ['anale','anone','anolo','anammina','anammide','anonitrile','anoico','anoato','olo','ammina','oico','oato','nitrile']:
        t={'anale':'ale','anone':'one','anolo':'olo','anammina':'ammina','anammide':'ammide','anonitrile':'nitrile','anoico':'oico','anoato':'oato'}.get(suffix,suffix);specs=[(1 if t=='one' else 0,t,1)]
    else:
        mt=re.fullmatch(r'ano?-(\d+(?:,\d+)*)-(di|tri)?(olo|one|ammina)',suffix)
        if not mt:raise ValueError('unsupported suffix '+suffix)
        specs=[(int(p)-1,mt[3],1) for p in mt[1].split(',')]
    functional=None
    for pos,t,order in specs:
        if arom and t in ['oico','oato']:pos=attach(pos,'C')
        if t=='bond':rw.GetBondBetweenAtoms(pos,pos+1).SetBondType(BT[order]);continue
        if t=='olo':functional=attach(pos,'O')
        elif t in ['ale','one']:functional=attach(pos,'O',2)
        elif t in ['oico','oato','ammide']:attach(pos,'O',2);functional=attach(pos,'N' if t=='ammide' else 'O')
        elif t=='ammina':functional=attach(pos,'N')
        elif t=='nitrile':
            if arom:c=attach(pos,'C');functional=attach(c,'N',3)
            else:functional=attach(pos,'N',3)
        else:raise ValueError(t)
    if ester:
        prev=functional
        for _ in range({'metile':1,'etile':2,'propile':3}[ester]):prev=attach(prev,'C')
    prefix=prefix.rstrip('-');cursor=0
    pat=r'(?:(N(?:,N)*|\d+(?:,\d+)*)-)?(tri|di)?(metossi|etossi|idrossi|ammino|cloro|bromo|fluoro|osso|fenil|metil|etil)'
    for mt in re.finditer(pat,prefix):
        if prefix[cursor:mt.start()].strip('-'):raise ValueError('unsupported prefix '+prefix)
        loc,mult,group=mt.groups();positions=loc.split(',') if loc else ['1']
        if mult and len(positions)!={'di':2,'tri':3}[mult]:raise ValueError('quantity mismatch')
        for p in positions:
            pos=functional if p=='N' else int(p)-1
            if pos is None or pos>=n and p!='N':raise ValueError('locant '+prefix)
            if group in ['metil','etil']:
                prev=pos
                for _ in range(1 if group=='metil' else 2):prev=attach(prev,'C')
            elif group in ['metossi','etossi']:
                prev=attach(pos,'O')
                for _ in range(1 if group=='metossi' else 2):prev=attach(prev,'C')
            elif group=='fenil':
                ids=[atom('C') for _ in range(6)];rw.AddBond(pos,ids[0],BT[1])
                for i in range(6):rw.GetAtomWithIdx(ids[i]).SetIsAromatic(True);rw.AddBond(ids[i],ids[(i+1)%6],Chem.BondType.AROMATIC)
            else:attach(pos,{'idrossi':'O','ammino':'N','cloro':'Cl','bromo':'Br','fluoro':'F','osso':'O'}[group],2 if group=='osso' else 1)
        cursor=mt.end()
    if prefix[cursor:].strip('-'):raise ValueError('unsupported prefix '+prefix)
    m=rw.GetMol();Chem.SanitizeMol(m)
    if arom:Chem.Kekulize(m,clearAromaticFlags=True);Chem.SetAromaticity(m)
    return m,stereo,list(range(n))
def priorities(m,idx):
    mh=Chem.AddHs(m);label(mh);a=mh.GetAtomWithIdx(idx)
    if not a.HasProp('_CIPNeighborOrder'):a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW);label(mh)
    if not a.HasProp('_CIPNeighborOrder'):raise ValueError('CIP unsupported/non-stereogenic')
    order=list(a.GetPropsAsDict()['_CIPNeighborOrder'])
    if len(order)!=4:raise ValueError('CIP incomplete')
    return mh,order
FCOORD={'top':(0,1,-1),'bottom':(0,-1,-1),'left':(-1,0,1),'right':(1,0,1)}
def geometric_rs(sides):
    p=np.array([FCOORD[x] for x in sides]);det=np.linalg.det(p[:3]-p[3])
    if abs(det)<1e-8:raise ValueError('degenerate geometry')
    return 'R' if det<0 else 'S'
def read_fischer(svg):
    root=ET.fromstring(svg);texts=[(float(t.get('x')),float(t.get('y')),t.text) for t in root.iter() if t.tag.endswith('text')]
    ends=sorted((y,t) for x,y,t in texts if x==180 and t!='C')
    if len(ends)!=2:raise ValueError('unsupported Fischer ends')
    rows=sorted(set(y-6 for x,y,t in texts if x in [67,293]));n=len(rows)+2;rw=Chem.RWMol()
    for _ in range(n):rw.AddAtom(Chem.Atom('C'))
    for i in range(n-1):rw.AddBond(i,i+1,BT[1])
    def attach(pos,sym,order=1):
        j=rw.AddAtom(Chem.Atom(sym));rw.AddBond(pos,j,BT[order]);return j
    for pos,(_,txt) in zip([0,n-1],ends):
        if txt=='COOH':attach(pos,'O',2);attach(pos,'O')
        elif txt!='CH₃':raise ValueError('unsupported Fischer end '+txt)
    sideatoms={}
    for i,y in enumerate(rows,1):
        sides={'top':i-1,'bottom':i+1}
        for x,side in [(67,'left'),(293,'right')]:
            txt=next(t for xx,yy,t in texts if xx==x and yy==y+6);sym={'H':'H','F':'F','Cl':'Cl','Br':'Br','OH':'O','NH₂':'N','CH₃':'C'}.get(txt)
            if not sym:raise ValueError('unsupported substituent '+txt)
            sides[side]=attach(i,sym)
        sideatoms[i]=sides
    m=rw.GetMol();Chem.SanitizeMol(m);results={}
    for idx,sides in sideatoms.items():
        if sum(m.GetAtomWithIdx(a).GetSymbol()=='H' for a in sides.values())==2:continue
        m.GetAtomWithIdx(idx).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW);mh,order=priorities(m,idx);inv={v:k for k,v in sides.items()};ordered=[inv[a] for a in order];answer=geometric_rs(ordered)
        results[str(idx+1)]={'answer':answer,'sides':ordered,'atom':idx}
        if label(m)[idx]!=answer:m.GetAtomWithIdx(idx).InvertChirality()
    m=Chem.RemoveHs(m);label(m);return m,results
def name_stereo_check(d):
    ref=Chem.MolFromSmiles(d['smiles']);codes=label(ref);nm,sts,chain=parse_name(d['name'])
    assert canonical(ref,False)==canonical(nm,False),'name graph differs'
    for match in ref.GetSubstructMatches(nm,uniquify=False):
        ok=True
        for p,v in sts:
            pos=int(p) if p else None
            if v in 'RS':val=codes.get(match[pos-1]) if pos else next(iter(codes.values())) if len(codes)==1 else '?'
            else:
                bs=[b for b in ref.GetBonds() if b.HasProp('_CIPCode') and b.GetProp('_CIPCode') in 'EZ']
                if pos:bs=[b for b in bs if {b.GetBeginAtomIdx(),b.GetEndAtomIdx()}=={match[pos-1],match[pos]}]
                val=bs[0].GetProp('_CIPCode') if len(bs)==1 else '?'
            ok &=val==v
        if ok:return ref,nm,match,sts
    raise AssertionError('name descriptors differ from structure')
COLORS={'#FF0000':'O','#0000FF':'N','#00CC00':'Cl','#33CCCC':'F','#7F4C19':'Br','#A01EEF':'I','#CCCC00':'S'}
GLYPHS={}
def learn_glyphs(svgs):
    candidates=collections.defaultdict(set)
    for svg in svgs:
        root=ET.fromstring(svg[svg.index('<svg'):])
        for el in root.iter():
            if re.fullmatch(r'atom-\d+',el.get('class','')):
                sym=COLORS.get(el.get('fill','').upper())
                if sym:candidates[''.join(re.findall('[MLQZ]',el.get('d','')))].add(sym)
    GLYPHS.update({k:next(iter(v)) for k,v in candidates.items() if len(v)==1})
def svg_graph(svg):
    root=ET.fromstring(svg[svg.index('<svg'):]);groups=collections.defaultdict(list);symbols={};coords={};wedges={}
    for el in root.iter():
        mt=re.fullmatch(r'atom-(\d+)',el.get('class',''))
        if mt:
            sym=COLORS.get(el.get('fill','').upper()) or GLYPHS.get(''.join(re.findall('[MLQZ]',el.get('d',''))))
            if sym:symbols[int(mt[1])]=sym
    for el in root.iter():
        mt=re.match(r'bond-(\d+) atom-(\d+) atom-(\d+)',el.get('class',''))
        if not mt:continue
        bi,a,b=map(int,mt.groups());style=el.get('style','');pts=np.array([(float(x),float(y)) for x,y in re.findall(r'(-?[\d.]+),(-?[\d.]+)',el.get('d',''))])
        if len(pts)<2:raise ValueError('unsupported SVG path')
        groups[bi].append((a,b,pts,style))
    for bi,paths in groups.items():
        a,b=paths[0][:2];triangles=[p for _,_,p,sty in paths if len(p)>=3 and 'fill:none' not in sty]
        if triangles:
            p=triangles[0];coords[a]=p[0];last=triangles[-1];coords[b]=np.mean(last[1:3],axis=0);wedges[bi]=(a,b,1)
        elif len(paths)>6:
            mids=[p.mean(axis=0) for _,_,p,_ in paths];width=[np.linalg.norm(p[-1]-p[0]) for _,_,p,_ in paths];order=np.argsort(width);i,j=int(order[0]),int(order[-1]);step=(mids[j]-mids[i])/(width[j]-width[i]);origin=mids[i]-step*width[i];coords[a]=origin;coords.setdefault(b,mids[j]);wedges[bi]=(a,b,-1)
        else:coords.setdefault(a,paths[0][2][0]);coords.setdefault(b,paths[-1][2][-1])
    rw=Chem.RWMol();crossed=[];num=max(max(p[0],p[1]) for ps in groups.values() for p in ps)+1
    for i in range(num):rw.AddAtom(Chem.Atom(symbols.get(i,'C')))
    for bi,ps in groups.items():
        a,b=ps[0][:2]
        if bi in wedges:order=1
        else:
            p0=ps[0][2];v=p0[-1]-p0[0];v=v/np.linalg.norm(v);normal=np.array([-v[1],v[0]]);offsets=sorted(float(np.dot(p.mean(axis=0),normal)) for _,_,p,_ in ps);clusters=[]
            for off in offsets:
                if not clusters or abs(off-clusters[-1])>1:clusters.append(off)
            order=len(clusters)
            if len(ps)==2:
                v1=ps[0][2][-1]-ps[0][2][0];v2=ps[1][2][-1]-ps[1][2][0];cross=abs(float(v1[0]*v2[1]-v1[1]*v2[0]))/(np.linalg.norm(v1)*np.linalg.norm(v2))
                if cross>.1:order=2;crossed.append((a,b))
        if order not in BT:raise ValueError('unsupported SVG multiplicity')
        rw.AddBond(a,b,BT[order])
    m=rw.GetMol();Chem.SanitizeMol(m);cf=Chem.Conformer(num);cf.Set3D(False)
    for i in range(num):cf.SetAtomPosition(i,(float(coords[i][0]),-float(coords[i][1]),0))
    m.AddConformer(cf)
    for bi,(a,b,sign) in wedges.items():m.GetBondBetweenAtoms(a,b).SetBondDir(Chem.BondDir.BEGINWEDGE if sign==1 else Chem.BondDir.BEGINDASH)
    Chem.AssignChiralTypesFromBondDirs(m,replaceExistingTags=True);Chem.SetBondStereoFromDirections(m);Chem.DetectBondStereochemistry(m);Chem.AssignStereochemistry(m,cleanIt=True,force=True)
    for a,b in crossed:m.GetBondBetweenAtoms(a,b).SetStereo(Chem.BondStereo.STEREOANY)
    label(m);return m
def numbering_check(name):
    m,_,chain=parse_name(name)
    if any(a.IsInRing() for a in m.GetAtoms()):return []
    base=re.sub(r'\([^)]*\)-','',name).replace('acido ','').split(' di ')[0];typ='none'
    for ending in ['oico','oato','ammide','nitrile','ale','one','olo','ammina']:
        if base.endswith(ending):typ=ending;break
    main=set(chain);required=set()
    if typ in ['oico','oato','ammide','nitrile','ale']:required={0}
    elif typ in ['one','olo','ammina']:
        sym='N' if typ=='ammina' else 'O';required={a for a in chain if any(x.GetSymbol()==sym and (typ!='one' or m.GetBondBetweenAtoms(a,x.GetIdx()).GetBondTypeAsDouble()==2) for x in m.GetAtomWithIdx(a).GetNeighbors())}
    rb={(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in m.GetBonds() if b.GetBondTypeAsDouble()>1 and b.GetBeginAtomIdx() in main and b.GetEndAtomIdx() in main};cs=[a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='C']
    def paths(path):
        yield path
        for nb in m.GetAtomWithIdx(path[-1]).GetNeighbors():
            j=nb.GetIdx()
            if j in cs and j not in path:yield from paths(path+[j])
    candidates=[]
    for i in cs:
        for p in paths([i]):
            if not required.issubset(p) or not all(a in p and b in p and abs(p.index(a)-p.index(b))==1 for a,b in rb):continue
            if typ in ['oico','oato','ammide','nitrile','ale'] and p[0]!=0:continue
            candidates.append(p)
    def score(p):
        suffix=tuple(sorted(p.index(a)+1 for a in required));multi=tuple(sorted(min(p.index(a),p.index(b))+1 for a,b in rb));subs=[];alphabetical=[]
        for k,i in enumerate(p,1):
            for nb in m.GetAtomWithIdx(i).GetNeighbors():
                j=nb.GetIdx()
                if j in p:continue
                if i in required and nb.GetSymbol() in ['N','O'] and typ!='none':continue
                subs.append(k);sym=nb.GetSymbol()
                if sym=='C':
                    stack=[j];branch=set()
                    while stack:
                        node=stack.pop()
                        if node in branch or node in p:continue
                        branch.add(node);stack += [x.GetIdx() for x in m.GetAtomWithIdx(node).GetNeighbors() if x.GetSymbol()=='C']
                    grp={1:'metil',2:'etil'}.get(len(branch),'ramo')
                else:grp={'Br':'bromo','Cl':'cloro','F':'fluoro','O':'idrossi','N':'ammino'}.get(sym,sym)
                alphabetical.append((grp,k))
        return (-len(p),suffix,multi,tuple(sorted(subs)),tuple(k for grp,k in sorted(alphabetical)))
    best=min(score(p) for p in candidates);return [] if score(chain)==best else [score(chain),best]
