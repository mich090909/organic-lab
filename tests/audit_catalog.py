"""Exhaustive independent graph/name/geometry audit of all effective entries."""
import csv,json,re,collections
from rdkit import Chem,rdBase
from rdkit.Chem import rdMolDescriptors
from catalog import load,ROOT
from chemistry import *
COMMON={'isottano':'CC(C)CC(C)(C)C','toluene':'Cc1ccccc1','p-xilene':'Cc1ccc(C)cc1','para-xilene':'Cc1ccc(C)cc1','dimetiletere':'COC','etilmetiletere':'CCOC','formaldeide':'C=O','acetaldeide':'CC=O','acetone':'CC(=O)C','acido formico':'O=CO','acido acetico':'CC(=O)O','acido lattico':'CC(O)C(=O)O','acido glicolico':'OCC(=O)O','formiato di metile':'COC=O','acetato di metile':'COC(C)=O','acetato di etile':'CCOC(C)=O','metilammina':'CN','etilammina':'CCN','propilammina':'CCCN','isopropilammina':'CC(C)N','etilmetilammina':'CCNC','trimetilammina':'CN(C)C','anilina':'Nc1ccccc1','acetammide':'CC(N)=O','acetonitrile':'CC#N','glicina':'NCC(=O)O','alanina':'CC(N)C(=O)O','acetoacetato di metile':'COC(=O)CC(C)=O','lattato di etile':'CCOC(=O)C(C)O','p-diclorobenzene':'Clc1ccc(Cl)cc1','acetofenone':'CC(=O)c1ccccc1','propionitrile':'CCC#N','beta-alanina':'NCCC(=O)O','acido gamma-amminobutirrico':'NCCCC(=O)O','piruvato di metile':'COC(=O)C(C)=O'}
def alias_check(alias,ref):
    if alias in COMMON:assert canonical(Chem.MolFromSmiles(COMMON[alias]),False)==canonical(ref,False),alias;return
    alias=alias.replace('amino','ammino').replace('oxo','osso').replace('trans-','(E)-').replace('cis-','(Z)-');roots='|'.join(sorted(ROOTS,key=len,reverse=True))
    alias=re.sub(r'(\d+(?:,\d+)*)-('+roots+r')an(olo|one|ammina)',r'\2an-\1-\3',alias);alias=re.sub(r'(\d+(?:,\d+)*)-('+roots+r')(ene|ino)',r'\2-\1-\3',alias);alias=re.sub(r'(\d+(?:,\d+)*)-('+roots+r')andiolo',r'\2ano-\1-diolo',alias);alias=re.sub(r'(\d+)-('+roots+r')anolo',r'\2an-\1-olo',alias);alias=re.sub(r'(\d+)-('+roots+r')anone',r'\2an-\1-one',alias);alias=alias.replace('cloroetanolo','cloroetan-1-olo');name_stereo_check({'name':alias,'smiles':canonical(ref)})
def mirror(m):
    x=Chem.Mol(m)
    for a in x.GetAtoms():a.InvertChirality()
    return x
def carbon_skeleton(m):
    rw=Chem.RWMol();mp={}
    for a in m.GetAtoms():
        if a.GetSymbol()=='C':mp[a.GetIdx()]=rw.AddAtom(Chem.Atom('C'))
    for b in m.GetBonds():
        a,c=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if a in mp and c in mp:rw.AddBond(mp[a],mp[c],BT[1])
    mol=rw.GetMol();Chem.SanitizeMol(mol);return canonical(mol,False)
def functional_signature(m):
    patterns=['[CX3H1](=O)','[#6][CX3](=O)[#6]','[CX3](=O)[OX2H1]','[CX3](=O)[OX2][#6]','[OX2H1][CX4]','[#6][OX2][#6]','[NX3;!$(N-C=O)]','C#N','C=C','C#C','[F,Cl,Br,I]']
    return tuple(len(m.GetSubstructMatches(Chem.MolFromSmarts(s))) for s in patterns)
def run():
    D=load();M=D['molecules'];Q=D['isomerism'];learn_glyphs(m['svg'] for m in M);rows=[];fail=[];fcount=0;aliases=0;seen=set();refs={}
    for d in M:
        status='PASS';notes=''
        try:
            ref,nm,match,sts=name_stereo_check(d);refs[d['id']]=ref;assert not numbering_check(d['name']),'non-minimal numbering';assert rdMolDescriptors.CalcMolFormula(ref)==d['formula'],'formula';assert canonical(ref) not in seen,'duplicate structure';seen.add(canonical(ref));sm=svg_graph(d['svg']);assert canonical(sm)==canonical(ref),'visible skeletal graph/stereo'
            pot=Chem.FindPotentialStereo(ref);actual=sum(x.type==Chem.StereoType.Atom_Tetrahedral for x in pot);assert d['stereoCenters']==actual,'stereocenter count';assert d['unspecifiedStereo']==any(x.specified==Chem.StereoSpecified.Unspecified for x in pot),'unspecified stereo flag'
            if re.search(r'gruppo -NH[₂2]|principale è -NH[₂2]',d['why']) and 'non hanno necessariamente' not in d['why']:assert any(a.GetSymbol()=='N' and a.GetTotalNumHs()==2 for a in ref.GetAtoms()),'false NH2 explanation'
            assert d['why'].strip();g=D['guide'][str(d['id'])];assert g['length']==len(parse_name(d['name'])[2]);assert g['prefix']==prefix_for_name(d['name']);nums=re.findall(r'\d+',re.sub(r'\([^)]*\)-','',d['name']));assert g['locants']==(','.join(map(str,sorted(map(int,nums)))) or 'nessuno');assert [(str(x['position']),x['value']) for x in g['stereo']]==sts
            for alias in d.get('aliases',[]):alias_check(alias,ref);aliases+=1
            if d.get('fischerSvg'):
                fm,rs=read_fischer(d['fischerSvg']);assert canonical(fm)==canonical(ref),'Fischer graph/stereo';fcount+=1
                for key,calc in rs.items():
                    stored=d['cipFischer'][key];assert stored['answer']==calc['answer']==stored['geometryAnswer'];assert [v['side'] for v in stored['priorities']]==calc['sides'],'Fischer CIP order'
        except Exception as e:status='FAIL';notes=str(e);fail.append(('molecule',d['id'],notes))
        rows.append({'type':'molecule','id':d['id'],'category':d['family'],'name':d['name'],'status':status,'notes':notes})
    for q in Q:
        status='PASS';notes=''
        try:
            kind=q['kind']
            if kind in ['pair','stereo_pair']:
                a,b=refs[q['a']],refs[q['b']]
                if kind=='stereo_pair':
                    assert canonical(a,False)==canonical(b,False);expected='identici' if canonical(a)==canonical(b) else 'enantiomeri' if canonical(mirror(a))==canonical(b) else 'diastereoisomeri'
                elif rdMolDescriptors.CalcMolFormula(a)!=rdMolDescriptors.CalcMolFormula(b):expected='non_isomeri'
                elif functional_signature(a)!=functional_signature(b):expected='funzione'
                elif carbon_skeleton(a)!=carbon_skeleton(b):expected='catena'
                else:expected='posizione'
                assert q['answer']==expected,(q['answer'],expected)
            else:
                d=M[q['moleculeId']-1] if q.get('moleculeId') else q;actual=svg_graph(d['svg']);assert canonical(actual)==canonical(Chem.MolFromSmiles(d['smiles']));codes=label(actual)
                if kind=='meso':expected='meso' if canonical(mirror(actual))==canonical(actual) else 'chirale'
                elif kind=='EZ':
                    bonds=[b.GetProp('_CIPCode') for b in actual.GetBonds() if b.HasProp('_CIPCode')];assert len(bonds)==1;expected=bonds[0]
                elif q.get('centers'):
                    _,nm,_,sts=name_stereo_check(d);mt=actual.GetSubstructMatch(nm);assert mt;expected=','.join(str(c['position'])+codes[mt[int(c['position'])-1]] for c in q['centers'])
                    for c in q['centers']:
                        key=str(c['position']);idx=mt[int(c['position'])-1];saved=q['cipDetails'][key];assert saved['answer']==codes[idx];ref=Chem.MolFromSmiles(d['smiles']);_,_,match,_=name_stereo_check(d);_,ro=priorities(ref,match[int(c['position'])-1]);assert [v['atomIndex'] for v in saved['priorities']]==ro
                else:assert len(codes)==1;expected=next(iter(codes.values()));assert q['cipDetails']['answer']==expected
                assert q['answer']==expected,(q['answer'],expected)
            assert q['auditStatus']=='checked'
        except Exception as e:status='FAIL';notes=str(e);fail.append(('question',q['id'],notes))
        rows.append({'type':'question','id':q['id'],'category':q['kind'],'name':q.get('name',''),'status':status,'notes':notes})
    assert len(M)==912 and len(Q)==142
    summary={'engine':'RDKit '+rdBase.rdkitVersion+' modern CIPLabeler','molecules':len(M),'questions':len(Q),'families':dict(collections.Counter(m['family'] for m in M)),'questionCategories':dict(collections.Counter(q['kind'] for q in Q)),'fischer':fcount,'aliasesChecked':aliases,'failures':fail,'rows':len(rows),'excluded':[]}
    (ROOT/'audit/summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    with (ROOT/'audit/inventory.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(summary,ensure_ascii=False,indent=2));assert not fail,fail
if __name__=='__main__':run()
