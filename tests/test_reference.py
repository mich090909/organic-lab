import unittest,xml.etree.ElementTree as ET
from catalog import load
from chemistry import *
class ReferenceChemistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.old=load(True);cls.new=load();learn_glyphs(m['svg'] for m in cls.new['molecules'])
    def test_reported_fischer_is_2S_3S(self):
        _,c=read_fischer(self.old['molecules'][732]['fischerSvg']);self.assertEqual([c[k]['answer'] for k in ['2','3']],['S','S']);self.assertEqual(c['3']['sides'],['left','top','bottom','right'])
    def test_carbon_branches_and_swapped_substituents(self):
        _,c=read_fischer(self.old['molecules'][732]['fischerSvg']);self.assertEqual(c['3']['sides'][1:3],['top','bottom']);root=ET.fromstring(self.old['molecules'][732]['fischerSvg']);ts=[x for x in root.iter() if x.tag.endswith('text')];a=next(x for x in ts if x.get('x')=='67' and x.text=='Br');y=a.get('y');b=next(x for x in ts if x.get('x')=='293' and x.get('y')==y);a.text,b.text=b.text,a.text;_,c=read_fischer(ET.tostring(root,encoding='unicode'));self.assertEqual([c[k]['answer'] for k in ['2','3']],['S','R'])
    def test_fischer_rotation_parity(self):
        sides=['left','top','bottom','right'];self.assertEqual(geometric_rs(sides),'S');r180={'left':'right','right':'left','top':'bottom','bottom':'top'};r90={'left':'top','top':'right','right':'bottom','bottom':'left'};self.assertEqual(geometric_rs([r180[s] for s in sides]),'S');self.assertEqual(geometric_rs([r90[s] for s in sides]),'R')
    def test_modern_cip_multiple_bonds(self):
        m=Chem.MolFromSmiles('N[C@H](C(=O)O)CO');mh,order=priorities(m,1);self.assertEqual(order[:3],[0,2,5]);self.assertEqual(mh.GetAtomWithIdx(order[3]).GetSymbol(),'H')
    def test_modern_cip_ring(self):
        m=Chem.MolFromSmiles('F[C@H](C1CC1)CC');_,order=priorities(m,1);self.assertEqual(order[:3],[0,2,5])
    def test_modern_cip_first_difference_beyond_first_shell(self):
        m=Chem.MolFromSmiles('F[C@H](CCCl)CCC');_,order=priorities(m,1);self.assertEqual(order[:3],[0,2,5])
    def test_visible_skeletal_mirror_not_saved_descriptor(self):
        a=svg_graph(self.new['molecules'][46]['svg']);b=svg_graph(self.new['molecules'][47]['svg']);self.assertNotEqual(canonical(a),canonical(b));self.assertEqual(list(label(a).values()),['R']);self.assertEqual(list(label(b).values()),['S'])
    def test_all_graphs_and_ids_preserved(self):
        self.assertEqual([m['id'] for m in self.old['molecules']],[m['id'] for m in self.new['molecules']]);self.assertEqual([canonical(Chem.MolFromSmiles(m['smiles'])) for m in self.old['molecules']],[canonical(Chem.MolFromSmiles(m['smiles'])) for m in self.new['molecules']])
    def test_unknown_names_fail_closed(self):
        with self.assertRaises(ValueError):parse_name('nome non supportato')
    def test_finite_ring_nomenclature_reference(self):
        expected={31:'ciclopropano',32:'ciclobutano',33:'cicloesano',34:'metilcicloesano',35:'1,2-dimetilcicloesano',36:'ciclopentene',37:'benzene',38:'metilbenzene',39:'etilbenzene',40:'1,4-dimetilbenzene',53:'fenolo',65:'cicloesanone',73:'acido benzoico',87:'benzenammina',121:'clorobenzene',122:'1,4-diclorobenzene',137:'benzoato di metile',138:'benzoato di etile',145:'ciclopentanone',146:'1-feniletan-1-one',154:'benzonitrile'};actual={}
        for d in self.new['molecules']:
            m=Chem.MolFromSmiles(d['smiles'])
            if any(a.IsInRing() for a in m.GetAtoms()):actual[d['id']]=d['name']
        self.assertEqual(actual,expected)
    def test_explanations_distinguish_N_substitution_and_terminal_locants(self):
        for d in self.new['molecules'][530:562]:self.assertIn('anche quando è terminale',d['why'])
        for d in self.new['molecules'][612:632]:self.assertIn('N-sostituite',d['why'])
        for d in self.new['molecules'][641:664]:self.assertIn('C(=O)-N',d['why'])
if __name__=='__main__':unittest.main()
