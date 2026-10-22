import json,re,subprocess
from pathlib import Path

BASE='640c01c42d4666cc01ecaa2d660796c0af80639b'
ROOT=Path(__file__).resolve().parents[1]
def load(baseline=False):
    def read(f):return subprocess.check_output(['git','show',BASE+':'+f],cwd=ROOT,text=True) if baseline else (ROOT/f).read_text()
    html=read('index.html');ext=json.loads(read('organic-v25-extension.js').split('=',1)[1].strip().rstrip(';'))
    raw=json.loads(re.search(r'const DATA = (.*);\n',html)[1]);g=json.loads(re.search(r'const GUIDE = (.*);\n',html)[1])
    for i,fix in ext['legacyFixes'].items():raw['molecules'][int(i)-1].update(fix)
    raw['molecules']+=ext['molecules'];raw['isomerism']+=ext['extraRS']
    for q in raw['isomerism']:
        if q['kind']=='RS' and q['id'] in ext['existingRS']:q['cipDetails']=ext['existingRS'][q['id']]
    g.update(ext['guide']);g.update(ext['legacyGuide']);raw['guide']=g
    return raw
