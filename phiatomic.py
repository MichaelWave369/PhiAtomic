"""Exploratory nuclear topology; not a physical dynamics simulator."""
from __future__ import annotations
import argparse, csv, json, math, random
from dataclasses import dataclass
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
MAGIC = (2, 8, 20, 28, 50, 82, 126)
@dataclass(frozen=True)
class Nuclide:
    z: int
    n: int
    binding: float
    name: str = ''
    def __post_init__(self):
        if self.z < 1 or self.n < 0 or not math.isfinite(self.binding):
            raise ValueError('Invalid nuclide')
    @property
    def a(self): return self.z + self.n
    @property
    def ratio(self): return self.z / self.a

def orbit(c: complex, steps: int = 64):
    if steps < 1 or not math.isfinite(c.real) or not math.isfinite(c.imag):
        raise ValueError('Invalid orbit parameters')
    w = 0j
    for i in range(1, steps + 1):
        w = w*w+c
        if abs(w) > 2: return i, abs(w), True
    return steps, abs(w), False

def parameters(z, n):
    a = z+n
    return complex(z/a,0), complex(z/a,(n-z)/a), complex(2*z/a-1,(n-z)/a)

def load(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        rows = [Nuclide(int(r['z']),int(r['n']),float(r['binding_mev_per_nucleon']),r.get('name','')) for r in csv.DictReader(f)]
    if len({(r.z,r.n) for r in rows}) != len(rows): raise ValueError('Duplicate nuclides')
    return rows

def liquid_drop_features(r):
    a,z,n = r.a,r.z,r.n
    pairing = (1 if z%2==0 and n%2==0 else -1 if z%2 and n%2 else 0)/a**1.5
    return [1., -a**(-1/3), -z*(z-1)/a**(4/3), -(n-z)**2/a**2, pairing]

def shell_features(r):
    return [math.exp(-min(abs(k-m) for m in MAGIC)**2/18) for k in (r.z,r.n)]

def recursion_features(r):
    out=[]
    for c in parameters(r.z,r.n):
        i,m,e = orbit(c)
        out.extend([i/64, math.log1p(m), float(e)])
    return out

def fit_predict(train_x, train_y, test_x, alpha=0.1):
    x=np.asarray(train_x,dtype=float); t=np.asarray(test_x,dtype=float); y=np.asarray(train_y,dtype=float)
    mean=x.mean(axis=0); sd=x.std(axis=0); sd[sd<1e-12]=1
    x=(x-mean)/sd; t=(t-mean)/sd
    x=np.column_stack([np.ones(len(x)),x]); t=np.column_stack([np.ones(len(t)),t]); reg=np.eye(x.shape[1])*alpha; reg[0,0]=0
    return t@np.linalg.solve(x.T@x+reg,x.T@y)

def benchmark(rows, folds=3, seed=369):
    # Split by element, not isotope; no held-out labels used in fitting/scaling.
    groups=sorted({r.z for r in rows}); random.Random(seed).shuffle(groups)
    if folds<2 or len(groups)<folds: raise ValueError('Insufficient distinct elements')
    models={'liquid_drop':liquid_drop_features,'liquid_drop_shell':lambda r:liquid_drop_features(r)+shell_features(r),'liquid_drop_shell_recursive':lambda r:liquid_drop_features(r)+shell_features(r)+recursion_features(r)}
    errors={k:[] for k in models}; detail=[]
    for fold in range(folds):
        held=set(groups[fold::folds]); train=[r for r in rows if r.z not in held]; test=[r for r in rows if r.z in held]
        if not train or not test: raise ValueError('Empty fold')
        fold_result={'held_out_z':sorted(held),'train_n':len(train),'test_n':len(test),'mae':{}}
        for name,fn in models.items():
            pred=fit_predict([fn(r) for r in train],[r.binding for r in train],[fn(r) for r in test]); e=np.abs(pred-np.array([r.binding for r in test])); errors[name].extend(e.tolist()); fold_result['mae'][name]=float(e.mean())
        detail.append(fold_result)
    return {'seed':seed,'n':len(rows),'folds':detail,'pooled_mae_mev_per_nucleon':{k:float(np.mean(v)) for k,v in errors.items()},'warning':'Small illustrative datasets cannot establish predictive validity; mappings are heuristic.'}

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--data',type=Path,default=ROOT/'data/demo.csv'); p.add_argument('--folds',type=int,default=3); p.add_argument('--seed',type=int,default=369); p.add_argument('--output',type=Path); args=p.parse_args(); rows=load(args.data)
    result=benchmark(rows,args.folds,args.seed); result['literal_escape']=[{'name':r.name,'ratio':r.ratio,'iteration':orbit(complex(r.ratio))[0],'escaped':orbit(complex(r.ratio))[2]} for r in rows]; result['input']=str(args.data); s=json.dumps(result,indent=2,allow_nan=False); print(s)
    if args.output: args.output.write_text(s+'\n',encoding='utf-8')
if __name__=='__main__': main()
