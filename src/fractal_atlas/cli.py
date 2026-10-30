"""Generic dataset-driven Julia atlas. Artistic encoding, not inference."""
from __future__ import annotations
import argparse
import base64
import io
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from PIL import Image

def read_table(path):
    suffix = path.suffix.lower()
    if suffix in (".csv", ".tsv"):
        return pd.read_csv(path, sep="\t" if suffix == ".tsv" else ",")
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError("Supported input types: CSV, TSV, Parquet")

def prepare(config, config_path):
    d = config["data"]
    path = (config_path.parent / d["path"]).resolve()
    df = read_table(path)
    kind = d.get("kind", "features")
    if kind == "features":
        id_col = d.get("id_column", "id")
        ids = df[id_col].astype(str).tolist()
        if len(set(ids)) != len(ids):
            raise ValueError("Entity identifiers must be unique")
        columns = d.get("columns") or [x for x in df.select_dtypes(include="number").columns if x != id_col]
        if not columns:
            raise ValueError("No numeric features selected")
        values = df[columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        policy = d.get("missing", "median")
        if policy == "drop":
            keep = np.isfinite(values).all(axis=1)
            values = values[keep]; ids = [x for x, k in zip(ids, keep) if k]
        elif policy in ("median", "zero"):
            fill = np.nanmedian(np.where(np.isfinite(values), values, np.nan), axis=0) if policy == "median" else np.zeros(values.shape[1])
            fill = np.nan_to_num(fill)
            values = np.where(np.isfinite(values), values, fill)
        elif not np.isfinite(values).all():
            raise ValueError("Missing values found: use missing: median, zero, or drop")
        if len(ids) == 0:
            raise ValueError("No entities remain")
        normal = d.get("normalization", "zscore")
        if normal == "zscore":
            std = values.std(axis=0); values = (values - values.mean(axis=0)) / np.where(std > 0, std, 1)
        elif normal == "minmax":
            rng = np.ptp(values, axis=0); values = (values - values.min(axis=0)) / np.where(rng > 0, rng, 1)
        elif normal != "none":
            raise ValueError("Unknown normalization: " + normal)
        if len(ids)**2 > int(d.get("max_cells", 4000000)):
            raise ValueError("Input exceeds configured max_cells; no automatic downsampling")
        method = d.get("similarity", "cosine")
        if method == "cosine":
            norm = np.linalg.norm(values, axis=1)
            v = values / np.where(norm[:, None] > 0, norm[:, None], 1)
            scores = v @ v.T
        elif method == "pearson":
            z = values - values.mean(axis=1, keepdims=True)
            norm = np.linalg.norm(z, axis=1)
            v = z / np.where(norm[:, None] > 0, norm[:, None], 1)
            scores = v @ v.T
        elif method == "rbf":
            norm2 = np.sum(values**2, axis=1)
            distances = np.maximum(norm2[:,None] + norm2[None,:] - 2 * values @ values.T, 0)
            positive = distances[distances > 0]
            scale = np.median(positive) if positive.size else 1.0
            scores = np.exp(-distances / max(scale, 1e-10))
        else:
            raise ValueError("Unknown similarity: " + method)
        features = values.astype(np.float32)
    elif kind == "relationships":
        a = d.get("pair_source_column", "source")
        b = d.get("pair_target_column", "target")
        v = d.get("pair_value_column", "value")
        ids = list(dict.fromkeys(df[a].astype(str).tolist() + df[b].astype(str).tolist()))
        if len(ids)**2 > int(d.get("max_cells", 4000000)):
            raise ValueError("Input exceeds configured max_cells")
        lookup = {k:i for i,k in enumerate(ids)}
        scores = np.full((len(ids),len(ids)), np.nan, dtype=np.float32)
        symmetric = bool(d.get("symmetric", True))
        for record in df[[a,b,v]].itertuples(index=False, name=None):
            i,j = lookup[str(record[0])],lookup[str(record[1])]
            val = float(record[2])
            if np.isfinite(scores[i,j]):
                raise ValueError("Duplicate relationship")
            scores[i,j] = val
            if symmetric:
                scores[j,i] = val
        features = None
        columns = []
    else:
        raise ValueError("Unknown data kind")
    return ids, np.asarray(scores, dtype=np.float32), features, columns, kind

def region_stats(scores, rect):
    r0,r1,c0,c1=rect
    part=scores[r0:r1,c0:c1]
    valid=part[np.isfinite(part)]
    if not valid.size:
        return {"mean":0.0,"std":0.0,"coverage":0.0,"c_real":-0.745,"c_imag":0.186}
    m=float(valid.mean()); s=float(valid.std())
    # Stable bounded mapping to parameters. Fixed across datasets.
    real = -0.745 + 0.17*np.tanh(m)
    imag = 0.186 + 0.17*np.tanh(s)
    return {"mean":m,"std":s,"coverage":float(valid.size/part.size),"c_real":float(real),"c_imag":float(imag)}

def julia(c, size=450, max_iter=100, center=(0.0,0.0), scale=1.5):
    yy,xx=np.mgrid[0:size,0:size]
    z=(xx/size*2-1)*scale + center[0]+1j*((yy/size*2-1)*scale+center[1])
    iterations=np.full(z.shape,max_iter,dtype=np.uint16)
    alive=np.ones(z.shape,dtype=bool)
    for n in range(max_iter):
        z[alive]=z[alive]**2 + c
        escaped=alive & (np.abs(z)>4)
        iterations[escaped]=n
        alive[escaped]=False
        if not alive.any():
            break
    return iterations

def build(config_path, output):
    config_path=Path(config_path).resolve()
    cfg=yaml.safe_load(config_path.read_text()) if config_path.suffix.lower() in (".yaml",".yml") else json.loads(config_path.read_text())
    ids,scores,features,columns,kind=prepare(cfg,config_path)
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    size=int(cfg.get("render",{}).get("max_iter",100))
    size=max(10,min(size,1000))
    root=region_stats(scores,(0,len(ids),0,len(ids)))
    image=julia(complex(root["c_real"],root["c_imag"]),size=500,max_iter=size)
    pic=(image / size * 255).astype("uint8")
    Image.fromarray(pic,mode="L").save(out/"preview.png")
    payload={"name":cfg.get("name","Fractal Atlas"),"ids":ids,"scores":scores.tolist(),"features":features.tolist() if features is not None else None,"columns":columns,"kind":kind,"max_iter":size}
    # Embed JSON as a non-executable script tag, replacing HTML-active chars.
    safe_json=json.dumps(payload,allow_nan=False).replace("<","\\u003c").replace(">","\\u003e").replace("&","\\u0026")
    # Convert missing scores to JSON null.
    safe_json=safe_json.replace("NaN","null")
    page=HTML.replace("__ATLAS_DATA__",safe_json)
    (out/"index.html").write_text(page,encoding="utf8")
    summary={"name":payload["name"],"entities":len(ids),"features":columns,"kind":kind,"parameters":root}
    (out/"summary.json").write_text(json.dumps(summary,indent=2))
    print("Created",out/"index.html")

HTML = r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Fractal Atlas</title><style>
body{background:#10141b;color:#e7eaf0;font:15px system-ui;margin:0;padding:22px}
main{max-width:1080px;margin:auto}button,input{font:inherit;border-radius:6px;padding:8px;margin:4px}
button{background:#233345;color:white;border:1px solid #5c7087;cursor:pointer}
button:hover{background:#345571}canvas{background:#000;width:100%;height:auto;max-width:540px;cursor:crosshair;border-radius:8px}
.layout{display:grid;grid-template-columns:minmax(300px,1fr) minmax(270px,.8fr);gap:24px}
.panel{background:#1b2530;padding:14px;border-radius:10px}small{color:#b3bdcb}
@media(max-width:700px){.layout{grid-template-columns:1fr}}
</style><main><h1 id="heading">Fractal Atlas</h1>
<p>Data zoom selects a quadrant of the relationship matrix. Mathematical zoom magnifies the Julia set without changing data.</p>
<div><button id="back">Back a data level</button><button id="reset">Reset</button><button id="mathin">Math zoom +</button><button id="mathout">Math zoom −</button></div>
<div class="layout"><section><canvas id="canvas" width="480" height="480" aria-label="Julia set for selected data region"></canvas><small>Click a quadrant in the grid to descend; fractal canvas shows the selected region.</small></section>
<section class="panel"><h2>Relationship grid</h2><canvas id="grid" width="320" height="320" aria-label="Data region quadrant selection"></canvas><div id="detail"></div></section></div>
</main><script id="data" type="application/json">__ATLAS_DATA__</script><script>
const a=JSON.parse(document.getElementById('data').textContent);
document.getElementById('heading').textContent=a.name;
const N=a.ids.length;
let stack=[[0,N,0,N]], mag=1;
const canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d');
const grid=document.getElementById('grid'),g=grid.getContext('2d');
function stats(b){let [r0,r1,c0,c1]=b,x=[],total=(r1-r0)*(c1-c0);
for(let i=r0;i<r1;i++)for(let j=c0;j<c1;j++)if(a.scores[i][j]!==null)x.push(a.scores[i][j]);
let m=x.reduce((v,n)=>v+n,0)/Math.max(x.length,1);
let v=x.reduce((v,n)=>v+(n-m)**2,0)/Math.max(x.length,1);
return {mean:m,std:Math.sqrt(v),coverage:x.length/total,c:[-.745+.17*Math.tanh(m),.186+.17*Math.tanh(Math.sqrt(v))]};
}
function children(b){let [r0,r1,c0,c1]=b,rm=Math.floor((r0+r1)/2),cm=Math.floor((c0+c1)/2);
return [[r0,rm,c0,cm],[r0,rm,cm,c1],[rm,r1,c0,cm],[rm,r1,cm,c1]].filter(b=>b[0]<b[1]&&b[2]<b[3]);}
function drawJulia(){
let st=stats(stack.at(-1)),w=canvas.width,im=ctx.createImageData(w,w), max=a.max_iter;
for(let y=0;y<w;y++)for(let x=0;x<w;x++){
let zr=(x/w*2-1)*1.5/mag,zi=(y/w*2-1)*1.5/mag,n=0;
for(;n<max;n++){let r=zr*zr-zi*zi+st.c[0],i=2*zr*zi+st.c[1];zr=r;zi=i;if(zr*zr+zi*zi>16)break;}
let p=(y*w+x)*4,t=n/max;im.data[p]=n===max?5:Math.floor(24+210*t);
im.data[p+1]=n===max?7:Math.floor(30+140*Math.sqrt(t));
im.data[p+2]=n===max?13:Math.floor(60+190*(1-t));
im.data[p+3]=255;}
ctx.putImageData(im,0,0);
}
function draw(){let b=stack.at(-1),cs=children(b);
g.fillStyle='#182332';g.fillRect(0,0,320,320);
for(const c of cs){let s=stats(c);
let x=(c[2]-b[2])/(b[3]-b[2])*320,y=(c[0]-b[0])/(b[1]-b[0])*320,
ww=(c[3]-c[2])/(b[3]-b[2])*320,hh=(c[1]-c[0])/(b[1]-b[0])*320;
g.fillStyle='hsl('+Math.round(210-145*(s.mean+1)/2)+' 62% 36%)';
g.fillRect(x,y,ww,hh);g.strokeStyle='#ddd';g.strokeRect(x,y,ww,hh);
}
let s=stats(b),detail='Level '+(stack.length-1)+' · Rows '+b[0]+'–'+(b[1]-1)+' · Columns '+b[2]+'–'+(b[3]-1)+'<p>Mean: '+s.mean.toFixed(5)+' · Coverage: '+(100*s.coverage).toFixed(1)+'%</p>';
if((b[1]-b[0])===1&&(b[3]-b[2])===1){detail+='<p>Pair: '+a.ids[b[0]]+' × '+a.ids[b[2]]+'</p><p>Score: '+a.scores[b[0]][b[2]]+'</p>';if(a.features){detail+='<p>Features available: '+a.columns.length+'</p>';}}
document.getElementById('detail').textContent='';
document.getElementById('detail').innerText=detail.replace(/<[^>]*>/g,' ');
drawJulia();}
grid.addEventListener('click',ev=>{let b=stack.at(-1),cs=children(b),rect=grid.getBoundingClientRect(),xx=(ev.clientX-rect.left)/rect.width,yy=(ev.clientY-rect.top)/rect.height;
let row=b[0]+Math.min(b[1]-b[0]-1,Math.floor(yy*(b[1]-b[0]))),col=b[2]+Math.min(b[3]-b[2]-1,Math.floor(xx*(b[3]-b[2])));
let c=cs.find(v=>row>=v[0]&&row<v[1]&&col>=v[2]&&col<v[3]);
if(c){stack.push(c);mag=1;draw();}});
document.getElementById('back').onclick=()=>{if(stack.length>1)stack.pop();mag=1;draw();};
document.getElementById('reset').onclick=()=>{stack=[[0,N,0,N]];mag=1;draw();};
document.getElementById('mathin').onclick=()=>{mag=Math.min(mag*2,1024);drawJulia();};
document.getElementById('mathout').onclick=()=>{mag=Math.max(mag/2,1);drawJulia();};draw();
</script></html>'''
def main(argv=None):
    parser=argparse.ArgumentParser(description="Build a dataset-independent fractal atlas")
    parser.add_argument("--config","-c",required=True)
    parser.add_argument("--out","-o",required=True)
    args=parser.parse_args(argv)
    build(args.config,args.out)
if __name__=="__main__":
    main()
