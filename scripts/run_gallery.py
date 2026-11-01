"""Reproducible biological Fractal Atlas gallery. Publishes derived images only."""
from __future__ import annotations
import csv
import json
import shutil
import tempfile
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps
from sklearn.datasets import load_iris, load_breast_cancer
from fractal_atlas.cli import build

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "results" / "gallery"
WEB = ROOT / "docs" / "gallery"
ITER = 110

SOURCES = {
    "iris": "https://archive.ics.uci.edu/dataset/53/iris",
    "human": "https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic",
    "maize": "https://github.com/corn2code/TWAStutorial/blob/main/merged_gene_tpms.csv",
    "mouse": "https://archive.ics.uci.edu/dataset/342/mice+protein+expression",
    "ecoli": "https://github.com/SBRG/precise-db",
    "virus": "https://github.com/nextstrain/ncov/blob/master/data/references_sequences.fasta",
}
RAW = {
    "maize": "https://raw.githubusercontent.com/corn2code/TWAStutorial/main/merged_gene_tpms.csv",
    "mouse": "https://raw.githubusercontent.com/kathleenmagbual/MiceProtein_MLClassification/master/Data_Cortex_Nuclear.csv",
    "ecoli": "https://raw.githubusercontent.com/SBRG/precise-db/master/data/A.csv",
    "virus": "https://raw.githubusercontent.com/nextstrain/ncov/master/data/references_sequences.fasta",
}

def download(which, directory):
    target=directory/(which + (".fasta" if which=="virus" else ".csv"))
    if target.exists(): return target
    for attempt in range(3):
        try:
            request=urllib.request.Request(RAW[which],headers={"User-Agent":"FractalAtlas/0.1"})
            with urllib.request.urlopen(request,timeout=75) as inp, target.open("wb") as out:
                shutil.copyfileobj(inp,out)
            if target.stat().st_size < 100: raise ValueError("Unexpectedly small download")
            return target
        except Exception:
            if target.exists(): target.unlink()
            if attempt == 2: raise
            time.sleep(attempt+1)

def dataset(items,group,slug,label,frame,source,description):
    if "id" not in frame.columns or frame["id"].isna().any() or frame["id"].duplicated().any():
        raise ValueError("Invalid entity identifiers: "+slug)
    if len(frame)<20 or len(frame.columns)<3: raise ValueError("Insufficient data: "+slug)
    items.append(dict(group=group,slug=slug,label=label,frame=frame.reset_index(drop=True),
                      source=source,description=description))

def make_inputs(temp):
    items=[]
    data=load_iris(as_frame=True)
    frame=data.data.copy()
    frame.insert(0,"id",["iris-%04d"%i for i in range(len(frame))])
    for index,name in enumerate(data.target_names):
        dataset(items,"Plants","iris_"+name,"Iris — "+name,frame.loc[data.target==index].copy(),
                SOURCES["iris"],"Flower morphology, one iris species")

    data=load_breast_cancer(as_frame=True)
    frame=data.data.copy()
    frame.insert(0,"id",["nucleus-%04d"%i for i in range(len(frame))])
    for value,label in [(0,"malignant"),(1,"benign")]:
        dataset(items,"Human","human_"+label,"Human breast nuclei — "+label,
                frame.loc[data.target==value].head(180).copy(),SOURCES["human"],
                "30 measured nucleus features; diagnostic class")

    maize=pd.read_csv(download("maize",temp))
    values=[c for c in maize if c!="GeneID"]
    maize[values]=maize[values].apply(pd.to_numeric,errors="coerce")
    maize["_mean"]=maize[values].mean(axis=1)
    maize=maize[maize["_mean"].gt(0.02) & np.isfinite(maize["_mean"])]
    maize=maize.drop_duplicates("GeneID").sort_values(["_mean","GeneID"])
    for label,sub in [("lower",maize.head(180)),("higher",maize.tail(180))]:
        frame=sub[["GeneID"]+values].rename(columns={"GeneID":"id"}).copy()
        frame[values]=np.log1p(frame[values].clip(lower=0))
        dataset(items,"Plants","maize_"+label,"Maize RNA expression — "+label,frame,SOURCES["maize"],
                "Transcript TPM across six samples; lower/higher expression cohorts")

    mouse=pd.read_csv(download("mouse",temp))
    proteins=[c for c in mouse.columns if c.endswith("_N")]
    if "Genotype" not in mouse or len(proteins)<40:
        raise ValueError("Unexpected mouse cortex schema")
    for code,label in [("c","control"),("t","trisomy")]:
        values=mouse["Genotype"].astype(str).str.strip().str.lower()
        sub=mouse.loc[values.isin(([code,"control"] if code=="c" else [code,"ts65dn","trisomy"]))].head(180)
        frame=sub[["MouseID"]+proteins].rename(columns={"MouseID":"id"}).copy()
        dataset(items,"Animals","mouse_"+label,"Mouse cortex proteins — "+label,frame,SOURCES["mouse"],
                "77 protein measurements in cortex; mouse genotype group")

    ecoli=pd.read_csv(download("ecoli",temp),index_col=0)
    ecoli=ecoli.apply(pd.to_numeric,errors="coerce")
    for label,columns in [("all",list(ecoli.columns)),
                          ("iron",[c for c in ecoli.columns if c.startswith("fur__")])]:
        frame=ecoli[columns].reset_index().rename(columns={ecoli.index.name or "index":"id"})
        frame=frame.rename(columns={frame.columns[0]:"id"})
        frame["id"]=frame["id"].astype(str)
        dataset(items,"Bacteria","ecoli_"+label,"E. coli iModulon activity — "+label,frame,
                SOURCES["ecoli"],"iModulon activities across experimental conditions, not bacterial cells")

    seqs={}
    current=None
    for line in download("virus",temp).read_text().splitlines():
        if line.startswith(">"):
            current=line[1:].split()[0];seqs[current]=""
        elif current:seqs[current]+=line.strip().upper()
    if len(seqs)<2:raise ValueError("Expected two SARS-CoV-2 reference sequences")
    kmers=[a+b for a in "ACGT" for b in "ACGT"]
    for i,(name,sequence) in enumerate(list(seqs.items())[:2]):
        rows=[]
        for start in range(0,len(sequence)-159,160):
            segment=sequence[start:start+160]
            pairs=[segment[j:j+2] for j in range(len(segment)-1)]
            valid=[p for p in pairs if len(p)==2 and all(k in "ACGT" for k in p)]
            if len(valid)<120:continue
            row={"id":f"{name}:{start+1}-{start+160}"}
            row.update({"dinuc_"+x:valid.count(x)/len(valid) for x in kmers})
            rows.append(row)
        dataset(items,"Viruses","sars_ref_"+str(i+1),
                "SARS-CoV-2 reference — "+name,pd.DataFrame(rows),SOURCES["virus"],
                "160-nt windows with 16 dinucleotide frequencies; viruses have no cell types")
    return items

def metrics(filename):
    a=np.asarray(Image.open(filename).convert("L"),dtype=np.uint8)
    hist=np.bincount(a.ravel(),minlength=256)/a.size
    active=hist[hist>0]
    return dict(mean_gray=round(float(a.mean()),4),
                entropy_bits=round(float(-(active*np.log2(active)).sum()),4),
                black_pct=round(float((a<=5).mean()*100),3))

def contact_sheet(records):
    width,height,cols=310,360,4
    sheet=Image.new("RGB",(cols*width,((len(records)+cols-1)//cols)*height),"#141d28")
    d=ImageDraw.Draw(sheet)
    for k,row in enumerate(records):
        x=(k%cols)*width;y=(k//cols)*height
        pic=ImageOps.fit(Image.open(DEST/(row["slug"]+".png")).convert("RGB"),(290,290))
        sheet.paste(pic,(x+10,y+8))
        d.text((x+12,y+305),row["label"][:42],fill="white")
        d.text((x+12,y+328),f'{row["rows"]} entities | {row["features"]} features',fill="#adc7df")
    sheet.save(DEST/"contact_sheet.png",optimize=True)

def web_gallery(records):
    import html
    rec=[{"label":x["label"],"slug":x["slug"],"c":x["julia_c"]} for x in records]
    cards="".join(
        '<article><a href="images/'+x["slug"]+'.png"><img src="images/'+
        x["slug"]+'.png" alt="Fractal portrait" loading="lazy"></a><h3>'+
        html.escape(x["label"])+'</h3><small>'+str(x["rows"])+' entities · '+str(x["features"])+
        ' features</small><p><a href="'+x["source"]+'">Data source</a></p></article>'
        for x in records)
    template=r"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Fractal Atlas — gallery</title>
<style>
:root{color-scheme:dark}body{font:15px system-ui,Arial;background:#0d1420;color:#eef3f8;margin:0}
main{max-width:1200px;margin:auto;padding:24px}p{max-width:890px;line-height:1.55;color:#b8c7d6}
a{color:#95c8fd}.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(235px,1fr));gap:16px}
article{padding:12px;background:#172232;border:1px solid #334459;border-radius:8px}
article img{width:100%;aspect-ratio:1}article h3{font-size:1rem}small{color:#adbdcd}
.compare{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}
canvas{width:100%;max-width:420px;background:#000}select,input{font:inherit;margin:6px;padding:8px}
@media(max-width:600px){.compare{grid-template-columns:1fr}}
</style><main><h1>Fractal Atlas · biological dataset gallery</h1>
<p>Different real datasets, the same mathematical rule. These portraits explore how the data changes fractal appearance; they are not claims of biological similarity or causality.</p>
<h2>Compare and zoom mathematically</h2>
<label>Left <select id="a"></select></label><label>Right <select id="b"></select></label>
<label>Magnification <input type="range" id="zoom" min="1" max="5" value="1"> <output id="factor">1×</output></label>
<div class="compare"><canvas id="ca" width="280" height="280"></canvas><canvas id="cb" width="280" height="280"></canvas></div>
<p>Zooming here reveals more Julia-set mathematics, not new measurements. Download or run the source to navigate the underlying data hierarchy.</p>
<h2>All portraits</h2><div class="grid">__CARDS__</div>
<p>See <a href="https://github.com/H-Lupercal/fractal-atlas/tree/main/results/gallery">methods, sources, and comparisons</a>.</p></main>
<script id="data" type="application/json">__DATA__</script><script>
let datasets=JSON.parse(document.getElementById("data").textContent),
a=document.getElementById("a"),b=document.getElementById("b"),zoom=document.getElementById("zoom");
for(let s of [a,b])for(let i=0;i<datasets.length;i++){let o=document.createElement("option");o.value=i;o.textContent=datasets[i].label;s.appendChild(o)}
b.selectedIndex=Math.min(4,datasets.length-1);
function render(c,cs,scale){
let ctx=c.getContext("2d"),w=c.width,h=c.height,img=ctx.createImageData(w,h);
for(let y=0;y<h;y++)for(let x=0;x<w;x++){
let zr=(x/w*2-1)*1.5/scale,zi=(y/h*2-1)*1.5/scale,n=0;
for(;n<110;n++){let re=zr*zr-zi*zi+cs[0],im=2*zr*zi+cs[1];zr=re;zi=im;if(zr*zr+zi*zi>16)break}
let v=n===110?0:Math.floor(n/110*255),k=(y*w+x)*4;
img.data[k]=v;img.data[k+1]=v;img.data[k+2]=v;img.data[k+3]=255;}
ctx.putImageData(img,0,0)
}
function draw(){let s=Math.pow(2,Number(zoom.value)-1);
document.getElementById("factor").textContent=s+"×";
render(document.getElementById("ca"),datasets[Number(a.value)].c,s);
render(document.getElementById("cb"),datasets[Number(b.value)].c,s)}
a.addEventListener("change",draw);b.addEventListener("change",draw);
zoom.addEventListener("input",draw);draw();
</script></html>"""
    (WEB/"index.html").write_text(template.replace("__CARDS__",cards).replace(
        "__DATA__",json.dumps(rec).replace("<","\\u003c")),encoding="utf8")

def main():
    DEST.mkdir(parents=True,exist_ok=True);WEB.mkdir(parents=True,exist_ok=True)
    entries=[]
    with tempfile.TemporaryDirectory(prefix="fractal_atlas_gallery_") as tmp:
        root=Path(tmp)
        for item in make_inputs(root):
            slug=item["slug"];frame=item["frame"]
            file=root/(slug+".csv");frame.to_csv(file,index=False)
            config={"name":item["label"],"data":{
                "path":str(file),"kind":"features","id_column":"id",
                "normalization":"zscore","missing":"median",
                "similarity":"cosine","max_cells":4000000,
                },"render":{"max_iter":ITER}}
            cfg=root/(slug+".json");cfg.write_text(json.dumps(config))
            output=root/(slug+"_render")
            build(cfg,output)
            dst=DEST/(slug+".png");shutil.copyfile(output/"preview.png",dst)
            summary=json.loads((output/"summary.json").read_text())
            info={"slug":slug,"label":item["label"],"group":item["group"],
                  "source":item["source"],"description":item["description"],
                  "rows":len(frame),"features":len(frame.columns)-1,
                  "julia_c":[summary["parameters"]["c_real"],summary["parameters"]["c_imag"]],
                  "mean_similarity":summary["parameters"]["mean"],
                  "std_similarity":summary["parameters"]["std"],**metrics(dst)}
            entries.append(info)
            print(f'{slug}: {len(frame)} entities, {len(frame.columns)-1} features',flush=True)
    comparisons=[]
    for i,left in enumerate(entries):
        x=np.asarray(Image.open(DEST/(left["slug"]+".png")).convert("L"),dtype=np.int16)
        for right in entries[i+1:]:
            y=np.asarray(Image.open(DEST/(right["slug"]+".png")).convert("L"),dtype=np.int16)
            dif=np.abs(x-y)
            comparisons.append({"left":left["slug"],"right":right["slug"],
                "mean_absolute_grayscale_difference":round(float(dif.mean()),5),
                "pixels_differing_by_10_or_more_pct":round(float((dif>=10).mean()*100),3)})
    with (DEST/"pairwise_visual_differences.csv").open("w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(comparisons[0]))
        writer.writeheader();writer.writerows(comparisons)
    (DEST/"manifest.json").write_text(json.dumps({
        "experiment":"Real biological data gallery v1",
        "code":"scripts/run_gallery.py",
        "preprocessing":"median imputation; z-score; cosine similarity",
        "julia":"z[n+1] = z[n]^2 + c; c=(-0.745+0.17*tanh(mean)) + i*(0.186+0.17*tanh(std))",
        "max_iterations":ITER,"datasets":entries,
        "warning":"Visual differences are not evidence of biological similarity or causal interactions",
    },indent=2))
    contact_sheet(entries)
    site_images=WEB/"images"
    site_images.mkdir(parents=True,exist_ok=True)
    for entry in entries:
        shutil.copyfile(DEST/(entry["slug"]+".png"),site_images/(entry["slug"]+".png"))
    web_gallery(entries)
    lines=["# Real-data fractal gallery","","Only computed images and metadata are stored here, not redistributed raw source datasets.","",
        "## Reproduce","","    python -m pip install -e '.[gallery]'",
        "    python scripts/run_gallery.py","",
        "All portraits use the same preprocessing and quadratic Julia parameters. Results are illustrative visual comparisons, not biological inference.",
        "",
        "| Group | Dataset | Entities | Features | Mean gray | Entropy (bits) | Source |",
        "|---|---|---:|---:|---:|---:|---|"]
    for r in entries:
        lines.append(f'| {r["group"]} | [{r["label"]}]({r["slug"]}.png) | {r["rows"]} | {r["features"]} | {r["mean_gray"]} | {r["entropy_bits"]} | [Source]({r["source"]}) |')
    lines.extend(["",
      "## Interpretation","",
      "The image differences in pairwise_visual_differences.csv compare pixel intensities with fixed iteration count and palette. These are **not** biological comparisons.",
      "Input types include plant flower morphology, plant RNA expression, human nucleus morphology, mouse cortex protein abundance, bacterial iModulon activity, and viral nucleotide-window composition.",
      "Source matrices may have different numbers of features and records. Viral genomes are not cells. Feature similarity is not experimentally measured genetic interaction.",
      "The square Julia-set previews are only single-region portraits; the generic Fractal Atlas CLI supports full data and mathematical navigation.",
      "",
      "See [interactive gallery](../../docs/gallery/index.html), [montage](contact_sheet.png), [manifest](manifest.json), and [pairwise metrics](pairwise_visual_differences.csv)."])
    (DEST/"README.md").write_text("\n".join(lines)+"\n")
    print(f"Complete: {len(entries)} portraits; {len(comparisons)} pairwise comparisons.",flush=True)

if __name__=="__main__":main()
