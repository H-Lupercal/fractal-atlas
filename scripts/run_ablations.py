"""Controlled input-size variations for real biological fractal portraits.

Experiments vary either the feature count at fixed sample count (mouse cortex)
or the entity count at fixed feature schema (maize expression). They do not
assert mathematical or biological information increases monotonically.
"""
from __future__ import annotations
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps

from fractal_atlas.cli import build
from run_gallery import ROOT, download, metrics, SOURCES

OUT=ROOT/"results"/"ablations"
SEED=1741

def render(frame, slug, label, tmp, details):
    src=tmp/(slug+".csv");frame.to_csv(src,index=False)
    cfg={"name":label,"data":{"path":str(src),"kind":"features","id_column":"id",
          "normalization":"zscore","missing":"median","similarity":"cosine",
          "max_cells":4000000},"render":{"max_iter":110}}
    config=tmp/(slug+".json");config.write_text(json.dumps(cfg))
    output=tmp/(slug+"_out")
    build(config,output)
    shutil.copyfile(output/"preview.png",OUT/(slug+".png"))
    meta=json.loads((output/"summary.json").read_text())
    rec={"slug":slug,"label":label,"entities":len(frame),
         "features":len(frame.columns)-1,"julia_parameters":
         [meta["parameters"]["c_real"],meta["parameters"]["c_imag"]],
         "mean_similarity":meta["parameters"]["mean"],**details,
         **metrics(OUT/(slug+".png"))}
    print(slug,rec["entities"],rec["features"],flush=True)
    return rec

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    results=[]
    with tempfile.TemporaryDirectory(prefix="atlas_ablation_") as d:
        tmp=Path(d)
        mouse=pd.read_csv(download("mouse",tmp))
        proteins=sorted(c for c in mouse if c.endswith("_N"))
        control=mouse[mouse["Genotype"].astype(str).eq("Control")].head(180)
        if len(control)<180:raise ValueError("Insufficient mouse controls")
        for k in [8,24,55,len(proteins)]:
            frame=control[["MouseID"]+proteins[:k]].rename(columns={"MouseID":"id"}).copy()
            results.append(render(frame,f"mouse_features_{k:02d}",
                 f"Mouse cortex: {k} features, 180 observations",tmp,
                 {"experiment":"feature-count","source":SOURCES["mouse"],
                 "selection":"first k sorted protein feature names; same 180 mouse observations"}))

        maize=pd.read_csv(download("maize",tmp))
        cols=[c for c in maize if c!="GeneID"]
        maize[cols]=maize[cols].apply(pd.to_numeric,errors="coerce")
        maize["_mean"]=maize[cols].mean(axis=1)
        viable=maize[(maize["_mean"]>0.02)&np.isfinite(maize["_mean"])].drop_duplicates("GeneID")
        high=viable.sort_values(["_mean","GeneID"]).tail(180)
        order=np.random.default_rng(SEED).permutation(len(high))
        for n in [64,128,180]:
            sub=high.iloc[np.sort(order[:n])]
            frame=sub[["GeneID"]+cols].rename(columns={"GeneID":"id"}).copy()
            frame[cols]=np.log1p(frame[cols].clip(lower=0))
            results.append(render(frame,f"maize_genes_{n:03d}",
                 f"Maize: {n} genes, six expression features",tmp,
                 {"experiment":"entity-count","source":SOURCES["maize"],
                 "selection":"nested seeded subsets of the same 180-gene high-expression cohort"}))

    tile_w,tile_h,ncols=330,360,4
    canvas=Image.new("RGB",(tile_w*ncols,tile_h*2),"#101924")
    dr=ImageDraw.Draw(canvas)
    for k,x in enumerate(results):
        xp=(k%ncols)*tile_w;yp=(k//ncols)*tile_h
        im=ImageOps.fit(Image.open(OUT/(x["slug"]+".png")).convert("RGB"),(310,310))
        canvas.paste(im,(xp+10,yp+7))
        dr.text((xp+10,yp+321),x["label"][:44],fill="white")
    canvas.save(OUT/"contact_sheet.png",optimize=True)
    (OUT/"manifest.json").write_text(json.dumps({
        "title":"Controlled real-data input-size fractal experiments",
        "seed":SEED,"preprocessing":"median, zscore, cosine similarity",
        "fractal":"quadratic Julia, same c mapping and iteration count as gallery",
        "datasets":results,
        "warning":"Changes in the images do not prove additional visual complexity or biological insight."
    },indent=2))
    lines=["# Input-size experiments","","All images use real data and the same Julia rendering rule.",
       "For the **feature-count** test we hold the 180 mouse cortex samples fixed and introduce more measured protein features.",
       "For the **entity-count** test we hold six maize expression features fixed and grow nested, reproducible gene subsets.",
       "The results test visual sensitivity, **not** improved information recovery or biology.","",
       "| Experiment | Portrait | Entities | Features | Gray entropy (bits) | Source |",
       "|---|---|---:|---:|---:|---|"]
    for x in results:
        lines.append(f'| {x["experiment"]} | [{x["label"]}]({x["slug"]}.png) | {x["entities"]} | {x["features"]} | {x["entropy_bits"]} | [Source]({x["source"]}) |')
    lines.extend(["","## Reproduce","","    python -m pip install -e '.[gallery]'",
                  "    python scripts/run_ablations.py","",
                  "Only derived images and metadata are committed, never the third-party raw datasets.",
                  "Feature selection is deterministic but not a controlled biological experiment.",
                  "The order of rows and feature selection can affect results, as can normalization."])
    (OUT/"README.md").write_text("\n".join(lines)+"\n")
    print("Produced",len(results),"controlled portraits")

if __name__=="__main__":main()
