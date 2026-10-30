import json
from pathlib import Path
import numpy as np
import pandas as pd
from fractal_atlas.cli import prepare, region_stats, julia, build

def test_feature_similarity(tmp_path):
    pd.DataFrame({"id":["A","B","C"],"x":[1,1,-1],"y":[0,0,0]}).to_csv(tmp_path/"a.csv",index=False)
    ids,scores,features,cols,kind=prepare({"data":{"path":"a.csv","normalization":"none"}},tmp_path/"config.yml")
    assert ids==["A","B","C"]
    assert kind=="features"
    assert np.isclose(scores[0,1],1)
    assert np.isclose(scores[0,2],-1)

def test_julia_and_region():
    scores=np.ones((4,4),dtype=np.float32)
    s=region_stats(scores,(0,4,0,4))
    assert s["coverage"]==1
    assert julia(complex(s["c_real"],s["c_imag"]),size=40,max_iter=30).shape==(40,40)

def test_build(tmp_path):
    pd.DataFrame({"id":["a","b"],"x":[0.1,0.2],"y":[2,3]}).to_csv(tmp_path/"data.csv",index=False)
    (tmp_path/"config.json").write_text(json.dumps({"data":{"path":"data.csv"}}))
    build(tmp_path/"config.json",tmp_path/"out")
    assert (tmp_path/"out"/"index.html").is_file()
    assert (tmp_path/"out"/"preview.png").is_file()
    assert (tmp_path/"out"/"summary.json").is_file()
