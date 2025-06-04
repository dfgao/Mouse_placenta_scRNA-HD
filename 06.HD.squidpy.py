import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import cv2
import os
import scipy.sparse
import squidpy as sq

import bin2cell as b2c

import seaborn as sns
from pathlib import Path

path008 = "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A3/outs/binned_outputs/square_008um/"
path002 = "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A3/outs/binned_outputs/square_002um/"
source_image_path = "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/01.data/X101SC24102043-Z01-J001/Rawdata/A3/A3.jpg"

bdata = b2c.read_visium(path008, source_image_path = source_image_path)
bdata.var_names_make_unique()

path008 = "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A4/outs/binned_outputs/square_008um/"
path002 = "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A4/outs/binned_outputs/square_002um/"
source_image_path = "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/01.data/X101SC24102043-Z01-J001/Rawdata/A4/A4.jpg"

ddata = b2c.read_visium(path008, source_image_path = source_image_path)
ddata.var_names_make_unique()

spatial_data = bdata.obsm["spatial"]

plt.figure(figsize=(8, 6))
plt.scatter(spatial_data[:, 0], spatial_data[:, 1], alpha=0.5,s=0.2)
plt.title("Scatter Plot of Spatial 8um Coordinates")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

mask_young = bdata.obsm["spatial"][:, 0] < 13500  # 左侧区域 young
mask_age = ~mask_young  # 右侧区域 age

bdata_young = bdata[mask_young].copy()  # 左侧区域
bdata_age = bdata[mask_age].copy()  # 右侧区域

print(f"Young cells: {bdata_young.n_obs}")
print(f"Age cells: {bdata_age.n_obs}")

bdata_young.write_h5ad('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/A3_8um.young.h5ad')
bdata_age.write_h5ad('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/A3_8um.age.h5ad')

spatial_data = ddata.obsm["spatial"]

plt.figure(figsize=(8, 6))
plt.scatter(spatial_data[:, 0], spatial_data[:, 1], alpha=0.5,s=0.2)
plt.title("Scatter Plot of Spatial 8um Coordinates")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

mask_age = ddata.obsm["spatial"][:, 0] < 14000  # 左侧区域 age
mask_ko = ~mask_age  # 右侧区域 ko

ddata_ko = ddata[mask_ko].copy()  # 左侧区域
ddata_age = ddata[mask_age].copy()  # 右侧区域

print(f"KO cells: {ddata_ko.n_obs}")
print(f"Age cells: {ddata_age.n_obs}")

ddata_ko.write_h5ad('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/A4_8um.ko.h5ad')
ddata_age.write_h5ad('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/A4_8um.age.h5ad')

import pandas as pd
age_jzdb = pd.read_csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/age3_jz_db.csv')

barcode = age_jzdb['Barcode'].tolist()
bdata_age_jzdb = bdata_age[bdata_age.obs_names.isin(barcode), :]
bdata_age_jzdb

spatial_data = bdata_age_jzdb.obsm["spatial"]

plt.figure(figsize=(3, 6))
plt.scatter(spatial_data[:, 0], spatial_data[:, 1], alpha=0.5,s=0.2)
plt.title("Scatter Plot of Spatial 8um Coordinates")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

young_jzdb = pd.read_csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/young1_jz_db.csv')
barcode = young_jzdb['Barcode'].tolist()
bdata_young_jzdb = bdata_young[bdata_young.obs_names.isin(barcode), :]
bdata_young_jzdb

spatial_data = bdata_young_jzdb.obsm["spatial"]

plt.figure(figsize=(3, 6))
plt.scatter(spatial_data[:, 0], spatial_data[:, 1], alpha=0.5,s=0.2)
plt.title("Scatter Plot of Spatial 8um Coordinates")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

ko_jzdb = pd.read_csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/ko2_jz_db.csv')
barcode = ko_jzdb['Barcode'].tolist()
ddata_ko_jzdb = ddata_ko[ddata_ko.obs_names.isin(barcode), :]
ddata_ko_jzdb

spatial_data = ddata_ko_jzdb.obsm["spatial"]

plt.figure(figsize=(3, 6))
plt.scatter(spatial_data[:, 0], spatial_data[:, 1], alpha=0.5,s=0.2)
plt.title("Scatter Plot of Spatial 8um Coordinates")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

age_5ct_macrop = pd.read_csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.csv')
barcode = age_5ct_macrop['Barcode'].tolist()
bdata_age_5ct_macrop = bdata_age[bdata_age.obs_names.isin(barcode), :]
bdata_age_5ct_macrop

age_5ct_macrop = age_5ct_macrop.set_index('Barcode')
bdata_age_5ct_macrop.obs = bdata_age_5ct_macrop.obs.join(age_5ct_macrop, how="left")
print(bdata_age_5ct_macrop.obs.head())

sq.gr.spatial_neighbors(bdata_age_5ct_macrop,coord_type="generic",spatial_key="spatial",delaunay=True)
sq.gr.nhood_enrichment(bdata_age_5ct_macrop,cluster_key="guess.ct")
sq.pl.nhood_enrichment(
    bdata_age_5ct_macrop,cluster_key="guess.ct",
    method="average",cmap="inferno",
    vmin=-50,vmax=100,figsize=(5,5)
)

fig1 = plt.figure()
ax = fig1.add_subplot(111) 
sq.pl.nhood_enrichment(
    bdata_age_5ct_macrop,cluster_key="guess.ct",
    method="average",cmap="inferno",
    vmin=-50,vmax=100,figsize=(5,5),ax=ax
)
fig1.savefig("/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/age.neibor.box.pdf", format="pdf")
plt.show()

fig2 = plt.figure()
ax = fig2.add_subplot(111) 
sc.pl.spatial(bdata_age_5ct_macrop, img_key="hires", color="guess.ct",alpha_img=0.2,ax = ax)
# fig2.savefig("/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/age.RCTD.pdf", format="pdf")

plt.show()

bdata_age_5ct_macrop.obs.columns = bdata_age_5ct_macrop.obs.columns.str.replace('.', '_')

sq.gr.centrality_scores(bdata_age_5ct_macrop, cluster_key="guess_ct")
sq.pl.centrality_scores(bdata_age_5ct_macrop, cluster_key="guess_ct", figsize=(10, 10))

sq.pl.co_occurrence(
    bdata_age_5ct_macrop,
    cluster_key="guess_ct",
    clusters="DSCs",
)
plt.savefig('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/age.co_occurrence.svg', format='svg')
fig = sg.fromfile('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/age.co_occurrence.svg')
fig.save('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/age.co_occurrence.pdf')

SVG('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/age.co_occurrence.pdf')

sq.pl.co_occurrence(
    bdata_age_5ct_macrop,
    cluster_key="guess_ct",
    clusters="DSCs",
)

young_5ct_macrop = pd.read_csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.csv')
barcode = young_5ct_macrop['Barcode'].tolist()
bdata_young_5ct_macrop = bdata_young[bdata_young.obs_names.isin(barcode), :]
bdata_young_5ct_macrop
young_5ct_macrop = young_5ct_macrop.set_index('Barcode')
bdata_young_5ct_macrop.obs = bdata_young_5ct_macrop.obs.join(young_5ct_macrop, how="left")
print(bdata_young_5ct_macrop.obs.head())

sq.gr.spatial_neighbors(bdata_young_5ct_macrop,coord_type="generic",spatial_key="spatial",delaunay=True)
sq.gr.nhood_enrichment(bdata_young_5ct_macrop,cluster_key="guess.ct")
sq.pl.nhood_enrichment(
    bdata_age_5ct_macrop,cluster_key="guess.ct",
    method="ward",cmap="inferno",
    vmin=-100,vmax=100,figsize=(5,5)
)
sc.pl.spatial(bdata_young_5ct_macrop, img_key="hires", color="guess.ct",alpha_img=0.1,size=5,palette=bdata_young_5ct_macrop.uns['guess.ct_colors'])
bdata_young_5ct_macrop.obs.columns = bdata_young_5ct_macrop.obs.columns.str.replace('.', '_')

sq.gr.co_occurrence(
    bdata_young_5ct_macrop,
    cluster_key="guess_ct",
)
fig7 = plt.figure()
sq.pl.co_occurrence(
    bdata_young_5ct_macrop,
    cluster_key="guess_ct",
    clusters="DSCs",
)
plt.savefig('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/young.co_occurrence.svg', format='svg')
fig = sg.fromfile('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/young.co_occurrence.svg')
fig.save('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/young.co_occurrence.pdf')

SVG('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/young.co_occurrence.pdf')

ko_5ct_macrop = pd.read_csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.csv')
barcode = ko_5ct_macrop['Barcode'].tolist()
ddata_ko_5ct_macrop = ddata_ko[ddata_ko.obs_names.isin(barcode), :]
ddata_ko_5ct_macrop

fig8 = plt.figure()
ax = fig8.add_subplot(111) 
sc.pl.spatial(ddata_ko_5ct_macrop, img_key="hires", color="guess_ct",alpha_img=0.2,ax = ax)
fig8.savefig("/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/KO.RCTD.pdf", format="pdf")

plt.show()

fig4 = plt.figure()
ax = fig4.add_subplot(111) 
sq.pl.nhood_enrichment(
    ddata_ko_5ct_macrop,cluster_key="guess_ct",
    method="average",cmap="inferno",
    vmin=-50,vmax=100,figsize=(5,5),ax=ax
)
fig4.savefig("/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/ko.neibor.box.pdf", format="pdf")
plt.show()
sq.pl.co_occurrence(
    ddata_ko_5ct_macrop,
    cluster_key="guess_ct",
    clusters="DSCs",
)
plt.savefig('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/ko.co_occurrence.svg', format='svg')
fig = sg.fromfile('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/ko.co_occurrence.svg')
fig.save('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/ko.co_occurrence.pdf')

SVG('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/03.py.analysis/ko.co_occurrence.pdf')
