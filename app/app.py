import os
import streamlit as st
import pandas as pd
import csv
import collections
import random
st.set_page_config(layout="wide")

root = os.path.dirname(os.path.abspath(__file__))

# read data

@st.cache
def load_screening():
    df = pd.read_csv(os.path.join(root, "../data/screening.tsv"), sep="\t")
    return df

@st.cache
def load_screening_hits():
    db = pd.read_csv(os.path.join(root, "../data/screening_hits.tsv"), sep="\t")
    return db

@st.cache
def load_human_proteome():
    human_proteome = pd.read_csv(os.path.join(root, "../data/human_proteome_with_gene_names.tab"), sep="\t")
    return human_proteome

@st.cache
def load_hek_proteome():
    hek_proteome = []
    with open(os.path.join(root, "../data/hek293t_core.tsv"), "r") as f:
        reader = csv.reader(f)
        for r in reader:
            hek_proteome += [r[0]]
    hek_proteome = set(hek_proteome)
    return hek_proteome

def load_pid2name_primary():
    pid2name = {}
    with open(os.path.join(root, "../data/pid2name_primary.tsv"), "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            pid2name[r[0]] = r[1]
    return pid2name

hek_proteome = load_hek_proteome()
human_proteome = load_human_proteome()
df = load_screening()
db = load_screening_hits()

pid2name_primary = load_pid2name_primary()
for r in human_proteome[["Entry", "Gene names  (primary )"]].values:
    if r[0] in pid2name_primary:
        continue
    if str(r[1]) == "nan":
        pid2name_primary[r[0]] = r[0]
    else:
        pid2name_primary[r[0]] = r[1]
    
any2pid = {}
for k,v in pid2name_primary.items():
    any2pid[v] = k
for r in human_proteome[["Entry", "Gene names"]].values:
    g = str(r[1])
    any2pid[r[0]] = r[0]
    if g == "nan":
        continue
    for x in r[1].split(" "):
        any2pid[x] = r[0]

pid2fid = collections.defaultdict(list)
for r in db[["Accession", "FragID"]].values:
    pid2fid[r[0]] += [r[1]]

frequent_hitters = set()
normal_hitters = set()
specific_hitters = set()
for k,v in pid2fid.items():
    if len(v) >= 40:
        frequent_hitters.update([k])
        continue
    if len(v) >= 10:
        normal_hitters.update([k])
        continue
    specific_hitters.update([k])

st.title("Ligand Discovery Protein Profile Explorer")

example = list(db[db["FragID"] == "C391"]["Accession"])
example = sorted([pid2name_primary[x] for x in example])


options = set([k for k, _ in any2pid.items()])

user_input = st.multiselect(label="Input proteins", options = options, default=random.sample(example, 20))
user_pids = {}

for i in user_input:
    user_pids[i] = any2pid[i]

st.file_uploader(label="Or load file")

columns = st.columns([2, 2, 2, 1, 1])

done = set()

col = columns[0]
col.subheader("Frequent hitters")
R = []
for r in user_input:
    pid = user_pids[r]
    if pid in frequent_hitters:
        R += [[r, pid2fid[pid]]]
        done.update([r])
df = pd.DataFrame(R, columns=["Entry", "Fragments"])
col.dataframe(df)

col = columns[1]
col.subheader("Medium specificity")
R = []
for r in user_input:
    pid = user_pids[r]
    if pid in normal_hitters:
        R += [[r, len(pid2fid[pid]), pid2fid[pid]]]
        done.update([r])
df = pd.DataFrame(R, columns=["Entry", "Hits", "Fragments"])
col.dataframe(df)

col = columns[2]
col.subheader("High specificity")
R = []
for r in user_input:
    pid = user_pids[r]
    if pid in specific_hitters:
        R += [[r, len(pid2fid[pid]), pid2fid[pid]]]
        done.update([r])
df = pd.DataFrame(R, columns=["Entry", "Hits", "Fragments"])
col.dataframe(df)

col = columns[3]
col.subheader("Never hitted")
R = []
for r in user_input:
    if r in done:
        continue
    pid = user_pids[r]
    if pid in hek_proteome:
        R += [[r, len(pid2fid[pid]), pid2fid[pid]]]
        done.update([r])
df = pd.DataFrame(R, columns=["Entry", "Hits", "Fragments"])
col.dataframe(df["Entry"])

col = columns[4]
col.subheader("Not in HEK293T")
R = []
for r in user_input:
    if r in done:
        continue
    pid = user_pids[r]
    if pid in human_proteome:
        R += [[r, len(pid2fid[pid]), pid2fid[pid]]]
df = pd.DataFrame(R, columns=["Entry", "Hits", "Fragments"])
col.dataframe(df["Entry"])

@st.cache
def convert_df(df):
    return df.to_csv().encode('utf-8')

data = convert_df(df)
st.download_button(label="Download search results", data=data, file_name="ligand_discovery_search_results.csv", mime="text/csv")