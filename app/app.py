import os
import streamlit as st
import pandas as pd
import csv
import collections
import joblib
st.set_page_config(layout="wide")

root = os.path.dirname(os.path.abspath(__file__))

FREQUENT_CUTOFF = 40
MEDIUM_CUTOFF = 10

# read data

@st.cache(suppress_st_warning=True)
def load_screening_hits():
    db = pd.read_csv(os.path.join(root, "../data/screening_hits.tsv"), sep="\t")
    return db

@st.cache(suppress_st_warning=True)
def load_human_proteome():
    human_proteome = pd.read_csv(os.path.join(root, "../data/human_proteome_with_gene_names.tab"), sep="\t")
    return human_proteome

@st.cache(suppress_st_warning=True)
def load_hek_proteome():
    hek_proteome = []
    with open(os.path.join(root, "../data/hek293t_core.tsv"), "r") as f:
        reader = csv.reader(f)
        for r in reader:
            hek_proteome += [r[0]]
    hek_proteome = set(hek_proteome)
    return hek_proteome

@st.cache(suppress_st_warning=True)
def load_pid2name_primary():
    return joblib.load(os.path.join(root, "../data/pid2name_primary.joblib"))

@st.cache(suppress_st_warning=True)
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')

@st.cache(suppress_st_warning=True)
def convert_df_no_header(df):
    return df.to_csv(index=False, header=False).encode('utf-8')

@st.cache(suppress_st_warning=True)
def example_input_load():
    pids = []
    with open(os.path.join(root, "../data/example_input.csv"), "r") as f:
        reader = csv.reader(f)
        for r in reader:
            pids += [r[0]]
    return pids

db = load_screening_hits()
hek_proteome = load_hek_proteome()
pid2name_primary = load_pid2name_primary()
human_proteome = set(pid2name_primary.keys())
example_input = example_input_load()

any2pid = {}
for k,v in pid2name_primary.items():
    any2pid[v] = k
    any2pid[k] = k

pid2fid = collections.defaultdict(list)
fid2pid = collections.defaultdict(list)
for r in db[["Accession", "FragID"]].values:
    pid2fid[r[0]] += [r[1]]
    fid2pid[r[1]] += [r[0]]

frequent_hitters = set()
normal_hitters = set()
specific_hitters = set()
for k,v in pid2fid.items():
    if len(v) >= FREQUENT_CUTOFF:
        frequent_hitters.update([k])
        continue
    if len(v) >= MEDIUM_CUTOFF:
        normal_hitters.update([k])
        continue
    specific_hitters.update([k])

options = sorted([x for k,v in pid2name_primary.items() for x in [k,v]])

# layout

st.title("Ligand Discovery Protein Profile Explorer")
st.write("Welcome to the Ligand Discovery Protein Profile Explorer! We have screened 412 fully-functionalized small molecule fragments in HEK293t cells. For {0} of the fragments, we found at least one protein enriched. In total, we enriched {1} proteins at least once. Query your protein sets of interest and explore them in light of our dataset!".format(len(fid2pid), len(pid2fid)))

cols = st.columns([2,1,1])

col = cols[0]

manual_input = col.multiselect(label="Input proteins manually", options = [""] + sorted(options), default=[], help="Select proteins by UniProt Accession code or Gene Symbol")
user_pids = {}
user_input = []
for i in manual_input:
    user_pids[i] = any2pid[i]
    user_input += [i]

col = cols[1]

fids = sorted(set(db["FragID"]))
fid_input = col.selectbox(label="Select pre-screened fragment by identifier", options = [""] + fids, help="Select an already profiled fragment in our primary screening. Use the Ligand Discovery fragment identifier")
if fid_input != "":
    user_input = fid2pid[fid_input]
    user_pids = dict((r,r) for r in user_input)

col = cols[2]

example_file = db
file_input = col.file_uploader(label="Upload a file", help="Provide a file containing one UniProt Accession code or Gene Symbol per row.")
if file_input:
    user_input = list(pd.read_csv(file_input, header=None)[0])
    for i in user_input:
        user_pids[i] = any2pid[i]

col.download_button(label="Download example file", data=convert_df_no_header(pd.DataFrame({"uniprot_ac": example_input})), file_name="protein_profile_example.csv", mime="text/csv")

# checks

if not manual_input:
    manual_input = None

if not fid_input:
    fid_input = None

if not file_input:
    file_input = None

if not manual_input and not file_input and not fid_input:
    st.info("Use any of the options above to explore a protein profile...")
    query_is_available = False
else:
    c = 0
    for x in [manual_input, fid_input, file_input]:
        if x is not None:
            c += 1
    if c > 1:
        st.error("More than one input type has been provided! Please only choose one of the options, i.e. input proteins manually, or select a pre-screened fragment, or upload a file. Refresh this window to get started again.")
        query_is_available = False
    else:
        query_is_available = True


def serialize_s(cat, r):
    s = [cat] + r[:-1] + [" ".join(r[-1])]
    return s

if query_is_available:
    columns = st.columns([2, 2, 2, 1, 1])

    done = set()

    col = columns[0]
    cat_name = "Frequently enriched"
    col.subheader(cat_name)
    S = []
    R = []
    for r in user_input:
        pid = user_pids[r]
        if pid in frequent_hitters:
            R += [[pid, pid2name_primary[pid], len(pid2fid[pid]), pid2fid[pid]]]
            S += [serialize_s(cat_name, R[-1])]
            done.update([r])
    df = pd.DataFrame(R, columns=["UniProt", "GeneName", "Hits", "Fragments"])
    col.metric(label="Counts", value=df.shape[0])
    col.dataframe(df, use_container_width=True)

    col = columns[1]
    cat_name = "Medium specificity"
    col.subheader(cat_name)
    R = []
    for r in user_input:
        pid = user_pids[r]
        if pid in normal_hitters:
            R += [[pid, pid2name_primary[pid], len(pid2fid[pid]), sorted(pid2fid[pid])]]
            S += [serialize_s(cat_name, R[-1])]
            done.update([r])
    df = pd.DataFrame(R, columns=["UniProt", "GeneName", "Hits", "Fragments"])
    col.metric(label="Counts", value=df.shape[0])
    col.dataframe(df, use_container_width=True)

    col = columns[2]
    cat_name = "High specificity"
    col.subheader(cat_name)
    R = []
    for r in user_input:
        pid = user_pids[r]
        if pid in specific_hitters:
            R += [[pid, pid2name_primary[pid], len(pid2fid[pid]), sorted(pid2fid[pid])]]
            S += [serialize_s(cat_name, R[-1])]
            done.update([r])
    df = pd.DataFrame(R, columns=["UniProt", "GeneName", "Hits", "Fragments"])
    col.metric(label="Counts", value=df.shape[0])
    col.dataframe(df, use_container_width=True)

    col = columns[3]
    cat_name = "Never enriched"
    col.subheader(cat_name)
    R = []
    for r in user_input:
        if r in done:
            continue
        pid = user_pids[r]
        if pid in hek_proteome:
            R += [[pid, pid2name_primary[pid], len(pid2fid[pid]), sorted(pid2fid[pid])]]
            S += [serialize_s(cat_name, R[-1])]
            done.update([r])
    df = pd.DataFrame(R, columns=["UniProt", "GeneName", "Hits", "Fragments"])
    col.metric(label="Counts", value=df.shape[0])
    col.dataframe(df[["UniProt", "GeneName"]], use_container_width=True)

    col = columns[4]
    cat_name = "Not in HEK293t"
    col.subheader(cat_name)
    R = []
    for r in user_input:
        if r in done:
            continue
        pid = user_pids[r]
        if pid in human_proteome:
            fids_ = sorted(pid2fid[pid])
            R += [[pid, pid2name_primary[pid], len(pid2fid[pid]), fids_]]
            S += [serialize_s(cat_name, R[-1])]
    df = pd.DataFrame(R, columns=["UniProt", "GeneName", "Hits", "Fragments"])
    col.metric(label="Counts", value=df.shape[0])
    col.dataframe(df[["UniProt", "GeneName"]], use_container_width=True)
    
    data = pd.DataFrame(S, columns = ["Category", "UniProt", "GeneName", "Hits", "Fragments"])
    data = data.sort_values(by=["Hits", "GeneName", "Category"], ascending=[False, True, True]).reset_index(drop=True)
    data = convert_df(data)
    st.download_button(label="Download search results", data=data, file_name="ligand_discovery_search_results.csv", mime="text/csv")