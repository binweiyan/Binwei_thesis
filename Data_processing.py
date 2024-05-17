gene_2_gene = True
gene_2_cell = False
Merfish_dir = '/home/bineva/Merfish/'
Starmap_dir = '/home/bineva/STARmapPLUS.Nature2023/'
Slideseq_dir = '/home/bineva/imputation/'


# In[2]:


from collections import defaultdict
import pickle
import pandas as pd
import os
def process_Merfish_df(Merfish_dir = Merfish_dir):
    os.chdir(Merfish_dir)
    df = {}
    metadata = pd.read_pickle('metadata.pkl')
    path = 'data/'
    for i in os.listdir(path):

        if 'pkl' in i and path + i.split('.')[0] in metadata['label'].values:
            
            df[path + i.split('.')[0]] = pd.read_pickle(path + i)
            print(i, df[path + i.split('.')[0]].shape)
            if df[path + i.split('.')[0]].shape[0] ==0:
                del df[path + i.split('.')[0]]
            else:
                print(i, df[path + i.split('.')[0]])
    label = {}
    for i in range(len(metadata)):
        label[metadata['NAME'][i]] = metadata['label'][i]

    groups = defaultdict(list)
    for k, v in label.items():
        groups[v].append(k)

    num_check = defaultdict(dict)
    for i in groups:
        index = 0
        groups[i] = [j for j in groups[i]]
        groups[i].sort()
        for j in groups[i]:
            num_check[i][j] = index
            index += 1

    edges_2_from = defaultdict(list)
    edges_2_to = defaultdict(list)
    edges_dict = defaultdict(lambda: defaultdict(list))

    for i in df:
        for j in df[i].index:

            edges_2_from[i].append(num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_2_to[i].append(num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_dict[i][num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])]].append(num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_dict[i][num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])]].append(num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])])

    g = {i:i for i in df}
    with open('test_keys.pkl', 'rb') as f:
        test_keys = pickle.load(f)
    g_test = {key: g[key] for key in test_keys}
    for i in test_keys:
        del g[i]

    return df, edges_2_from, edges_2_to, edges_dict, g, g_test, groups

def process_Slideseq_df(Slideseq_dir = Slideseq_dir):
    os.chdir(Slideseq_dir)
    metadata = pd.read_pickle('metadata.pkl')
    df = {}
    path = 'data/'
    for i in os.listdir(path):
        if 'pkl' in i and path + i.split('.')[0] in metadata['label'].values:
            df[path + i.split('.')[0]] = pd.read_pickle(path + i)
            print(i, df[path + i.split('.')[0]].shape)
            if df[path + i.split('.')[0]].shape[0] ==0:
                del df[path + i.split('.')[0]]
    label = {}
    for i in range(len(metadata)):
        label[metadata['NAME'][i]] = metadata['label'][i]
    groups = defaultdict(list)
    for k, v in label.items():
        groups[v].append(k)
    num_check = defaultdict(dict)
    for i in groups:
        index = 0
        groups[i] = [j for j in groups[i]]
        groups[i].sort()
        for j in groups[i]:
            num_check[i][j] = index
            index += 1
    edges_2_from = defaultdict(list)
    edges_2_to = defaultdict(list)
    edges_dict = defaultdict(lambda: defaultdict(list))
    for i in df:
        for j in df[i].index:
            edges_2_from[i].append(num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_2_to[i].append(num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_dict[i][num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])]].append(num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_dict[i][num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])]].append(num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])])
    g = {i:i for i in df}
    with open('test_keys.pkl', 'rb') as f:
        test_keys = pickle.load(f)
    g_test = {key: g[key] for key in test_keys}
    for i in test_keys:
        del g[i]
    return df, edges_2_from, edges_2_to, edges_dict, g, g_test, groups

def process_Starmap_df(Starmap_dir = Starmap_dir):
    os.chdir(Starmap_dir)
    metadata = pd.read_pickle('metadata.pkl')
    df = {}
    for i in os.listdir('.'):
        if 'pkl' in i and i.split('.')[0] in metadata['label'].values:
            df[i.split('.')[0]] = pd.read_pickle(i)
            if df[i.split('.')[0]].shape[1] !=5:
                del df[i.split('.')[0]]
            else:
                print(i, df[i.split('.')[0]].shape)
    label = {}
    for i in range(len(metadata)):
        label[metadata['NAME'][i]] = metadata['label'][i]

    groups = defaultdict(list)
    for k, v in label.items():
        groups[v].append(k)

    num_check = defaultdict(dict)
    for i in groups:
        index = 0
        groups[i] = [j for j in groups[i]]
        groups[i].sort()
        for j in groups[i]:
            num_check[i][j] = index
            index += 1

    edges_2_from = defaultdict(list)
    edges_2_to = defaultdict(list)
    edges_dict = defaultdict(lambda: defaultdict(list))

    for i in df:
        for j in df[i].index:

            edges_2_from[i].append(num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_2_to[i].append(num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_dict[i][num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])]].append(num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])])
            edges_dict[i][num_check[i][max(df[i]['cell1'][j],df[i]['cell2'][j])]].append(num_check[i][min(df[i]['cell1'][j],df[i]['cell2'][j])])

    g = {i:i for i in df}
    with open('test_keys.pkl', 'rb') as f:
        test_keys = pickle.load(f)
    g_test = {key: g[key] for key in test_keys}
    for i in test_keys:
        del g[i]
    return df, edges_2_from, edges_2_to, edges_dict, g, g_test, groups

