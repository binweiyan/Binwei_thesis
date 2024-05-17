#!/usr/bin/env python
# coding: utf-8

# In[1]:


#read all txt files in the directory
import glob
import os
import pandas as pd

device = 'cuda:0'

gene_2_gene = True
gene_2_cell = False
Merfish_dir = '/home/bineva/Merfish/'
Starmap_dir = '/home/bineva/STARmapPLUS.Nature2023/'
Slideseq_dir = '/home/bineva/imputation/'


# In[2]:


from collections import defaultdict
import pickle
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


# In[3]:


metadata_Merfish = pd.read_pickle(Merfish_dir + 'metadata.pkl')
metadata_Starmap = pd.read_pickle(Starmap_dir + 'metadata.pkl')
metadata_Slideseq = pd.read_pickle(Slideseq_dir + 'metadata.pkl')
Merfish_df, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, Merfish_g, Merfish_g_test, Merfish_groups = process_Merfish_df()
Starmap_df, edges_2_from_Starmap, edges_2_to_Starmap, edges_dict_Starmap, Starmap_g, Starmap_g_test, Starmap_groups = process_Starmap_df()
Slideseq_df, edges_2_from_Slideseq, edges_2_to_Slideseq, edges_dict_Slideseq, Slideseq_g, Slideseq_g_test, Slideseq_groups = process_Slideseq_df()


# In[4]:



# ----------- 2. create model -------------- #
# build a two-layer GraphSAGE model
import dgl
import torch
import torch.nn as nn
import torch.nn.functional as F


from dgl.nn import SAGEConv
from dgl.nn import HeteroEmbedding
from dgl.nn import HeteroLinear
class GraphSAGE_gc(nn.Module):
    def __init__(self, in_feats, h_feats):
        super(GraphSAGE_gc, self).__init__()
        self.conv1 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(in_feats, h_feats, 'mean',feat_drop=0.8),
        'gene2cell': SAGEConv(in_feats, h_feats, 'mean',feat_drop=0.8),
        'gene2gene': SAGEConv(in_feats, h_feats, 'mean',feat_drop=0.8)
        },aggregate='sum')
        self.conv2 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'gcn'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'gcn'),
        'gene2gene': SAGEConv(h_feats, h_feats, 'gcn')        
        },aggregate='sum')
        self.conv3 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'pool'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'pool'),
        'gene2gene': SAGEConv(h_feats, h_feats, 'pool')
        },aggregate='sum')

        self.Lin1 = HeteroLinear({'cell': h_feats, 'gene': h_feats}, h_feats)

        self.conv4 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5),
        'gene2cell': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5),
        'gene2gene': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5)
        },aggregate='sum')
        self.conv5 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'gcn'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'gcn'),
        'gene2gene': SAGEConv(h_feats, h_feats, 'gcn')
        },aggregate='sum')
        self.conv6 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'pool'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'pool'),
        'gene2gene': SAGEConv(h_feats, h_feats, 'pool')
        },aggregate='sum')

        self.Lin2 = HeteroLinear({'cell': h_feats, 'gene': h_feats}, h_feats)

        # self.conv7 = dgl.nn.HeteroGraphConv({
        # 'cell2gene': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5),
        # 'gene2cell': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5)
        # },aggregate='sum')
        # self.conv8 = dgl.nn.HeteroGraphConv({
        # 'cell2gene': SAGEConv(h_feats, h_feats, 'gcn'),
        # 'gene2cell': SAGEConv(h_feats, h_feats, 'gcn')
        # },aggregate='sum')
        # self.conv9 = dgl.nn.HeteroGraphConv({
        # 'cell2gene': SAGEConv(h_feats, h_feats, 'pool'),
        # 'gene2cell': SAGEConv(h_feats, h_feats, 'pool')
        # },aggregate='sum')
        self.embed = HeteroEmbedding({'gene': in_feats}, in_feats)

    def forward(self, g, in_feat):
        input = {}
        input['gene'] = in_feat['gene']
        h = self.embed(input)
        in_feat['gene'] = h['gene']

        #check nan
        h = self.conv1(g, in_feat, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']},"gene2gene": {'edge_weight': g.edges['gene2gene'].data['w']}})
        h = {k: F.relu(v) for k, v in h.items()}
        h1 = self.conv2(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']},"gene2gene": {'edge_weight': g.edges['gene2gene'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h1 = self.conv3(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']},"gene2gene": {'edge_weight': g.edges['gene2gene'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h = self.Lin1(h)

        h1 = self.conv4(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']},"gene2gene": {'edge_weight': g.edges['gene2gene'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h1 = self.conv5(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']},"gene2gene": {'edge_weight': g.edges['gene2gene'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h1 = self.conv6(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']},"gene2gene": {'edge_weight': g.edges['gene2gene'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h = self.Lin2(h)

        # h1 = self.conv7(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        # h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        # h1 = self.conv8(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        # h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        # h1 = self.conv9(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        # h = {k: h[k] + F.relu(v) for k, v in h1.items()}
        
        return h
    


# In[5]:


import torch
import numpy as np
import scipy.sparse as sp
def batch_negative_edge(g, i, batchsize, edges_2_from, edges_2_to, edges_dict, metadata, groups):
    u, v = torch.tensor(edges_2_from[i] + edges_2_to[i]).int() , torch.tensor(edges_2_to[i] + edges_2_from[i]).int()

    eids = np.arange(len(edges_2_from[i] + edges_2_to[i]))
    eids = np.random.permutation(eids)
    pos_u, pos_v = u[eids], v[eids]

    # Find all negative edges and split them for training and testing
    # if not too many nodes, you can use the following code to generate all negative edges.
    if g.number_of_nodes('cell') < 100:
        adj = sp.coo_matrix((np.ones(len(u)), (u.cpu().numpy(), v.cpu().numpy())), shape = (len(metadata[metadata['label'] == i]), len(metadata[metadata['label'] == i])))

        adj_neg = 1 - adj.todense() - np.eye(len(metadata[metadata['label'] == i]))
        neg_u, neg_v = np.where(adj_neg != 0)

        neg_eids = np.random.choice(len(neg_u), len(edges_2_from[i] + edges_2_to[i]))
        neg_u, neg_v = neg_u[neg_eids], neg_v[neg_eids]
    else:
        neg_u = np.array([])
        neg_v = np.array([])
        #redo the sampling if not enough
        while len(neg_u) < len(edges_2_from[i] + edges_2_to[i]):
            neg_u = np.concatenate((neg_u, np.random.choice(g.number_of_nodes('cell'), len(edges_2_from[i] + edges_2_to[i]) - len(neg_u))))
            neg_v = np.concatenate((neg_v, np.random.choice(g.number_of_nodes('cell'), len(edges_2_from[i] + edges_2_to[i]) - len(neg_v))))
            # remove existing edges, use edges_dict to check
            for j in range(len(neg_u)):
                while neg_v[j] in edges_dict[i][neg_u[j]] or neg_u[j] == neg_v[j]:
                    neg_u[j] = np.random.choice(g.number_of_nodes('cell'))
                    neg_v[j] = np.random.choice(g.number_of_nodes('cell'))

    #divide the train and test sets into batches, and train the model on each batch.
    #the positive batch:
    if batchsize > len(pos_u):
        batchsize = len(pos_u)
    pos_size = batchsize // 2
    neg_size = batchsize - pos_size
    neg_size *= 1
    pos_split_u, pos_split_v = torch.split(pos_u, pos_size), torch.split(pos_v, pos_size)
    #turn neg_u and neg_v into tensors
    neg_u, neg_v = torch.tensor(neg_u), torch.tensor(neg_v)
    neg_split_u, neg_split_v = torch.split(neg_u, neg_size), torch.split(neg_v, neg_size)

    return  pos_split_u, pos_split_v, neg_split_u, neg_split_v


# In[6]:



import dgl.function as fn

class DotPredictor(nn.Module):
    def forward(self, h_u, h_v):
        # h_u is the representation of the source nodes
        # h_v is the representation of the destination nodes
        #return cosine similarity
        return torch.sum(h_u * h_v, dim=1) / (torch.norm(h_u, dim=1) * torch.norm(h_v, dim=1))

class MLPPredictor(nn.Module):
    def __init__(self, h_feats):
        super().__init__()
        self.W1 = nn.Linear(h_feats * 2, h_feats)
        self.W2 = nn.Linear(h_feats, 1)

    def forward(self, h_u, h_v):
        score = torch.cat([h_u, h_v], 1)
        score = F.relu(self.W1(score))
        score = self.W2(score)
        return score.squeeze()

from dgl.nn import EdgePredictor
class EdgePred(nn.Module):
    def __init__(self, in_features, h_features):
        super().__init__()
        self.predictor = EdgePredictor('cos')
        self.W1 = nn.Linear(in_features, h_features)
        self.W2 = nn.Linear(in_features, h_features)
        self.t = nn.Parameter(torch.Tensor([1.0]))
        self.cutoff = 0.5

    def forward(self, g, h):
        u, v = g.edges()
        h_u = self.W1(F.relu(h[u.long()]))
        h_v = self.W2(F.relu(h[v.long()]))
        score = self.act(self.predictor(h_u, h_v))
        return score.squeeze()

    def act(self, score):
        #parameterized sigmoid
        return torch.sigmoid(self.t * score)

    def compute_auc(self, pos_score, neg_score):
        scores = torch.cat([pos_score, neg_score]).cpu().detach().numpy()
        labels = torch.cat(
            [torch.ones(pos_score.shape[0]), torch.zeros(neg_score.shape[0])]).numpy()

        precision, recall, thresholds = precision_recall_curve(labels, scores)
        # Use AUC function to calculate the area under the curve of precision recall curve
        auc_precision_recall = auc(recall, precision)
        
        return accuracy_score(labels, scores > self.cutoff), auc_precision_recall

from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, accuracy_score
def compute_loss(pos_score, neg_score):
    scores = torch.cat([pos_score, neg_score])
    labels = torch.cat([torch.ones(pos_score.shape[0]), torch.zeros(neg_score.shape[0])]).to(device)
    return F.binary_cross_entropy_with_logits(scores, labels)

def compute_auc(pos_score, neg_score):
    if isinstance(pos_score, np.ndarray):
        pos_score = torch.from_numpy(pos_score)
    if isinstance(neg_score, np.ndarray):
        neg_score = torch.from_numpy(neg_score)
    scores = torch.cat([pos_score, neg_score]).cpu().detach().numpy()
    labels = torch.cat(
        [torch.ones(pos_score.shape[0]), torch.zeros(neg_score.shape[0])]).numpy()
    
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    # Use AUC function to calculate the area under the curve of precision recall curve
    auc_precision_recall = auc(recall, precision)
    return accuracy_score(labels, scores > 0.5), auc_precision_recall


# In[7]:



# ----------- 2. create model -------------- #
# build a two-layer GraphSAGE model
from dgl.nn import HeteroEmbedding
from dgl.nn import HeteroLinear
class GraphSAGE(nn.Module):
    def __init__(self, in_feats, h_feats):
        super(GraphSAGE, self).__init__()
        self.conv1 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(in_feats, h_feats, 'mean',feat_drop=0.8),
        'gene2cell': SAGEConv(in_feats, h_feats, 'mean',feat_drop=0.8)
        },aggregate='sum')
        self.conv2 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'gcn'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'gcn')   
        },aggregate='sum')
        self.conv3 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'pool'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'pool')
        },aggregate='sum')

        self.Lin1 = HeteroLinear({'cell': h_feats, 'gene': h_feats}, h_feats)

        self.conv4 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5),
        'gene2cell': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5)
        },aggregate='sum')
        self.conv5 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'gcn'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'gcn')
        },aggregate='sum')
        self.conv6 = dgl.nn.HeteroGraphConv({
        'cell2gene': SAGEConv(h_feats, h_feats, 'pool'),
        'gene2cell': SAGEConv(h_feats, h_feats, 'pool')
        },aggregate='sum')

        self.Lin2 = HeteroLinear({'cell': h_feats, 'gene': h_feats}, h_feats)

        # self.conv7 = dgl.nn.HeteroGraphConv({
        # 'cell2gene': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5),
        # 'gene2cell': SAGEConv(h_feats, h_feats, 'mean', feat_drop=0.5)
        # },aggregate='sum')
        # self.conv8 = dgl.nn.HeteroGraphConv({
        # 'cell2gene': SAGEConv(h_feats, h_feats, 'gcn'),
        # 'gene2cell': SAGEConv(h_feats, h_feats, 'gcn')
        # },aggregate='sum')
        # self.conv9 = dgl.nn.HeteroGraphConv({
        # 'cell2gene': SAGEConv(h_feats, h_feats, 'pool'),
        # 'gene2cell': SAGEConv(h_feats, h_feats, 'pool')
        # },aggregate='sum')
        self.embed = HeteroEmbedding({'gene': in_feats}, in_feats)

    def forward(self, g, in_feat):
        input = {}
        input['gene'] = in_feat['gene']
        h = self.embed(input)
        in_feat['gene'] = h['gene']

        #check nan
        h = self.conv1(g, in_feat, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        h = {k: F.relu(v) for k, v in h.items()}
        h1 = self.conv2(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h1 = self.conv3(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h = self.Lin1(h)

        h1 = self.conv4(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h1 = self.conv5(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h1 = self.conv6(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        h = self.Lin2(h)

        # h1 = self.conv7(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        # h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        # h1 = self.conv8(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        # h = {k: h[k] + F.relu(v) for k, v in h1.items()}

        # h1 = self.conv9(g, h, mod_kwargs={'cell2gene': {'edge_weight': g.edges['cell2gene'].data['w']}, 'gene2cell': {'edge_weight': g.edges['gene2cell'].data['w']}})
        # h = {k: h[k] + F.relu(v) for k, v in h1.items()}
        
        return h
    


# In[8]:


with open('/home/bineva/index/total_genes.pkl', 'rb') as f:
    total_genes = pickle.load(f)
feature_size = len(total_genes)

if gene_2_gene:
    model = GraphSAGE_gc(feature_size, 256).to(device)
else:   
    model = GraphSAGE(feature_size, 256).to(device)
# model.load_state_dict(torch.load('model_hetro2.pt'))
pred = MLPPredictor(256).to(device)


# In[9]:


cell_connected_gene = pd.read_csv('/home/bineva/mouse.genes_link2cells.txt', header = None)
cell_connected_gene[0] = cell_connected_gene[0].str.upper()
pairs = pd.read_csv('/home/bineva/KEGG_2019_Mouse.signaling.gene_gene.txt', delimiter='\t')
pairs1 = pd.read_csv('/home/bineva/ChEA_2022.mouse.gene_gene.txt', delimiter='\t')
pairs = pd.concat([pairs, pairs1], ignore_index=True)
pairs['gene2'] = pairs['gene2'].str.upper()
pairs['gene1'] = pairs['gene1'].str.upper()


# In[10]:


with open('/home/bineva/index/bad.pkl', 'rb') as f:
    bad = pickle.load(f)


# In[30]:


import scanpy
import itertools
expression_new = {}
import warnings
os.chdir('/home/bineva/train_all')
if gene_2_gene:
    file = 'log_hetrogenes_g_link_total.txt'
else:
    file = 'log_hetrogenes_original_total.txt'
resume_True = True
epoch = 50
model_path = ''
if  resume_True:
    if gene_2_gene:
        model.load_state_dict(torch.load(model_path + 'model_hetro_g_link_total_' + str(epoch) + '.pt'))
        pred.load_state_dict(torch.load(model_path + 'pred_hetro_g_link_total_' + str(epoch) + '.pt'))
    else:
        model.load_state_dict(torch.load(model_path + 'model_hetro_original_total_' + str(epoch) + '.pt'))
        pred.load_state_dict(torch.load(model_path + 'pred_hetro_original_total_' + str(epoch) + '.pt'))
else:
    epoch = -1
if os.path.exists(file):
    f = open(file, 'a')
else:
    f = open(file, 'w')


warnings.filterwarnings("ignore")

optimizer = torch.optim.Adam(itertools.chain(model.parameters(), pred.parameters()), lr=1e-4)

#scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.7, patience=5, min_lr=1e-5)

for e in range(epoch+1,300):
    # forward
    model.train()
    pred.train()
    
    # if e == 40:

    #     optimizer = torch.optim.Adam(itertools.chain(model.parameters(), pred.parameters()), lr=1e-5)
    #     model.conv1.mods.cell2gene.feat_drop.p = 0.5
    #     model.conv1.mods.gene2cell.feat_drop.p = 0.5
    #     model.conv4.mods.cell2gene.feat_drop.p = 0
    #     model.conv4.mods.gene2cell.feat_drop.p = 0
    for i in Slideseq_g:
        l = []
        torch.cuda.empty_cache()
        all_pos = np.array([])
        all_neg = np.array([])
        X = scanpy.read_h5ad(Slideseq_dir + i[:i.find('dist10')] + '_total_genes.h5ad')
        
        scanpy.pp.normalize_total(X, target_sum=1e4)
        scanpy.pp.log1p(X)
        expression = X.to_df().T
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]
        matrix = expression[[str(j) for j in Slideseq_groups[i]]]

        for j in cell_connected_gene[0]:
            if j not in expression.index:
                cell_connected_gene = cell_connected_gene[cell_connected_gene[0] != j]
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
            
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)

        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Slideseq, edges_2_to_Slideseq, edges_dict_Slideseq, metadata_Slideseq, Slideseq_groups)
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):

            h = model(graph, graph.ndata['feat'])['cell']
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = compute_loss(pos_score, neg_score)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            l.append(loss.cpu().detach().numpy())
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)), file = f,flush = True)
        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)))

    for i in Merfish_g:
        if i in bad:
            continue
        #load i.h5ad
        name = 'data_transmit/' + i.split('/')[-1]
        X = scanpy.read_h5ad('/home/bineva/Merfish/' + name + '_total_genes.h5ad')
        expression = X.to_df().T
        # expression.index = X.var['gene_symbol']
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]
        # for j in cell_connected_gene[0]:
        #     if j not in expression.index:
        #         cell_connected_gene = cell_connected_gene[cell_connected_gene[0] != j]
     
        # i = 'data/Zhuang-ABCA-1110'
        l = []
        torch.cuda.empty_cache()
        all_pos = np.array([])
        all_neg = np.array([])
        matrix = expression[[str(j) for j in Merfish_groups[i]]]
        
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
            
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish, Merfish_groups)

        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            h = model(graph, graph.ndata['feat'])['cell']
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            #sigmoid of the score
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = compute_loss(pos_score, neg_score)
            optimizer.zero_grad()
            loss.backward()

            optimizer.step()

            l.append(loss.cpu().detach().numpy())
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)), file = f,flush = True)
        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)))

        print(name)
        X = scanpy.read_h5ad('/home/bineva/Merfish/' + name + '_imputed_total_genes.h5ad')
        
        scanpy.pp.normalize_total(X, target_sum=1e4)
        scanpy.pp.log1p(X)
        expression = X.to_df().T
        # expression.index = X.var['gene_symbol']
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]
        # for j in cell_connected_gene[0]:
        #     if j not in expression.index:
        #         cell_connected_gene = cell_connected_gene[cell_connected_gene[0] != j]
     
        # i = 'data/Zhuang-ABCA-1110'
        l = []
        torch.cuda.empty_cache()
        all_pos = np.array([])
        all_neg = np.array([])
        matrix = expression[[str(j) for j in Merfish_groups[i]]]
        
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values

        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
            
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish, Merfish_groups)

        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            h = model(graph, graph.ndata['feat'])['cell']
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            #sigmoid of the score
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = compute_loss(pos_score, neg_score)
            optimizer.zero_grad()
            loss.backward()

            optimizer.step()

            l.append(loss.cpu().detach().numpy())
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)), file = f,flush = True)
        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)))

    for i in Starmap_g:
        l = []
        torch.cuda.empty_cache()
        all_pos = np.array([])
        all_neg = np.array([])
        expression = pd.read_csv(Starmap_dir + i + 'processed_total_expression_pd.csv', index_col = 0)
        expression = expression.loc[total_genes]
        matrix = expression[[str(j) for j in Starmap_groups[i]]]
        for j in cell_connected_gene[0]:
            if j not in expression.index:
                cell_connected_gene = cell_connected_gene[cell_connected_gene[0] != j]
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
            
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Starmap, edges_2_to_Starmap, edges_dict_Starmap, metadata_Starmap, Starmap_groups)


        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            h = model(graph, graph.ndata['feat'])['cell']
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = compute_loss(pos_score, neg_score)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            l.append(loss.cpu().detach().numpy())
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)), file = f,flush = True)
        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)))


    # backward
    model.eval()
    pred.eval()
    all_pos = np.array([])
    all_neg = np.array([])
    for i in Merfish_g_test:
        #load i.h5ad
        name = 'data_transmit/' + i.split('/')[-1]
        X = scanpy.read_h5ad('/home/bineva/Merfish/' + name + '_total_genes.h5ad')
        expression = X.to_df().T
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]
        for j in cell_connected_gene[0]:
            if j not in expression.index:
                cell_connected_gene = cell_connected_gene[cell_connected_gene[0] != j]


        torch.cuda.empty_cache()
        matrix = expression[[str(j) for j in Merfish_groups[i]]]
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        
        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish, Merfish_groups)
        h = model(graph, graph.ndata['feat'])['cell']
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

        X = scanpy.read_h5ad('/home/bineva/Merfish/' + name + '_imputed_total_genes.h5ad')
        #normalize and log1p
        scanpy.pp.normalize_total(X, target_sum=1e4)
        scanpy.pp.log1p(X)
        expression = X.to_df().T
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]

        for j in cell_connected_gene[0]:
            if j not in expression.index:
                cell_connected_gene = cell_connected_gene[cell_connected_gene[0] != j]


        torch.cuda.empty_cache()
        matrix = expression[[str(j) for j in Merfish_groups[i]]]
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        
        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish, Merfish_groups)
        h = model(graph, graph.ndata['feat'])['cell']
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))


    for i in Starmap_g_test:
        torch.cuda.empty_cache()
        expression = pd.read_csv(Starmap_dir + i + 'processed_total_expression_pd.csv', index_col = 0)
        expression = expression.loc[total_genes]
        matrix = expression[[str(j) for j in Starmap_groups[i]]]
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        
        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Starmap, edges_2_to_Starmap, edges_dict_Starmap, metadata_Starmap, Starmap_groups)
        h = model(graph, graph.ndata['feat'])['cell']
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

    for i in Slideseq_g_test:
        torch.cuda.empty_cache()
        X = scanpy.read_h5ad(Slideseq_dir + i[:i.find('dist10')] + '_total_genes.h5ad')
        scanpy.pp.normalize_total(X, target_sum=1e4)
        scanpy.pp.log1p(X)
        expression = X.to_df().T
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]
        matrix = expression[[str(j) for j in Slideseq_groups[i]]]
        nonzero = matrix.values.nonzero()
        if gene_2_cell:
            matrix.loc[[i for i in matrix.index if i not in cell_connected_gene[0].values]] = 0
        matrix = matrix.values
        index_gene1 = []
        index_gene2 = []
        for j in range(len(pairs)):
            if pairs['gene1'][j] in expression.index and pairs['gene2'][j] in expression.index:
                index_gene1.append(np.where(expression.index == pairs['gene1'][j])[0][0])
                index_gene2.append(np.where(expression.index == pairs['gene2'][j])[0][0])

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        
        if gene_2_gene:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int())),
                                    ('gene', 'gene2gene', 'gene'): ((torch.Tensor(index_gene1).to(device).int(),torch.Tensor(index_gene2).to(device).int()))}
        else:
            data_dict = {('gene', 'cell2gene', 'cell'): ((torch.Tensor(nonzero[0]).to(device).int(),torch.Tensor(nonzero[1]).to(device).int())),
                                    ('cell', 'gene2cell', 'gene'): ((torch.Tensor(nonzero[1]).to(device).int(),torch.Tensor(nonzero[0]).to(device).int()))}
        graph = dgl.heterograph(data_dict, num_nodes_dict={'gene': matrix.shape[0], 'cell': matrix.shape[1]}, device = device)

        graph.nodes['gene'].data['feat'] = torch.Tensor(list(range(matrix.shape[0]))).int().to(device)
        graph.nodes['cell'].data['feat'] = torch.Tensor(matrix.T).to(device)

        graph.edges['gene2cell'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        graph.edges['cell2gene'].data['w'] = torch.Tensor(matrix[nonzero]).to(device)
        if gene_2_gene:
            graph.edges['gene2gene'].data['w'] = torch.Tensor(np.ones(len(index_gene1))).to(device)
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(graph, i, 1024, edges_2_from_Slideseq, edges_2_to_Slideseq, edges_dict_Slideseq, metadata_Slideseq, Slideseq_groups)
        h = model(graph, graph.ndata['feat'])['cell']
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            h_p_u = h[p_u.long()].to(device)
            h_p_v = h[p_v.long()].to(device)
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = h[n_u.long()].to(device)
            h_n_v = h[n_v.long()].to(device)
            neg_score = pred(h_n_u, h_n_v)
            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

    print('test','AUC','AUPR', compute_auc(all_pos, all_neg), file = f,flush = True)
    print('test','AUC','AUPR', compute_auc(all_pos, all_neg))
    #draw the aupr curve
    if gene_2_gene:
        torch.save(model.state_dict(),'model_hetro_g_link_total_' + str(e) + '.pt')
        torch.save(pred.state_dict(),'pred_hetro_g_link_total_' + str(e) + '.pt')
    else:
        torch.save(model.state_dict(),'model_hetro_original_total_' + str(e) + '.pt')
        torch.save(pred.state_dict(),'pred_hetro_original_total_' + str(e) + '.pt')

