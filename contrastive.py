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

    return df, edges_2_from, edges_2_to, edges_dict, g, g_test

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
    return df, edges_2_from, edges_2_to, edges_dict, g, g_test

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
# In[3]:


metadata_Merfish = pd.read_pickle(Merfish_dir + 'metadata.pkl')
metadata_Starmap = pd.read_pickle(Starmap_dir + 'metadata.pkl')
metadata_Slideseq = pd.read_pickle(Slideseq_dir + 'metadata.pkl')

Merfish_df, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, Merfish_g, Merfish_g_test = process_Merfish_df()
Starmap_df, edges_2_from_Starmap, edges_2_to_Starmap, edges_dict_Starmap, Starmap_g, Starmap_g_test = process_Starmap_df()
Slideseq_df, edges_2_from_Slideseq, edges_2_to_Slideseq, edges_dict_Slideseq, Slideseq_g, Slideseq_g_test, Slideseq_groups = process_Slideseq_df()

# In[4]:



# ----------- 2. create model -------------- #
# build a two-layer GraphSAGE model
import torch
import torch.nn as nn
import torch.nn.functional as F

import numpy as np
import scipy.sparse as sp

# In[6]:


class contrastive(nn.Module):
    def __init__(self, in_features, out_features, hidden_dim = 256):
        super(contrastive, self).__init__()
        self.fc1 = nn.Linear(in_features, hidden_dim)
        # self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.fc3 = nn.Linear(hidden_dim, out_features)
        self.dropout = nn.Dropout(0.8)

        self.batchnorm1 = nn.BatchNorm1d(64)
        self.batchnorm2 = nn.BatchNorm1d(64)
        
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        # x = self.batchnorm1(x)
        # x = self.dropout(x)
        # x = F.relu(self.fc2(x))
        # x = self.batchnorm2(x)
        # x = self.dropout(x)
        x = self.fc3(x)
        return x
    
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



# In[8]:


with open('/home/bineva/index/total_genes.pkl', 'rb') as f:
    total_genes = pickle.load(f)
feature_size = len(total_genes)
out_features = 256
model = contrastive(feature_size, out_features).to(device)

# model.load_state_dict(torch.load('model_hetro2.pt'))
pred = MLPPredictor(256).to(device)


# In[10]:


with open('/home/bineva/index/bad.pkl', 'rb') as f:
    bad = pickle.load(f)

# In[30]:

import scanpy
import itertools
expression_new = {}
import warnings
os.chdir('/home/bineva/train_all')
file = 'log_contrastive_total.txt'
resume_True = True
epoch = 6
model_path = ''
if  resume_True:
    if gene_2_gene:
        model.load_state_dict(torch.load(model_path + 'model_contrastive_' + str(epoch) + '.pt'))
        pred.load_state_dict(torch.load(model_path + 'pred_contrastive_' + str(epoch) + '.pt'))
else:
    epoch = -1
if os.path.exists(file):
    f = open(file, 'a')
else:
    f = open(file, 'w')


warnings.filterwarnings("ignore")

optimizer = torch.optim.Adam(itertools.chain(model.parameters(), pred.parameters()), lr=1e-4)

#scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.7, patience=5, min_lr=1e-5)


def batch_negative_edge( i, batchsize, edges_2_from, edges_2_to, edges_dict, metadata):
    u, v = torch.tensor(edges_2_from[i] + edges_2_to[i]).int() , torch.tensor(edges_2_to[i] + edges_2_from[i]).int()

    eids = np.arange(len(edges_2_from[i] + edges_2_to[i]))
    eids = np.random.permutation(eids)
    pos_u, pos_v = u[eids], v[eids]

    # Find all negative edges and split them for training and testing
    # if not too many nodes, you can use the following code to generate all negative edges.
    if len(metadata[metadata['label'] == i]) < 100:
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
            neg_u = np.concatenate((neg_u, np.random.choice(len(metadata[metadata['label'] == i]), len(edges_2_from[i] + edges_2_to[i]) - len(neg_u))))
            neg_v = np.concatenate((neg_v, np.random.choice(len(metadata[metadata['label'] == i]), len(edges_2_from[i] + edges_2_to[i]) - len(neg_v))))
            # remove existing edges, use edges_dict to check
            for j in range(len(neg_u)):
                while neg_v[j] in edges_dict[i][neg_u[j]] or neg_u[j] == neg_v[j]:
                    neg_u[j] = np.random.choice(len(metadata[metadata['label'] == i]))
                    neg_v[j] = np.random.choice(len(metadata[metadata['label'] == i]))

    #divide the train and test sets into batches, and train the model on each batch.
    #the positive batch:
    if batchsize > len(pos_u):
        batchsize = len(pos_u)
    pos_size = batchsize // 2
    neg_size = batchsize - pos_size
    pos_split_u, pos_split_v = torch.split(pos_u, pos_size), torch.split(pos_v, pos_size)
    #turn neg_u and neg_v into tensors
    neg_u, neg_v = torch.tensor(neg_u), torch.tensor(neg_v)
    neg_split_u, neg_split_v = torch.split(neg_u, neg_size), torch.split(neg_v, neg_size)

    return  pos_split_u, pos_split_v, neg_split_u, neg_split_v



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
        expression = expression[[str(j) for j in Slideseq_groups[i]]]
        matrix = expression

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()

        matrix = matrix.values


        pos_u, pos_v, neg_u, neg_v = batch_negative_edge( i, 1024, edges_2_from_Slideseq, edges_2_to_Slideseq, edges_dict_Slideseq, metadata_Slideseq)
        
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            
            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
            neg_score = pred(h_n_u, h_n_v)

            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = torch.mean(neg_score - pos_score)
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
        matrix = expression
        
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()

        matrix = matrix.values
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish)


        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            
            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
            neg_score = pred(h_n_u, h_n_v)

            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = torch.mean(neg_score - pos_score)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            l.append(loss.cpu().detach().numpy())
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)), file = f,flush = True)
        print(e, i,'AUC','AUPR', compute_auc(all_pos, all_neg),'loss: ' + str(np.mean(l)))


        X = scanpy.read_h5ad('/home/bineva/Merfish/' + name + '_imputed_total_genes.h5ad')
        #normalize and log1p
        scanpy.pp.normalize_total(X, target_sum=1e4)
        scanpy.pp.log1p(X)
        expression = X.to_df().T
        expression.index = expression.index.str.upper()
        expression = expression.loc[total_genes]
        # i = 'data/Zhuang-ABCA-1110'
        l = []
        torch.cuda.empty_cache()
        all_pos = np.array([])
        all_neg = np.array([])
        matrix = expression
        
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()

        matrix = matrix.values
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge(i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish)


        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            
            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
            neg_score = pred(h_n_u, h_n_v)

            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = torch.mean(neg_score - pos_score)
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
        matrix = expression
        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
        nonzero = matrix.values.nonzero()
        matrix = matrix.values
        
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge( i, 1024, edges_2_from_Starmap, edges_2_to_Starmap, edges_dict_Starmap, metadata_Starmap)

        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):

            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
            neg_score = pred(h_n_u, h_n_v)

            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            loss = torch.mean(neg_score - pos_score)
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
        torch.cuda.empty_cache()

        matrix = expression
        nonzero = matrix.values.nonzero()

        matrix = matrix.values

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
    
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge( i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish)

        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):

            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
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

        torch.cuda.empty_cache()

        matrix = expression
        nonzero = matrix.values.nonzero()

        matrix = matrix.values

        #sigmoid of the matrix
        # matrix = 1 / (1 + np.exp(-matrix))
    
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge( i, 1024, edges_2_from_Merfish, edges_2_to_Merfish, edges_dict_Merfish, metadata_Merfish)

        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):

            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
            neg_score = pred(h_n_u, h_n_v)

            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)

            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))


    for i in Starmap_g_test:
        torch.cuda.empty_cache()
        expression = pd.read_csv(Starmap_dir + i + 'processed_total_expression_pd.csv', index_col = 0)
        expression = expression.loc[total_genes]
        matrix = expression
        nonzero = matrix.values.nonzero()

        matrix = matrix.values

        pos_u, pos_v, neg_u, neg_v = batch_negative_edge( i, 1024, edges_2_from_Starmap, edges_2_to_Starmap, edges_dict_Starmap, metadata_Starmap)

        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):

            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
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

        expression = expression[[str(j) for j in Slideseq_groups[i]]]
        matrix = expression
        nonzero = matrix.values.nonzero()

        matrix = matrix.values
        pos_u, pos_v, neg_u, neg_v = batch_negative_edge( i, 1024, edges_2_from_Slideseq, edges_2_to_Slideseq, edges_dict_Slideseq, metadata_Slideseq)
        
        for p_u, p_v, n_u, n_v in zip( pos_u, pos_v, neg_u, neg_v):
            
            p_u = torch.Tensor(expression.values[:, p_u].T).to(device)
            p_v = torch.Tensor(expression.values[:, p_v].T).to(device)
            n_u = torch.Tensor(expression.values[:, n_u.int()].T).to(device)
            n_v = torch.Tensor(expression.values[:, n_v.int()].T).to(device)

            h_p_u = model(p_u.float().to(device))
            h_p_v = model(p_v.float().to(device))
            pos_score = pred(h_p_u, h_p_v)
            h_n_u = model(n_u.float().to(device))
            h_n_v = model(n_v.float().to(device))
            neg_score = pred(h_n_u, h_n_v)

            pos_score, neg_score = torch.sigmoid(pos_score), torch.sigmoid(neg_score)
            all_pos = np.concatenate((all_pos, pos_score.cpu().detach().numpy()))
            all_neg = np.concatenate((all_neg, neg_score.cpu().detach().numpy()))

    print('test','AUC','AUPR', compute_auc(all_pos, all_neg), file = f,flush = True)
    print('test','AUC','AUPR', compute_auc(all_pos, all_neg))
    #draw the aupr curve

    torch.save(model.state_dict(),'model_contrastive_' + str(e) + '.pt')
    torch.save(pred.state_dict(),'pred_contrastive_' + str(e) + '.pt')