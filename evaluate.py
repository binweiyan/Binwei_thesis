i = 'well01brain'
import pickle
gradients = {}
with open('gradients/gradients_ordinary_' + i + '.pkl', 'rb') as f:
    gradients[i] = pickle.load(f)
from scipy.stats import fisher_exact
import numpy as np
import matplotlib.pyplot as plt
ratio = []
overlap = []
total = []
fisher_exact_test = {}
for i in gradients:
    for j in gradients[i]:
        # print(j)
        if gradients1[i][j].shape[0] < 20: continue
        # if gradients1[i][j].shape != gradients[i][j].shape: continue
        
        cut = np.quantile( gradients[i][j], 0.80)
        new = gradients[i][j] > cut
        new.sum(axis = 0)
        cut_gene = np.quantile(new.sum(axis = 0), 0.80)
        good_gene = new.sum(axis = 0) > cut_gene
        # print(len(Starmap_index[good_gene[:1022]]) / len(Starmap_index[good_gene[1022:]]))
        # how many of Starmap_index[good_gene[:1022]] are in cell_connected_gene[0].values
        if len(Starmap_index[good_gene[:1022]]) == 0 or len(Starmap_index[good_gene[1022:]]) == 0:
            continue
        ratio.append(len([i for i in Starmap_index[good_gene[:1022]] if i in cell_connected_gene[0].values]) / len(Starmap_index[good_gene[:1022]]))
        ratio.append(len([i for i in Starmap_index[good_gene[1022:]] if i in cell_connected_gene[0].values]) / len(Starmap_index[good_gene[1022:]]))
        overlap.append(len([i for i in Starmap_index[good_gene[:1022]] if i in cell_connected_gene[0].values]) + len([i for i in Starmap_index[good_gene[1022:]] if i in cell_connected_gene[0].values]))
        total.append(len(Starmap_index[good_gene[:1022]]) + len(Starmap_index[good_gene[1022:]]))
        table = [[overlap[-1], len(cell_connected_gene[0].values) * 2 - overlap[-1] ], [total[-1] - overlap[-1], len(Starmap_index) * 2 - len(cell_connected_gene[0].values) * 2 - total[-1] + overlap[-1]]]
        # print(table)
        fisher_exact_test[j] = fisher_exact(table)[1]


# print(np.mean([i[1] for i in fisher_exact_test]))


