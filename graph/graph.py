#! /usr/bin/env python

from __future__ import division, absolute_import, print_function
import numpy as np
import math
import random
import networkx as nx
import matplotlib.pyplot as plt
import json

def Read_Reactions():
    mollist = []
    #etot = []

    # read molecules and energies
    with open('mols.txt','r') as f:
        for line in f:    
            dat = line.split()
            mollist.append(dat[1])
            #etot.append(float(dat[1]))

    reactions = []

    # read reactions 
    with open('react.txt','r') as f:
        for line in f:          
            dat = line.split()
            reactions.append([int(dat[0]),dat[1],dat[3],dat[5]])

    return mollist,reactions


print("reading reaction network")
mollist,reactions = Read_Reactions()

#print(mollist)
#print(reactions)

G = nx.DiGraph()
labels = {}
mols = []
re = []
edges = []

for m in mollist:
    G.add_node(m)
    labels[m] = m
    mols.append(m)

for r in reactions:
#    G.add_node(r[0],size=1)
#    G.add_edge(r[1],r[0],{'myweight':0.0})
    G.add_edge(r[1],r[2])
    G.add_edge(r[1],r[3])
#    labels[r[0]] =  1
#    labels[r[0]] = str(r[0])
    re.append(r[0])
#    edges.append((r[1],r[0]))
#    edges.append((r[0],r[2]))
#    edges.append((r[0],r[3]))

#for im1 in range(len(mollist)-2):
#    G.add_edge(mollist[im1],mollist[im1+1],{'myweight':100.0})        
#    G.add_edge(mollist[im1],mollist[im1+2],{'myweight':100.0})


for ix,deg in G.degree().items():
    G.node[ix]['degree'] = deg
    G.node[ix]['parity'] = (1-deg%2)

for ix,katz in nx.katz_centrality(G).items():
    G.node[ix]['katz'] = katz

from networkx.readwrite import json_graph

data = json_graph.node_link_data(G)

with open('graph.json', 'w') as f:
    json.dump(data, f, indent=4)


#print(G.nodes(data=True)[:5])

#pos = nx.spring_layout(G,weight='myweight',k=0.5)
#pos = nx.spectral_layout(G)
#pos = nx.random_layout(G)
#pos = nx.spring_layout(G,weight='myweight',k=10)
#pos = nx.spring_layout(G,k=1e-3,fixed=mols,pos=pos)

#nx.draw_networkx_nodes(G,pos,nodelist=mols,node_size=1300)
#nx.draw_networkx_labels(G,pos,labels=labels)

#pos = nx.spring_layout(G,weight='myweight',k=1e-2,fixed=mols,pos=pos)

#nx.draw_networkx_nodes(G,pos,nodelist=re,node_size=1,alpha=0.5)
#nx.draw_networkx_edges(G,pos,arrows=False,edgelist=edges)

#nx.draw_spring(G, labels=labels, with_labels=True,arrows=False)
#nx.draw_spring(G,arrows=False)
#nx.draw_random(g, node_size = node_sizes, labels=labels, with_labels=True)    

#plt.savefig("pic.png")



