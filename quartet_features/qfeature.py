from tqdist import *
from triproot import *
from treeswift import *
from math import sqrt
from utils import reroot_at_edge

def __reroot_and_prune__(nwk_tree,query_taxon,active_leafset=None):                       
# reroot the tree at the query_taxon and restrict it to the active_leafset
    tree_obj = read_tree_newick(nwk_tree)
    new_root = None
    for leaf in tree_obj.traverse_leaves():
        if leaf.label == query_taxon:
            new_root = leaf
            break     
    
    reroot_at_edge(tree_obj,new_root,new_root.edge_length/2 if new_root.edge_length is not None else None)  
    C = tree_obj.root.children
    c_star = None
    
    for c in C:
        if c.label == query_taxon:
            c_star = c
            break
    
    tree_obj.root.remove_child(c_star) 
    tree_obj.suppress_unifurcations()
    
    if active_leafset is not None:
        tree_obj = tree_obj.extract_tree_with(active_leafset,suppress_unifurcations=True)
    n = len([leaf for leaf in tree_obj.traverse_leaves()])
    return tree_obj.newick().replace("[&R] ",""),n

def get_qfeatures(backboneTree,refTrees,query_taxon):
# assumptions: + all nodes in backboneTree have unique labels
#              + query_taxon presents in each of the refTrees but not in backboneTree
    
    tree_obj = read_tree_newick(backboneTree)
    active_leafset = set(leaf.label for leaf in tree_obj.traverse_leaves())
    id2qScores = {}

    for refTree in refTrees:
        refTree_pruned, n = __reroot_and_prune__(refTree,query_taxon,active_leafset)         
        nquarts_b = n*(n-1)*(n-2)*(n-3)/24 # number of quartets before placement
        nquarts_a = (n+1)*n*(n-1)*(n-2)/24 # number of quartets after placement
        ntrips = n*(n-1)*(n-2)/6
        qscore = (1-quartet_distance(refTree_pruned,backboneTree))*nquarts_b/nquarts_a
        mystr = tripRootScore(refTree_pruned,backboneTree)
        if not mystr: # empty string returned by tripRootScore; usually because the refTree is too small after pruning to pair with the query tree
            continue
        for item in mystr.split(','):
            ID,s = item.split(':')
            score = float(s)*ntrips/nquarts_a + qscore
            if ID not in id2qScores:
                id2qScores[ID] = [score]
            else:
                id2qScores[ID].append(score)
                    
    return id2qScores

def get_qfeatures_training(backboneTree,branches,refTrees,background_qscores):
    tree_obj = read_tree_newick(backboneTree)
    k = len(refTrees)
    
    X = []
    for leaf in tree_obj.traverse_leaves():
        backbone_pruned = tree_obj.extract_tree_without([leaf.label],suppress_unifurcations=False).newick().replace("[&R] ","")
        id2qScores = get_qfeatures(backbone_pruned,refTrees,leaf.label)
        id2qScores[leaf.label] = background_qscores
        for b in branches:
            X.append((leaf.label,b,id2qScores[b]))
    return X        

def create_feature_set(backboneTree,refTrees,query_taxon):
    backboneTree_pruned, placement_pos, branches = prepare_backbone(backboneTree,query_taxon)
    background_qscores = compute_background_qscores(backboneTree_pruned,refTrees)

    X_train = get_qfeatures_training(backboneTree_pruned,branches,refTrees,background_qscores)           
    id2qScores = get_qfeatures(backboneTree_pruned,refTrees,query_taxon)
    x_test = []
    for ID in branches:
        x_test.append((query_taxon,ID,id2qScores[ID]))

    return X_train, x_test, backboneTree_pruned, placement_pos

def prepare_backbone(backboneTree,query_taxon):
    backboneTree_pruned,_ = __reroot_and_prune__(backboneTree,query_taxon, active_leafset=None)
    tree_obj = read_tree_newick(backboneTree_pruned)
    C = [node.label for node in tree_obj.root.children]
    tree_obj.deroot()
    placement_pos = None
    for node in tree_obj.root.children:
        if node.label in C:
            placement_pos = node.label
            break
    branches = [node.label for node in tree_obj.traverse_preorder() if not node.is_root()]
            
    return tree_obj.newick(),placement_pos,branches

def compute_background_qscores(backboneTree_pruned,refTrees):
    tree_obj = read_tree_newick(backboneTree_pruned)
    active_leafset = set(leaf.label for leaf in tree_obj.traverse_leaves())
    qscores = []
    for ref in refTrees:
        tree_obj = read_tree_newick(ref)
        tree_obj = tree_obj.extract_tree_with(active_leafset,suppress_unifurcations=True)
        qscores.append(1 - quartet_distance(backboneTree_pruned,tree_obj.newick().replace("[&R] ","")))
    return qscores

if __name__ == "__main__":
    from sys import argv
    from os import mkdir

    backbone_file = argv[1]
    reftrees_file = argv[2]
    query_taxon = argv[3]
    outdir = argv[4]
    
    with open(backbone_file,'r') as f:
        backboneTree = f.read().strip()
    with open(reftrees_file,'r') as f:
        refTrees = f.read().strip().split("\n")    
 
    X_train,x_test,backboneTree_pruned,placement_pos = create_feature_set(backboneTree,refTrees,query_taxon)

    mkdir(outdir)
    with open(outdir+"/X_train.txt",'w') as fout:
        for sp,br,qscores in X_train:
            fout.write(str(sp) + " " + (br) + " " + str(qscores) + "\n")
    with open(outdir+"/x_test.txt",'w') as fout:
        for sp,br,qscores in x_test:
            fout.write(str(sp) + " " + (br) + " " + str(qscores) + "\n")        
    with open(outdir+"/info.txt",'w') as fout:
        fout.write("Backbone tree: " + backboneTree_pruned + "\n")
        fout.write("True placement: " + placement_pos)        
