# ......
#
#
#


from qiskit.circuit.quantumcircuit import QuantumCircuit
import numpy as np


def matching_operator_positions(pauli_str_x, pauli_str_y):
    positions = []
    for i,  x_op in enumerate(pauli_str_x[0]):
        if x_op == pauli_str_y[0][i] and x_op != "I":
            positions.append(i)
    return positions


def cancellation_potential(pauli_layers):
    layer_num = len(pauli_layers)
    cancellation_potential_array = [0]*layer_num
    for i in range(1, layer_num):
        positions = matching_operator_positions(pauli_layers[i][0][0], pauli_layers[i-1][0][0])
        cancellation_potential_array[i - 1] = len(positions)
    return cancellation_potential_array


def finding_max_matching_layer_pairs(cancellation_potential_array):
    layer_num = len(cancellation_potential_array)
    unmatched_layer_list = list(range(layer_num))
    pairs = []
    while True:
        max_val = -1
        max_id = 0
        for i in range(len(unmatched_layer_list) - 1):
            if unmatched_layer_list[i+1] == unmatched_layer_list[i] + 1:
                if cancellation_potential_array[unmatched_layer_list[i]] > max_val:
                    max_id = i
                    max_val = cancellation_potential_array[unmatched_layer_list[i]]
        if max_val == -1:
            for i in unmatched_layer_list:
                pairs.append((i,i))
            break
        else:
            pairs.append((unmatched_layer_list[max_id], unmatched_layer_list[max_id]+1))
            unmatched_layer_list = unmatched_layer_list[:max_id]+unmatched_layer_list[max_id+2:]
    
    pairs = sorted(pairs, key=lambda pairs : pairs[0])
    return pairs

def non_identity_positions(pauli_str):
    positions = []
    for i, op in enumerate(pauli_str[0]):
        if op != 'I':
            positions.append(i)
    return positions


def cnot_tree(non_id_positions, outer_qubits):
    if len(non_id_positions) == 0:
        return [], -1
    elif len(non_id_positions) == 1:
        return [], non_id_positions[0]
    
    inner_qubits = []
    for position in non_id_positions:
        if position not in outer_qubits:
            inner_qubits.append(position)
    if outer_qubits != []:
        r1 = []
        r2 = []
        for qubit in inner_qubits:
            if qubit > outer_qubits[-1]:
                r2.append(qubit)
            else:
                r1.append(qubit)
    else:
        r1 = []
        r2 = inner_qubits
    outercnot_qubits = outer_qubits + r2
    left_cnotset = []
    for i in range(len(outercnot_qubits)-1):
        left_cnotset.append((outercnot_qubits[i],outercnot_qubits[i+1]))
    root = outercnot_qubits[-1]
    if r1 != []:
        for i in range(len(r1)-1):
            left_cnotset.append((r1[i],r1[i+1]))
        left_cnotset.append((r1[-1], root))
    right_cnotset = left_cnotset[::-1]
    return right_cnotset, root


def syn_pauli_string(qc, qubit_num, string, right_cnotset, root, coeff):

    for k in range(qubit_num):
        if string[k] == 'X':
            qc.h(k)
        if string[k] == 'Y':
            qc.s(k)
            qc.h(k)
    for k in reversed(right_cnotset):
        qc.cx(k[0], k[1])
        
    if root >= 0:
        qc.rz(2*coeff, root)
        
    for k in right_cnotset:
        qc.cx(k[0], k[1])
    for k in range(qubit_num):
        if string[k] == 'X':
            qc.h(k)
        if string[k] == 'Y':                    
            qc.h(k)
            qc.sdg(k)

def break_layer(pauli_layers):
    broken_pauli_layers = []
    for pauli_layer in pauli_layers:
        max_pauli_block_size = 2
        for pauli_block in pauli_layer:
            if len(pauli_block) > max_pauli_block_size:
                max_pauli_block_size = len(pauli_block)
                
        if max_pauli_block_size == 2:    
            broken_pauli_layers.append(pauli_layer)
        else:
            temp_layer = [0] * (max_pauli_block_size-1)
            for i in range(max_pauli_block_size-1):
                temp_layer[i] = []
            
            for pauli_block in pauli_layer:
                for i, pauli_str in enumerate(pauli_block[:-1]):
                    temp_layer[i].append([pauli_str, pauli_block[-1]])
            broken_pauli_layers += temp_layer
    return broken_pauli_layers        


def opt_ft_backend(pauli_layers):
    qubit_num = len(pauli_layers[0][0][0][0])
    
    pauli_layers = break_layer(pauli_layers)
    
    qc = QuantumCircuit(qubit_num)
    
    cancellation_potential_array = cancellation_potential(pauli_layers)
    
    pairs = finding_max_matching_layer_pairs(cancellation_potential_array)
    
    for pair in pairs:
        if pair[0] == pair[1]:
            pauli_block = pauli_layers[pair[0]][0]
            pauli_str = pauli_block[0]
            
            non_id_positions = non_identity_positions(pauli_str)
            right_cnotset, root = cnot_tree(non_id_positions, [])
            syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])
            
            for pauli_block in pauli_layers[pair[0]][1:]:
                pauli_str = pauli_block[0]
                non_id_positions = non_identity_positions(pauli_str)
                right_cnotset, root = cnot_tree(non_id_positions, [])
                syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])
        else:
            matching_positions = matching_operator_positions(pauli_layers[pair[0]][0][0], pauli_layers[pair[1]][0][0])
            non_id_positions = non_identity_positions(pauli_layers[pair[0]][0][0])
            right_cnotset, root = cnot_tree(non_id_positions, matching_positions)
            syn_pauli_string(qc, qubit_num, pauli_layers[pair[0]][0][0][0], right_cnotset, root, pauli_layers[pair[0]][0][0][1]*pauli_layers[pair[0]][0][-1])
            
            for pauli_block in pauli_layers[pair[0]][1:]:
                pauli_str = pauli_block[0]
                non_id_positions = non_identity_positions(pauli_str)
                right_cnotset, root = cnot_tree(non_id_positions, [])
                syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])
                    
            
            non_id_positions = non_identity_positions(pauli_layers[pair[1]][0][0])
            right_cnotset, root = cnot_tree(non_id_positions, matching_positions)
            syn_pauli_string(qc, qubit_num, pauli_layers[pair[1]][0][0][0], right_cnotset, root, pauli_layers[pair[1]][0][0][1]*pauli_layers[pair[1]][0][-1])
            
            for pauli_block in pauli_layers[pair[1]][1:]:
                pauli_str = pauli_block[0]
                non_id_positions = non_identity_positions(pauli_str)
                right_cnotset, root = cnot_tree(non_id_positions, [])
                syn_pauli_string(qc, qubit_num, pauli_str[0], right_cnotset, root, pauli_str[1]*pauli_block[-1])  

    return qc












class pNode:
    def __init__(self, idx):
        # self.child = []
        self.idx = idx
        self.adj = []
        self.lqb = None # logical qubit
        # self.parent = []
    def add_adjacent(self, idx):
        self.adj.append(idx)
    # def add_child(self, idx):
    #     if idx not in self.child:
    #         self.child.append(idx)
    # def add_parent(self, idx):
    #     if idx not in self.parent:
    #         self.parent.append(idx)

class pGraph:
    def __init__(self, G, C):
        n = G.shape[0]
        self.leng = n
        self.G = G # adj matrix
        self.C = C # cost matrix
        self.data = []
        self.coup_map = []
        for i in range(n):
            nd = pNode(i)
            for j in range(n):
                if G[i, j] == 1:
                    nd.add_adjacent(j)
                    self.coup_map.append([i, j])
            self.data.append(nd)
    def __getitem__(self, idx):
        return self.data[idx]
    def __len__(self):
        return self.leng

# sptSet: shortest path tree
def minDistance(dist, sptSet):
    # Initilaize minimum distance for next node
    minv = max_size
    min_index = -1
    n = len(dist)
    # Search not nearest vertex not in the
    for v in range(n):
        if dist[v] < minv and sptSet[v] == False:
            minv = dist[v]
            min_index = v
    return min_index
 
# Funtion that implements Dijkstra's single source
# shortest path algorithm for a graph represented
# using adjacency matrix representation
def dijkstra(dist_matrix, src):
    # print(src)
    n = dist_matrix.shape[0]
    dist = [dist_matrix[src, i] for i in range(n)]
    # print('Initial:', dist)
    sptSet = [False]*n # [dist_matrix[src, i] != max_size for i in range(n)]

    for cout in range(n):
        # Pick the minimum distance vertex from
        # the set of vertices not yet processed.
        # u is always equal to src in first iteration
        u = minDistance(dist, sptSet)
        if u == -1:
            break
        # Put the minimum distance vertex in the
        # shotest path tree
        sptSet[u] = True
        # Update dist value of the adjacent vertices
        # of the picked vertex only if the current
        # distance is greater than new distance and
        # the vertex in not in the shotest path tree
        for v in range(n):
            if dist_matrix[u, v] > 0 and sptSet[v] == False and \
                dist[v] > dist[u] + dist_matrix[u, v]:
                dist[v] = dist[u] + dist_matrix[u, v]
    # print(dist)
    for i in range(n):
        dist_matrix[src, i] = dist[i]
        dist_matrix[i, src] = dist[i]

def load_graph(code, dist_comp=False):
    pth = os.path.join(package_directory, 'data', 'ibmq_'+code+'_calibrations.csv')
    cgs = []
    n = 0
    with open(pth, 'r') as cf:
        g = csv.DictReader(cf, delimiter=',', quotechar='\"')
        for i in g:
            n += 1
            for j in i['CNOT error'].split(','):
                cgs.append(j.strip())
    G = np.zeros((n, n))
    C = np.ones((n, n))*max_size
    for i in range(n):
        C[i, i] = 0
    for i in cgs:
        si1 = i.find('_')
        si2 = i.find(':')
        iq1 = int(i[2:si1])
        iq2 = int(i[si1+1:si2])
        G[iq1, iq2] = 1
        acc = float(i[si2+2:])*1000
        C[iq1, iq2] = acc
    if dist_comp == True:
        for i in range(n):
            dijkstra(C, i)
    return G, C

def load_coupling_map(code):
    pth = os.path.join(package_directory, 'data', 'ibmq_'+code+'_calibrations.csv')
    cgs = []
    n = 0
    with open(pth, 'r') as cf:
        g = csv.DictReader(cf, delimiter=',', quotechar='\"')
        for i in g:
            n += 1
            for j in i['CNOT error'].split(','):
                cgs.append(j.strip())
    coupling = []
    for i in cgs:
        si1 = i.find('_')
        si2 = i.find(':')
        iq1 = int(i[2:si1])
        iq2 = int(i[si1+1:si2])
        coupling.append([iq1, iq2])
    return coupling




# def synthesis_initial(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
#     lnq = len(pauli_layers[0][0][0]) # logical qubits
#     if graph == None:
#         G, C = load_graph(arch, dist_comp=True) # G is adj, C is dist
#         graph = pGraph(G, C)
#     if pauli_map == None:
#         pauli_map = dummy_qubit_mapping(graph, lnq)
#     else:
#         add_pauli_map(graph, pauli_map)
#     pnq = len(graph) # physical qubits
#     if qc == None:
#         qc = QuantumCircuit(pnq)
#     return pauli_map, graph, qc


def dummy_qubit_mapping(graph, logical_qubit_num):
    for i in range(logical_qubit_num):
        graph[i].lqb = i
    return list(range(logical_qubit_num))


def compute_block_cover(pauli_block):
    block_cover = []
    for pauli_str in pauli_block[:-1]:
        for position in non_identity_positions(pauli_str):
            if position not in block_cover:
                block_cover.append(position)
    return block_cover

def compute_block_interior(pauli_block):
    block_interior = non_identity_positions(pauli_block[0])
    for pauli_str in pauli_block[1:-1]:
        block_interior = list(set(block_interior) & set(non_identity_positions(pauli_str)))
    return block_interior
                              
def logical_to_physical_cover(logical_to_physical_mapping, logical_cover):
    return [logical_to_physical_mapping[i] for i in logical_cover]

# def physical_list_logical(l):
    # return [i.lqb for i in l]
def physical_to_logical_cover(graph, physical_cover):
    return [graph[i].lqb for i in physical_cover]                              
                              

def max_dfs_tree(graph, cover, start, path=[]):
    max_connect_component = []
    for i in start.adj:
        if i not in cover:
            continue
        if i in path:
            continue
        else:
            max_connect_component += max_dfs_tree(graph, cover, graph[i], path=path + max_connect_component + [start.idx])
    return [start.idx] + max_connect_component



def find_short_node(graph, logical_to_physical_mapping, logical_qubits_to_move, max_connected_component):
    qubit_to_move = -1
    qubit_in_connected_component = -1
    mindist = 10000000
    for i in logical_qubits_to_move:
        for j in max_connected_component:
            dist = graph.C[logical_to_physical_mapping[i], j]
            if dist < mindist:
                mindist = dist
                qubit_to_move = i
                qubit_in_connected_component = graph[j].lqb
    return qubit_to_move, qubit_in_connected_component

def swap_nodes(logical_to_physical_mapping, graph_node_a, graph_node_b):
    t = graph_node_a.lqb
    graph_node_a.lqb = graph_node_b.lqb
    graph_node_b.lqb = t
    if graph_node_a.lqb != None:
        logical_to_physical_mapping[graph_node_a.lqb] = graph_node_a.idx
    if graph_node_b.lqb != None:
        logical_to_physical_mapping[graph_node_b.lqb] = graph_node_b.idx


def connect_node(graph, logical_to_physical_mapping, physical_qubit_to_move, physical_qubit_in_connected_component,  swap_list):
    physical_qubit_in_path = -1
    mindist = 10000000
    for i in graph[physical_qubit_to_move].adj:
        if graph.C[i, physical_qubit_in_connected_component] < mindist:
            physical_qubit_in_path = i
            mindist = graph.C[i, physical_qubit_in_connected_component]
    if physical_qubit_in_path == physical_qubit_in_connected_component:
        return
    else:
        swap_list.append(['swap',(physical_qubit_to_move, physical_qubit_in_path)])
        swap_nodes(logical_to_physical_mapping, graph[physical_qubit_to_move], graph[physical_qubit_in_path])
        connect_node(graph, logical_to_physical_mapping, physical_qubit_in_path, physical_qubit_in_connected_component,  swap_list)
        
        
class tree:
    def __init__(self, graph, dp, parent=None, depth=0):
        self.childs = []
        self.leaf = []
        self.depth = depth
        self.pid = dp[0]
        if len(dp) == 1:
            self.status = 1
            self.leaf = [self]
        else:
            self.status = 0   
            st = []
            for i in range(1, len(dp)):
                if dp[i] in graph[self.pid].adj:
                    st.append(i)
            st.append(len(dp))
            for i in range(len(st)-1):
                child = tree(graph, dp[st[i]:st[i+1]], parent=self, depth=self.depth+1)
                self.childs.append(child)
                self.leaf += child.leaf
        if parent != None:
            self.parent = parent
            
def pauli_single_gates(qc, logical_to_physical_mapping,  pauli_str, left=True):
    if left == True:
        for i in range(len( pauli_str[0])):
            if pauli_str[0][i] == 'X':
                qc.h(logical_to_physical_mapping[i])
            elif pauli_str[0][i] == 'Y':
                qc.s(logical_to_physical_mapping[i])
                qc.h(logical_to_physical_mapping[i])
    else:
        for i in range(len( pauli_str[0])):
            if pauli_str[0][i] == 'X':
                qc.h(logical_to_physical_mapping[i])
            elif pauli_str[0][i] == 'Y':
                qc.h(logical_to_physical_mapping[i])
                qc.sdg(logical_to_physical_mapping[i])
                

def tree_synthesis(qc, graph, logical_to_physical_mapping, physical_qubit_tree, pauli_str, block_coeff):
    string = pauli_str[0]
    non_id_positions = non_identity_positions(pauli_str)
    pauli_single_gates(qc, logical_to_physical_mapping, pauli_str, left=True)
    physical_qubit_tree_leaves = physical_qubit_tree.leaf # first in, first out
    swaps = {}
#    cnot_num = len(non_id_positions) - 1
#    lc = 0
    while physical_qubit_tree_leaves != []:
        physical_qubit_tree_leaves = sorted(physical_qubit_tree_leaves, key=lambda x: -x.depth)
        first_leaf = physical_qubit_tree_leaves[0]
        if first_leaf.depth == 0:
            # psd.real may be zero
            qc.rz(2 * pauli_str[1] * block_coeff, first_leaf.pid)
            break
        # actually, if psn is empty in the middle, we can stop it first
        # and the choice of root is also important
        if graph[first_leaf.pid].lqb in non_id_positions:
            if graph[first_leaf.parent.pid].lqb in non_id_positions:
                qc.cx(first_leaf.pid, first_leaf.parent.pid)
#                lc += 1
            else:
                qc.swap(first_leaf.pid, first_leaf.parent.pid)
                swaps[first_leaf.parent.pid] = first_leaf
                swap_nodes(logical_to_physical_mapping, graph[first_leaf.pid], graph[first_leaf.parent.pid])
        else:
            pass #lfs.remove(l)
        if first_leaf.parent not in physical_qubit_tree_leaves:
                physical_qubit_tree_leaves.append(first_leaf.parent)
        physical_qubit_tree_leaves.remove(first_leaf)
        # print(lfs)
#    if lc != cnot_num:
#        print('lala left:',psd.ps, cnot_num, lc)
    physical_qubit_tree_leaves = [physical_qubit_tree]
#    rc = 0
    while physical_qubit_tree_leaves != []:
        first_leaf = physical_qubit_tree_leaves[0]
        for i in first_leaf.childs:
            if graph[i.pid].lqb in non_id_positions:
                qc.cx(i.pid, first_leaf.pid)
#                rc += 1
                physical_qubit_tree_leaves.append(i)
        if first_leaf.pid in swaps.keys():
            qc.swap(first_leaf.pid, swaps[first_leaf.pid].pid)
            swap_nodes(logical_to_physical_mapping, graph[first_leaf.pid], graph[swaps[first_leaf.pid].pid])
            physical_qubit_tree_leaves.append(swaps[first_leaf.pid])
        physical_qubit_tree_leaves = physical_qubit_tree_leaves[1:]
#    if rc != cnot_num:
#        print('lala left:',psd.ps, cnot_num, rc)
    pauli_single_gates(qc, logical_to_physical_mapping, pauli_str, left=False)
    return qc

            
            
            
def opt_sc_backend(pauli_layers):
    
    logical_qubit_num = len(pauli_layers[0][0][0][0])
    
    # construct coupling graph (undirected), cost distance matrix
    adj_mat = 
    dist_mat = 
    sc_device_graph = pGraph(adj_mat, dist_mat)
    
    # find dense connected subset for initial mapping
    logical_to_physical_mapping = dummy_qubit_mapping(sc_device_graph, logical_qubit_num)
   
    physical_qubit_num = len(graph)
    
    qc = QuantumCircuit(physical_qubit_num)
    
    
    for current_layer in pauli_layers:
        for current_block in current_layer:
            logical_block_cover = compute_block_cover(current_block)
            logical_block_interior = compute_block_interior(current_block)
            physical_block_cover = logical_to_physical_cover(logical_to_physical_mapping, logical_block_cover)
            physical_block_interior = logical_to_physical_cover(logical_to_physical_mapping, logical_block_interior)
            
            max_connected_component_size = -1
            root_physical_qubit = -1
            max_connected_component = []
            
            for physical_qubit in physical_block_interior:
                connected_component = max_dfs_tree(sc_device_graph, physical_block_cover, sc_device_graph[physical_qubit])
                if len(connected_component) > max_connected_component_size:
                    max_connected_component_size = len(connected_component)
                    root_physical_qubit = physical_qubit
                    max_connected_component = connected_component
                    
            logical_connected_component_cover = physical_to_logical_cover(sc_device_graph, max_connected_component)
            
            logical_qubits_to_move = []
            for logical_qubit in logical_block_cover:
                if logical_qubit not in logical_connected_component_cover:
                    logical_qubits_to_move.append(logical_qubit)
                    
            swap_list = []
            while logical_qubits_to_move != []:
                qubit_to_move, qubit_in_connected_component = find_short_node(sc_device_graph, logical_to_physical_mapping, logical_qubits_to_move, max_connected_component)
                connect_node(sc_device_graph, logical_to_physical_mapping, logical_to_physical_mapping[qubit_to_move], logical_to_physical_mapping[qubit_in_connected_component], swap_list)
                max_connected_component.append(logical_to_physical_mapping[qubit_to_move])
                logical_qubits_to_move.remove(qubit_to_move)
                
            for swap_step in swap_list:
                qc.swap(swap_step[1][0], swap_step[1][1])
            
            physical_block_cover = logical_to_physical_cover(logical_to_physical_mapping, logical_block_cover)  
            
            connected_component = max_dfs_tree(sc_device_graph, physical_block_cover, sc_device_graph[root_physical_qubit])
            physical_qubit_tree = tree(sc_device_graph, connected_component)
            for pauli_str in current_block[:-1]:
                tree_synthesis1(qc, sc_device_graph, logical_to_physical_mapping, physical_qubit_tree, pauli_str, current_block[-1])
    return qc

# '''-------------'''
# def connected_tree_synthesis1(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
#     pauli_map, graph, qc = synthesis_initial(pauli_layers, pauli_map, graph, qc, arch)
#     for i1 in pauli_layers:
#         for i2 in i1:
#             lcover = compute_block_cover(i2)
#             itir = compute_block_interior(i2)
#             pcover = logical_list_physical(pauli_map, lcover)
#             ptir = logical_list_physical(pauli_map, itir)
#             lmc = -1
#             lmi = -1
#             lmt = []
#             for i3 in ptir: # pauli_map, pcover
#                 dp = max_dfs_tree(graph, pcover, graph[i3])
#                 if len(dp) > lmc:
#                     lmc = len(dp)
#                     lmi = i3
#                     lmt = dp
#             lcover1 = physical_list_logical(graph, lmt)
#             nc = []
#             for i3 in lcover:
#                 if i3 not in lcover1:
#                     nc.append(i3)
#             # print(nc)
#             ins = []
#             while nc != []:
#                 id0, id1 = find_short_node(graph, pauli_map, nc, lmt)
#                 # print(id0, id1)
#                 connect_node(graph, pauli_map, pauli_map[id0], pauli_map[id1], ins)
#                 # do we need to update lmt?
#                 lmt.append(pauli_map[id0])
#                 nc.remove(id0)
#             for i3 in ins:
#                 if i3[0] == 'swap':
#                     qc.swap(i3[1][0], i3[1][1])
#             pcover = logical_list_physical(pauli_map, lcover)
#             # root is lmi
#             dp = max_dfs_tree(graph, pcover, graph[lmi])
#             dt = tree(graph, dp)
#             for i3 in i2:
#                 tree_synthesis1(qc, graph, pauli_map, dt, i3)
#     return qc