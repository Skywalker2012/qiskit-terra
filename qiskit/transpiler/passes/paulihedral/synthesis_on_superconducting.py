# ......
#
#
#




def opt_sc_backend(pauli_layers):
    circuit = 0
    return circuit


def connected_tree_synthesis1(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    pauli_map, graph, qc = synthesis_initial(pauli_layers, pauli_map, graph, qc, arch)
    for i1 in pauli_layers:
        for i2 in i1:
            lcover = compute_block_cover(i2)
            itir = compute_block_interior(i2)
            pcover = logical_list_physical(pauli_map, lcover)
            ptir = logical_list_physical(pauli_map, itir)
            lmc = -1
            lmi = -1
            lmt = []
            for i3 in ptir: # pauli_map, pcover
                dp = max_dfs_tree(graph, pcover, graph[i3])
                if len(dp) > lmc:
                    lmc = len(dp)
                    lmi = i3
                    lmt = dp
            lcover1 = physical_list_logical(graph, lmt)
            nc = []
            for i3 in lcover:
                if i3 not in lcover1:
                    nc.append(i3)
            # print(nc)
            ins = []
            while nc != []:
                id0, id1 = find_short_node(graph, pauli_map, nc, lmt)
                # print(id0, id1)
                connect_node(graph, pauli_map, pauli_map[id0], pauli_map[id1], ins)
                # do we need to update lmt?
                lmt.append(pauli_map[id0])
                nc.remove(id0)
            for i3 in ins:
                if i3[0] == 'swap':
                    qc.swap(i3[1][0], i3[1][1])
            pcover = logical_list_physical(pauli_map, lcover)
            # root is lmi
            dp = max_dfs_tree(graph, pcover, graph[lmi])
            dt = tree(graph, dp)
            for i3 in i2:
                tree_synthesis1(qc, graph, pauli_map, dt, i3)
    return qc