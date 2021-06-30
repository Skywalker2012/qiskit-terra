
#from qiskit.transpiler.exceptions import 


#from qiskit.circuit.quantumcircuit import QuantumCircuit
#from qiskit.circuit.quantumregister import Qubit


def block_len(pauli_block):
    '''example block: [XYX, p1], [XZZ, p2], ..., [ZZZ, p5], p]'''
    qubit_num = len(pauli_block[0][0])
    length = 0
    for i in range(qubit_num):
        for pauli_str in pauli_block[0:-1]:
            if pauli_str[0][i] != 'I':
                length += 1
                break
    return length
    
def str_lex_key(weighted_pauli_str):
    value = 0
    '''q0 corresponds to the right-most pauli op'''
    for op in weighted_pauli_str[0]:
        value *= 4
        if op == 'I':
            value += 0
        elif op == 'X':
            value += 1
        elif op == 'Y':
            value += 2
        elif op == 'Z':
            value += 3
    return value

def block_lex_key(pauli_block):
    value = 0
    '''q0 corresponds to the right-most pauli op'''
    for op in pauli_block[0][0]:
        value *= 4
        if op == 'I':
            value += 0
        elif op == 'X':
            value += 1
        elif op == 'Y':
            value += 2
        elif op == 'Z':
            value += 3
    return value

def block_latency(pauli_block, qubit_num):
    latency = 0
    for pauli_str in pauli_block[:-1]:
        latency += max(2*(qubit_num - pauli_str[0].count('I'))-1, 0)
        if qubit_num - pauli_str[0].count('I') - pauli_str[0].count('Z') > 0:
            latency += 1 # or 2
            
    print(latency)
    return latency

def block_qubit_occupation(pauli_block, qubit_num):
    qubit_occupied = [0] * qubit_num
    for i in range(qubit_num):
        for pauli_str in pauli_block[0:-1]:
            if pauli_str[0][i] != 'I':
                qubit_occupied[i] = 1
            break
    print(qubit_occupied)
    return qubit_occupied
    

def qubit_occupation_template(qubit_occupied, qubit_num, latency):
    template = []
    start_index = 0
    end_index = 0
    
    '''0 -> qubit not occupied (I), 1 -> qubit occupied (X,Y,Z)'''
    while start_index < qubit_num:
        while end_index < qubit_num and qubit_occupied[end_index] == 0:
            end_index += 1
        if(end_index-start_index) >= 2:
            window = [1]*start_index + [0]*(end_index-start_index) + [1]*(qubit_num-end_index)
            template.append([window, latency])
        while end_index < qubit_num and qubit_occupied[end_index] != 0:
            end_index += 1
        start_index = end_index
        
    print(template)
    return template

def occupation_embedding_check(qubit_occupied_x, qubit_occupied_y):
    for i, occupied_x in enumerate(qubit_occupied_x):
        if occupied_x == 1 and qubit_occupied_y[i] == 1:
            return False
    return True


def qubit_occupation_combine(qubit_occupied_x, qubit_occupied_y):
    combined_occupation = [0]*len(qubit_occupied_x)
    for i, occupied_x in enumerate(qubit_occupied_x):
        if occupied_x == 1 or qubit_occupied_y[i] == 1:
            combined_occupation[i] = 1
    return combined_occupation

            
            
# def generate_templates1(ps, nq, latency, lt=2):
#     l = latency
#     idx0 = 0
#     idx1 = 0
#     t = []
#     while idx0 < nq:
#         while idx1 < len(ps) and ps[idx1] == "I":
#             idx1 += 1
#         if (idx1 - idx0) >= 2:
#             ts = (idx0-0)*'I'+(idx1-idx0)*'X'+(nq-idx1)*'I'
#             t.append([ts, l])
#         while idx1 < len(ps) and ps[idx1] != "I":
#             idx1 += 1
#         idx0 = idx1
#     return t

        
def lex_ordering(pauliIRprogram):
    for i, block in enumerate(pauliIRprogram):
        param = block[-1]
        pauliIRprogram[i] = sorted(block[0:-1], key=str_lex_key)
        pauliIRprogram[i].append(param)
    pauliIRprogram = sorted(pauliIRprogram, key=block_lex_key)
    return pauliIRprogram
    
    
def gco_ordering(pauliIRprogram, **kwargs):
    pauliIRprogram = lex_ordering(pauliIRprogram)
    pauli_layers = [0] * len(pauliIRprogram)
    for i in range(len(pauliIRprogram)):
        pauli_layers[i] = [pauliIRprogram[i]]
        
    print(pauli_layers)
    return pauli_layers
    
def do_ordering(pauliIRprogram, max_iteration = 20):
    qubit_num = len(pauliIRprogram[0][0][0])
    
    if qubit_num >= 4:
        large_block_threshold = qubit_num/2 - 2
    else:
        large_block_threshold = 2
    
    large_block_list = []
    small_block_list = []
    
    for block in pauliIRprogram:
        if block_len(block) > large_block_threshold:
            large_block_list.append(block)
        else:
            small_block_list.append(block)
            
    large_block_list = lex_ordering(large_block_list)
    small_block_list = lex_ordering(small_block_list)
    
    print("large_block: ", large_block_list)
    print("small_block: ", small_block_list)
    
    block_num = len(pauliIRprogram)
    
    pauli_layers = []
    
    while block_num > 0:
        if len(large_block_list) > 0:
            current_block = large_block_list[0]
            large_block_list = large_block_list[1:]
        elif len(small_block_list) > 0:
            current_block = small_block_list[0]
            small_block_list = small_block_list[1:]
            
        current_layer = [current_block]
        
        block_num -= 1
        
        latency = block_latency(current_block, qubit_num)
        layer_qubit_occupation = block_qubit_occupation(current_block, qubit_num)
        current_template = qubit_occupation_template(layer_qubit_occupation, qubit_num, latency)
        
        iteration_count = 0
        while len(current_template) > 0 and iteration_count < max_iteration:
            iteration_count += 1
            qubit_occupation_list = []
            for template in current_template:
                small_block_not_selected_list = []
                for small_block in small_block_list:
                    small_latency = block_latency(small_block, qubit_num)
                    small_occupation = block_qubit_occupation(small_block, qubit_num)
                    if small_latency <= template[1] and occupation_embedding_check(template[0], small_occupation):
                        current_layer.append(small_block)
                        qubit_occupation_list.append(small_occupation)
                        template[1] -= small_latency
                        block_num -= 1
                    else:
                        small_block_not_selected_list.append(small_block)
                small_block_list = small_block_not_selected_list
            for qubit_occupation in qubit_occupation_list:
                layer_qubit_occupation = qubit_occupation_combine(qubit_occupation, layer_qubit_occupation)
            current_template = qubit_occupation_template(layer_qubit_occupation, qubit_num, latency)
            if len(qubit_occupation_list) == 0:
                break
        pauli_layers.append(current_layer)
                        
                        
    print(pauli_layers)
    
    return pauli_layers

# '''--------'''
# def parallel_order_size_bl1(parr, length=0, maxiter=1):
#     nq = len(parr[0][0]) # num_qubits
#     parr_flat = []
#     for i in parr:
#         parr_flat.append(pauli_block(i))    
#     parr_big = []
#     parr_small = []
#     for i in parr_flat:
#         if i.len > length:
#             parr_big.append(i)
#         else:
#             parr_small.append(i)
            
#     def _key(pb):
#         s = 0
#         ps = pb.mstr
#         for i in ps:
#             s *= 4
#             if i == 'I':
#                 s += 0
#             elif i == 'X':
#                 s += 1
#             elif i == 'Y':
#                 s += 2
#             elif i == 'Z':
#                 s += 3
#         return -s
#     parr_big = sorted(parr_big, key=_key)
#     parr_small = sorted(parr_small, key=_key)
#     # print(parr_small)
#     ne = len(parr_flat)
#     ps_layers = []
#     while ne > 0:
#         if len(parr_big) > 0:
#             cpb = parr_big[0]
#             parr_big = parr_big[1:]
#         elif len(parr_small) > 0:
#             cpb = parr_small[0]
#             parr_small = parr_small[1:]
#         pl = [cpb.block]
#         ne -= 1
#         pso = cpb.mstr
#         latency = cpb.latency
#         pts = generate_templates1(pso, nq, latency)
#         cnt = 0
#         while len(pts) > 0 and cnt < maxiter:
#             cnt += 1
#             tmp1 = []
#             for i in pts:
#                 tmp2 = []
#                 for j in range(len(parr_small)):
#                     if parr_small[j].latency <= i[1] and pINC(i[0], parr_small[j].mstr):
#                         pl.append(parr_small[j].block)
#                         tmp1.append(parr_small[j].mstr)
#                         i[1] -= parr_small[j].latency
#                         ne -= 1
#                     else:
#                         tmp2.append(parr_small[j])
#                 parr_small = tmp2
#             for i in tmp1:
#                 pso = pOR(pso, i)
#             pts = generate_templates1(pso, nq, latency)
#             if len(tmp1) == 0:
#                 break
#         ps_layers.append(pl)
#     return ps_layers
# '''--------'''


def opt_ft_backend(pauli_layers):
    circuit = 0
    return circuit

def opt_sc_backend(pauli_layers):
    circuit = 0
    return circuit
    


class Paulihedral:
    
    def __init__(self, pauliIRprogram=None, ordering_method=None, backend=None, coupling_map=None):
        if pauliIRprogram == None:
            '''need to raise some error'''
            print("raise some error")
        self.pauliIRprogram = pauliIRprogram
        
        if ordering_method == 'gco' or ordering_method == 'gate-count-oriented':
            self.ordering_method = gco_ordering
        elif ordering_method == 'do' or ordering_method == 'depth-oriented':
            self.ordering_method = do_ordering
        else:
            '''need to raise some error'''
            print('unknown ordering method')
            
        if backend == 'ft' or backend == 'fault-tolerant':
            self.block_opt_method = opt_ft_backend
        elif backend == 'sc' or backend == 'superconducting':
            self.block_opt_method = opt_sc_backend
            if coupling_map == None:
                '''need to raise some error'''
                print('coupling_map missing')
            else:
                self.coupling_map = coupling_map
        else:
            '''need to raise some error'''
            print('backend not support')
        self.output_circuit = None
                
            
    def opt(self, do_max_iteration = 30):
        self.pauli_layers = self.ordering_method(self.pauliIRprogram, do_max_iteration)
        self.output_circuit = self.block_opt_method(self.pauli_layers)
        return self.output_circuit
            
    def to_circuit(self):
        if self.output_circuit == None:
            self.output_circuit = self.opt()
        return self.output_circuit
        
        

        
if __name__ == '__main__':
    pauliIRprogram = [[['IIXY', 0.1], ['IIYX', 0.1], 1],
                      [['XYII', 0.1], ['XYII', 0.1], 1]]
    test = Paulihedral(pauliIRprogram, ordering_method = 'gco')
    test.opt()
    
    
    
    
#[[XYX, p1], [XZZ, p2], ..., [ZZZ, p5], p]
#[[XYX, p1], [XZZ, p2], ..., [ZZZ, p5], p]
#[[XYX, p1], [XZZ, p2], ..., [ZZZ, p5], p]