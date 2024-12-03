import numpy as np
import edlib
from Bio.Seq import Seq
import time
import random
ref = ""
query = ""
tuples = []
best_tuples = []

def simulated_annealing(ref, query):
    initial_temperature = 5.0
    final_temperature = 0.01
    alpha = 0.92
    num_iterations = 1000

    current_state = best_tuples
    current_value = mycalc(current_state, ref, query)

    best_state = current_state
    best_value = current_value

    temperature = initial_temperature

    try:
        while temperature > final_temperature:
            for _ in range(num_iterations):
                new_state = perturb_state(current_state)
                new_value = mycalc(new_state, ref, query)
                if new_value > current_value or random.uniform(0, 0.5) < np.exp((new_value - current_value) / temperature):
                    current_state = new_state
                    current_value = new_value
                    if new_value > best_value:
                        best_value = new_value
                        best_state = new_state

            temperature *= alpha
            print("Current temperature:", temperature)
            print("Best value:", best_value)
            print("Current value:", current_value)
            print("Best state:", best_state)
    except KeyboardInterrupt:
        print("Process interrupted by user.")
        return best_state

    return best_state

def perturb_state(current_state):
    new_state = list(current_state)
    change = random.choice(range(len(new_state)))
    new_state[change] = min(max(0, new_state[change] + random.choice([-200, 200])),145063)
    return new_state
def calculate_value(tuples_str, ref, query):
    def rc(seq):
        return str(Seq(seq).reverse_complement())
    def get_points(tuples_str):
        data = []
        num = 0
        for c in tuples_str:
            if(ord('0') <= c <= ord('9')):
                num = num*10 + c - ord('0')
            elif(ord(',') == c):
                data.append(num)
                num = 0
        if(num != 0):
            data.append(num)
        return data
    def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
        A = ref[ref_st: ref_en]
        a = query[query_st: query_en]
        _a = rc(query[query_st: query_en])
        return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])
    def get_first(x):
        return x[0]
    try:
        slicepoints = np.array(get_points(tuples_str.encode()))
        if(len(slicepoints) > 0 and len(slicepoints) % 4 == 0):
            editdistance = 0
            aligned = 0
            preend = 0
            points = np.array(slicepoints).reshape((-1, 4)).tolist()
            points.sort(key = get_first)

            for onetuple in points:
                query_st, query_en, ref_st, ref_en = onetuple[0], onetuple[1], onetuple[2], onetuple[3]
                if(preend > query_st):#检测重叠m 
                    return 0
                preend = query_en
                editdistance += calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)
                aligned += query_en - query_st
            return max(aligned - editdistance - len(points) * 30, 0)#额外的惩罚碎片化的输出
        else:
            return 0
    except:
        return 0
def mycalc(tuples,ref,query,ifstr=0):
    tuples_str2 = "("
    for i in range(len(tuples)):
        t = (tuples[i] if i < len(tuples) else 0)
        tuples_str2 += str(t)
        if(i % 4 == 3):
            tuples_str2 += "),("
        else:
            tuples_str2 += ","
    tuples_str2 = tuples_str2[:-2]
    return calculate_value(tuples_str2, ref, query) if ifstr == 0 else tuples_str2
def main():
    with open('query.txt') as f:
        query = f.read().split('\n')[0]
    with open('ref.txt') as f:
        ref = f.read().split('\n')[0]
    global best_tuples
    # best_tuples = [134528, 145063, 160642200, 160652815, 128938, 134528, 160642141, 160647754, 122703, 128938, 160641429, 160647695, 117113, 122703, 160641360, 160646977, 116038, 117113, 160645817, 160646908, 109803, 116038, 160628466, 160634724, 103998, 109803, 160628202, 160634011, 102063, 103998, 160631812, 160633748, 101418, 102063, 160625620, 160626266, 100343, 101418, 160630097, 160631164, 94968, 100343, 160641326, 160646737, 92173, 94968, 160616341, 160619143, 90668, 92173, 160620380, 160621889, 90453, 90668, 160614613, 160614828, 90238, 90453, 160601965, 160602180, 87873, 90238, 160628258, 160630651, 82498, 87873, 160628376, 160633804, 76908, 82498, 160628291, 160633922, 73683, 76908, 160630610, 160633837, 60353, 73683, 160633820, 160647256, 58203, 60353, 160620535, 160622730, 56913, 58203, 160641417, 160642717, 43798, 56913, 160633765, 160646965, 38423, 43798, 160633909, 160639312, 33048, 38423, 160634044, 160639456, 28103, 33048, 160634612, 160639592, 22083, 28103, 160634131, 160640156, 17138, 22083, 160634723, 160639679, 0, 17138, 160600870, 160618084]
    best_tuples = [134528, 145063, 160642200, 160652815, 128938, 134528, 160642141, 160647754, 122703, 128938, 160641429, 160647696, 117113, 122703, 160641360, 160646977, 116038, 117113, 160645817, 160646908, 109803, 116038, 160628466, 160634726, 103998, 109803, 160628202, 160634011, 102278, 103998, 160632019, 160633748, 100773, 102278, 160624981, 160626481, 100343, 100773, 160630097, 160630525, 96258, 100343, 160642619, 160646737, 92173, 96258, 160616341, 160620437, 90668, 92173, 160620380, 160621889, 90453, 90668, 160614613, 160614828, 90238, 90453, 160601965, 160602180, 89808, 90238, 160624677, 160625107, 89378, 89808, 160613146, 160613579, 87873, 89378, 160628258, 160629776, 82498, 87873, 160628376, 160633804, 76908, 82498, 160628291, 160633922, 73468, 76908, 160630394, 160633837, 60353, 73468, 160633820, 160647049, 58203, 60353, 160620535, 160622730, 56913, 58203, 160641417, 160642717, 43798, 56913, 160633765, 160646965, 38423, 43798, 160633909, 160639312, 33048, 38423, 160634044, 160639456, 27888, 33048, 160634395, 160639592, 22083, 27888, 160634131, 160639937, 17138, 22083, 160634723, 160639679, 0, 17138, 160600871, 160618084]
    for i in range(len(best_tuples)):
        if (best_tuples[i]>160600000):
            best_tuples[i] = best_tuples[i] - 160600000
    best_state = best_tuples
    # best_state = simulated_annealing(ref, query)
    print(mycalc(best_state, ref, query, 0))
    
if __name__ == '__main__':
    pretime = time.time()
    main()
    nowtime = time.time()
    print("Time:", nowtime - pretime)