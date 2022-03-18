import numpy as np
import random
import glob
import pickle
from math import log, ceil

# ECS_MeDiv algorithm

def calculate_counters(data_sets, c_threshold=None, w_factor="fraction"):
    """Calculate 1-similarity, 0-similarity, and dissimilarity counters

    Arguments
    ---------
    data_sets : np.ndarray
        Array of arrays. Each sub-array contains m + 1 elements,
        with m being the length of the fingerprints. The first
        m elements are the column sums of the matrix of fingerprints.
        The last element is the number of fingerprints.

    c_threshold : {None, 'dissimilar', int, float}
        Coincidence threshold.
        None : Default, c_threshold = n_fingerprints % 2
        'dissimilar' : c_threshold = ceil(n_fingerprints / 2)
        int : Integer number < n_fingerprints
        float: Real number in the (0, 1) interval, c_threshold *= n_fingerprints

    w_factor : {"fraction", "power_n"}
        Type of weight function that will be used.
        'fraction' : similarity = d[k]/n
                     dissimilarity = 1 - (d[k] - n_fingerprints % 2)/n_fingerprints
        'power_n' : similarity = n**-(n_fingerprints - d[k])
                    dissimilarity = n**-(d[k] - n_fingerprints % 2)
        other values : similarity = dissimilarity = 1

    Returns
    -------
    counters : dict
        Dictionary with the weighted and non-weighted counters.

    Notes
    -----
    Please, cite the original papers on the n-ary indices:
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00505-3
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00504-4
    """
    # Setting matches
    total_data = np.sum(data_sets, axis=0)
    n_fingerprints = int(total_data[-1])
    c_total = total_data[:-1]
    
    # Assign c_threshold
    if not c_threshold or c_threshold == 'min':
        c_threshold = n_fingerprints % 2
    if isinstance(c_threshold, str):
        if c_threshold != 'dissimilar':
            raise TypeError("c_threshold must be None, 'dissimilar', or an integer.")
        else:
            c_threshold = ceil(n_fingerprints / 2)
    if isinstance(c_threshold, int):
        if c_threshold >= n_fingerprints:
            raise ValueError("c_threshold cannot be equal or greater than n_fingerprints.")
        c_threshold = c_threshold
    if 0 < c_threshold < 1:
        c_threshold *= n_fingerprints
    
    # Set w_factor
    if w_factor:
        if "power" in w_factor:
            power = int(w_factor.split("_")[-1])
            def f_s(d):
                return power**-float(n_fingerprints - d)
    
            def f_d(d):
                return power**-float(d - n_fingerprints % 2)
        elif w_factor == "fraction":
            def f_s(d):
                return d/n_fingerprints
    
            def f_d(d):
                return 1 - (d - n_fingerprints % 2)/n_fingerprints
        else:
            def f_s(d):
                return 1
    
            def f_d(d):
                return 1
    else:
        def f_s(d):
            return 1
    
        def f_d(d):
            return 1
    
    # Calculate a, d, b + c
    a = 0
    w_a = 0
    d = 0
    w_d = 0
    total_dis = 0
    total_w_dis = 0
    for s in c_total:
        if 2 * s - n_fingerprints > c_threshold:
            a += 1
            w_a += f_s(2 * s - n_fingerprints)
        elif n_fingerprints - 2 * s > c_threshold:
            d += 1
            w_d += f_s(abs(2 * s - n_fingerprints))
        else:
            total_dis += 1
            total_w_dis += f_d(abs(2 * s - n_fingerprints))
    total_sim = a + d
    total_w_sim = w_a + w_d
    p = total_sim + total_dis
    w_p = total_w_sim + total_w_dis
    
    counters = {"a": a, "w_a": w_a, "d": d, "w_d": w_d,
                "total_sim": total_sim, "total_w_sim": total_w_sim,
                "total_dis": total_dis, "total_w_dis": total_w_dis,
                "p": p, "w_p": w_p}
    
    return counters
    
def gen_sim_dict(data_sets, c_threshold=None, w_factor="fraction"):
    counters = calculate_counters(data_sets, c_threshold=c_threshold, w_factor="fraction")
    # Indices
    # AC: Austin-Colwell, BUB: Baroni-Urbani-Buser, CTn: Consoni-Todschini n
    # Fai: Faith, Gle: Gleason, Ja: Jaccard, Ja0: Jaccard 0-variant
    # JT: Jaccard-Tanimoto, RT: Rogers-Tanimoto, RR: Russel-Rao
    # SM: Sokal-Michener, SSn: Sokal-Sneath n

    # Weighted Indices
    ac_w = (2/np.pi) * np.arcsin(np.sqrt(counters['total_w_sim']/
                                         counters['w_p']))
    bub_w = ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'])/\
            ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'] + counters['total_w_dis'])
    ct1_w = (log(1 + counters['w_a'] + counters['w_d']))/\
            (log(1 + counters['w_p']))
    ct2_w = (log(1 + counters['w_p']) - log(1 + counters['total_w_dis']))/\
            (log(1 + counters['w_p']))
    ct3_w = (log(1 + counters['w_a']))/\
            (log(1 + counters['w_p']))
    ct4_w = (log(1 + counters['w_a']))/\
            (log(1 + counters['w_a'] + counters['total_w_dis']))
    fai_w = (counters['w_a'] + 0.5 * counters['w_d'])/\
            (counters['w_p'])
    gle_w = (2 * counters['w_a'])/\
            (2 * counters['w_a'] + counters['total_w_dis'])
    ja_w = (3 * counters['w_a'])/\
           (3 * counters['w_a'] + counters['total_w_dis'])
    ja0_w = (3 * counters['total_w_sim'])/\
            (3 * counters['total_w_sim'] + counters['total_w_dis'])
    jt_w = (counters['w_a'])/\
           (counters['w_a'] + counters['total_w_dis'])
    rt_w = (counters['total_w_sim'])/\
           (counters['w_p'] + counters['total_w_dis'])
    rr_w = (counters['w_a'])/\
           (counters['w_p'])
    sm_w =(counters['total_w_sim'])/\
          (counters['w_p'])
    ss1_w = (counters['w_a'])/\
            (counters['w_a'] + 2 * counters['total_w_dis'])
    ss2_w = (2 * counters['total_w_sim'])/\
            (counters['w_p'] + counters['total_w_sim'])


    ## Non-Weighted Indices
    ac_nw = (2/np.pi) * np.arcsin(np.sqrt(counters['total_w_sim']/
                                          counters['p']))
    bub_nw = ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'])/\
             ((counters['a'] * counters['d'])**0.5 + counters['a'] + counters['total_dis'])
    ct1_nw = (log(1 + counters['w_a'] + counters['w_d']))/\
             (log(1 + counters['p']))
    ct2_nw = (log(1 + counters['w_p']) - log(1 + counters['total_w_dis']))/\
             (log(1 + counters['p']))
    ct3_nw = (log(1 + counters['w_a']))/\
             (log(1 + counters['p']))
    ct4_nw = (log(1 + counters['w_a']))/\
             (log(1 + counters['a'] + counters['total_dis']))
    fai_nw = (counters['w_a'] + 0.5 * counters['w_d'])/\
             (counters['p'])
    gle_nw = (2 * counters['w_a'])/\
             (2 * counters['a'] + counters['total_dis'])
    ja_nw = (3 * counters['w_a'])/\
            (3 * counters['a'] + counters['total_dis'])
    ja0_nw = (3 * counters['total_w_sim'])/\
             (3 * counters['total_sim'] + counters['total_dis'])
    jt_nw = (counters['w_a'])/\
            (counters['a'] + counters['total_dis'])
    rt_nw = (counters['total_w_sim'])/\
            (counters['p'] + counters['total_dis'])
    rr_nw = (counters['w_a'])/\
            (counters['p'])
    sm_nw =(counters['total_w_sim'])/\
           (counters['p'])
    ss1_nw = (counters['w_a'])/\
             (counters['a'] + 2 * counters['total_dis'])
    ss2_nw = (2 * counters['total_w_sim'])/\
             (counters['p'] + counters['total_sim'])

    # Dictionary with all the results
    Indices = {'nw': {'AC':ac_nw,
                      'BUB':bub_nw,
                      'CT1':ct1_nw,
                      'CT2':ct2_nw,
                      'CT3':ct3_nw,
                      'CT4':ct4_nw,
                      'Fai':fai_nw,
                      'Gle':gle_nw,
                      'Ja0':ja0_nw,
                      'Ja':ja_nw,
                      'JT':jt_nw,
                      'RT':rt_nw,
                      'RR':rr_nw,
                      'SM':sm_nw,
                      'SS1':ss1_nw,
                      'SS2':ss2_nw},
                'w': {'AC':ac_w,
                      'BUB':bub_w,
                      'CT1':ct1_w,
                      'CT2':ct2_w,
                      'CT3':ct3_w,
                      'CT4':ct4_w,
                      'Fai':fai_w,
                      'Gle':gle_w,
                      'Ja0':ja0_w,
                      'Ja':ja_w,
                      'JT':jt_w,
                      'RT':rt_w,
                      'RR':rr_w,
                      'SM':sm_w,
                      'SS1':ss1_w,
                      'SS2':ss2_w}}
    return Indices

def calculate_medoid(total_data, n_ary = 'RR', weight = 'nw'):
    """Calculate the medoid of a set"""
    index = len(total_data[0]) + 1
    min_sim = 3.08
    total_sum = np.sum(total_data, axis = 0)
    for i, pixel in enumerate(total_data):
        i_sum = total_sum - total_data[i]
        data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(data_sets)
        sim_index = Indices[weight][n_ary]
        if sim_index < min_sim:
            min_sim = sim_index
            index = i
        else:
            pass
    return index
    
def calculate_outlier(total_data, n_ary = 'RR', weight = 'nw'):
    """Calculate the outlier of a set"""
    index = len(total_data[0]) + 1
    max_sim = -3.08
    total_sum = np.sum(total_data, axis = 0)
    for i, pixel in enumerate(total_data):
        i_sum = total_sum - total_data[i]
        data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(data_sets)
        sim_index = Indices[weight][n_ary]
        if sim_index > max_sim:
            max_sim = sim_index
            index = i
        else:
            pass
    return index

def get_single_index(total_data, indices, selected_n, c_threshold=None, n_ary = 'RR', weight = 'nw'):
    """Binary tie-breaker selection criterion"""
    index = len(total_data[0]) + 1
    min_value = 3.08
    for i in indices:
        v = 0
        for j in selected_n:
            c_total = total_data[j] + total_data[i]
            data_sets = [np.append(c_total, 2)]
            Indices = gen_sim_dict(data_sets, c_threshold=c_threshold)
            sim_index = Indices[weight][n_ary]
            v += sim_index
        av_v = v/(len(selected_n) + 1)
        if av_v < min_value:
            index = i
            min_value = av_v
    return index

def get_new_index_n(total_data, selected_condensed, n, select_from_n, selected_n, c_threshold=None,
                    n_ary = 'RR', weight = 'nw'):
    """Select a diverse object using the ECS_MeDiv algorithm"""
    n_total = n + 1
    # min value that is guaranteed to be higher than all the comparisons
    min_value = 3.08
    
    # placeholder index
    indices = [len(total_data[0]) + 1]
    
    # for all indices that have not been selected
    for i in select_from_n:
        # column sum
        c_total = selected_condensed + total_data[i]
        # calculating similarity
        data_sets = [np.append(c_total, n_total)]
        Indices = gen_sim_dict(data_sets, c_threshold=c_threshold)
        sim_index = Indices[weight][n_ary]
        # if the sim of the set is less than the similarity of the previous diverse set, update min_value and index
        if sim_index < min_value:
            indices = [i]
            min_value = sim_index
        elif sim_index == min_value:
            indices.append(i)
    if len(indices) == 1:
        index = indices[0]
    else:
        # Use average of binary similarities as tie-breaker
        index = get_single_index(total_data, indices, selected_n, c_threshold=None, n_ary = n_ary, weight = 'nw')
    return index


n_arys = ['AC', 'BUB', 'CT1', 'CT2', 'CT3', 'CT4', 'Fai',
          'Gle', 'Ja', 'Ja0', 'JT', 'RT', 'RR', 'SM', 'SS1', 'SS2']

for n_ary in n_arys:
    DiversityDict = {}
    for c_threshold in ['min']:
        for file in glob.glob('*.npy'):
            if 'rank' in file:
                pass
            else:
                base_name = file.split('.')[0]
                # Seed selection:
                # medoid = start from medoid
                # random = select random initial seed
                # out = start from outlier
                start = 'medoid'
                # Numpy array with the data
                total_data = np.load(file)
                total_n = []
                
                # total number of fingerprints
                fp_total = len(total_data)
                
                # indices of all the fingerprints
                total_indices = np.array(range(fp_total))
                
                # starting point
                if start =='medoid':
                    seed = calculate_medoid(total_data, n_ary = n_ary)
                elif start == 'random':
                    seed = random.randint(0, fp_total - 1)
                elif start == 'out':
                    seed = calculate_outlier(total_data, n_ary = n_ary)
                else:
                    print('Select a correct starting point')
                selected_n = [seed]
                
                # vector with the column sums of all the selected fingerprints
                selected_condensed = total_data[seed]
                
                # number of fingerprints selected
                n = 1
                while len(selected_n) < 10:
                    # indices from which to select the new fingerprints
                    select_from_n = np.delete(total_indices, selected_n)
                    
                    # new index selected
                    new_index_n = get_new_index_n(total_data, selected_condensed, n, select_from_n, selected_n,
                                                  c_threshold=c_threshold, n_ary = n_ary)
                    
                    # updating column sum vector
                    selected_condensed += total_data[new_index_n]
                    
                    # updating selected indices
                    selected_n.append(new_index_n)
                    
                    # updating n
                    n = len(selected_n)
                    print(selected_n)
                DiversityDict[c_threshold][base_name] = selected_n
    f = open(n_ary + '_' + start + '_' + base_name + '_.pkl', 'wb')
    pickle.dump(DiversityDict, f)
    f.close()
