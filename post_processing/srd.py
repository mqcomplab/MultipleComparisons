import numpy as np
import random
from scipy.stats import rankdata


def read_data(sim_file):
    """Read the results of the comparisons and the similarity indices used.

    Arguments
    ---------
    sim_file : str
        Name of the .sim file with the results of the comparisons.

    Returns
    -------
    headers, content
    headers: similarity indices.
    content: results of the comparisons.
    """
    with open(sim_file, "r") as infile:
        lines = infile.readlines()
    for i in range(len(lines)):
        if "#" in lines[i]:
            start = i + 1
        elif len(lines[i].strip()) == 0:
            end = i - 1
    headers = lines[start - 1].strip().split()[1:]
    content = lines[start:end + 1]
    return headers, content


def header_modifier(headers):
    """Simplify the nomenclature of the headers (the names of the indices).

    Arguments
    ---------
    headers : list
        List of strings with the original headers as they appear in the .sim file.
    """
    new_headers = []
    for header in headers:
        index, s, d = header.split("_")
        if "1" in s:
            s = "1"
        else:
            s = "0"
        new_headers.append(index + "_" + s)
    return new_headers


def get_total_data(headers, content, ignore=[]):
    """Organize the comparison results from the data extracted from the .sim file.

    Arguments
    ---------
    headers : list
        Similarity indices.
    content : list
        List with the bulk of the results of the comparisons from the .sim file.
    ignore : list
        List of integers. These are the indices of the rows in content that won't be taken
        into account in the analysis. Default is an empty list (all data is considered).
        Used in the cross-validation analysis.

    Returns
    -------
    data : dict
        Dictionary where the keys are the similarity indices
        and the values are the results of the comparisons.
    """
    data = {}
    for header in headers:
        data[header] = []
    for i in range(len(content)):
        if i in ignore:
            pass
        else:
            comparisons = content[i].strip().split()[1:]
            for j in range(len(headers)):
                data[headers[j]].append(float(comparisons[j]))
    for header in headers:
        data[header] = np.array(data[header])
    return data


def gen_ref(content, ignore=[]):
    """Generate the reference that will be used in the SRD analysis.

    Arguments
    ---------
    content : list
        List with the values of the comparisons.
    ignore : list
        List of integers. These are the indices of the rows in content that won't be taken
        into account in the analysis. Default is an empty list (all data is considered).
        Used in the cross-validation analysis.

    Returns
    -------
    np.array
        Reference values for the SRD analysis.
    """
    ref = []
    for i in range(len(content)):
        if i in ignore:
            pass
        else:
            comparisons = content[i].strip().split()[1:]
            value = 0
            for j in range(len(comparisons)):
                value += (float(comparisons[j]))
            ref.append(value/len(comparisons))
    return np.array(ref)


def _rank_array(array):
    """Rank an array following the SRD convention.

    Arguments
    ---------
    array : {list, np.array}
        Array whose elements will be ranked.

    Returns
    -------
    Ranked array.

    Notes
    -----
    The SRD convention for ranking is as follows:
    1- Elements are ranked in increasing order (from smallest to biggest).
    2- The rank of elements with the same value is the average of the ranks.
    Example:
    _rank_array([0.1, 0.4, 0.2]) -> np.array([2, 3, 1])
    _rank_array([0.6, 0.2, 0.2, 0.1]) -> np.array([4, 2.5, 2.5, 1])
    """
    order = array.argsort()
    non_repeated_ranks = list(order.argsort() + 1)
    repeated_ranks = list(rankdata(array, method="min"))
    ranks = [0]*len(repeated_ranks)
    for i in range(len(non_repeated_ranks)):
        if repeated_ranks.count(repeated_ranks[i]) == 1:
            ranks[i] = non_repeated_ranks[i]
        else:
            indices = [j for j, x in enumerate(repeated_ranks) if x == repeated_ranks[i]]
            value = 0
            for index in indices:
                value += non_repeated_ranks[index]
            value = value/len(indices)
            for index in indices:
                ranks[index] = value
    return np.array(ranks)


def ranked_data(data):
    """Rank the results of the comparisons for each similarity index.

    Arguments
    ---------
    data : dict
        Dictionary with the results of the comparisons.

    Returns
    -------
    rk_data : dict
        Dictionary with the ranks for each index.
    """
    rk_data = {}
    for key in data:
        rk_data[key] = _rank_array(data[key])
    return rk_data


def diff_rank(rk_data, ref_rank):
    """Calculate the absolute value of the difference of the ranks of the indices and the reference.

    Arguments
    ---------
    rk_data : dict
        Dictionary with the ranks for each index.
    ref_rank: dict
        Dictionary with the ranks of the reference.

    Returns
    -------
    d_rank : dict
        Dictionary with the absolute value of the difference of the ranks of the indices and the
        reference.
    """
    d_rank = {}
    for key in rk_data:
        d_rank[key] = abs(rk_data[key] - ref_rank)
    return d_rank


def gen_srd(d_rank):
    """Calculate the SRD.

    Arguments
    ---------
    d_rank : dict
        Dictionary with the absolute value of the difference of the ranks of the indices and the
        reference.

    Returns
    -------
    srd : dict
        Dictionary with the SRD.
    """
    srd = {}
    for key in d_rank:
        srd[key] = np.sum(d_rank[key])
    return srd


def _srd_maximum(n_data):
    """Calculate the maximum possible SRD value.

    Arguments
    ---------
    n_data : int
        Number of comparisons.

    Returns
    -------
    max_srd : int
        Maximum possible SRD value.
    """
    k = (n_data - n_data % 2)/2
    if n_data % 2 == 0:
        max_srd = 2 * k ** 2
    else:
        max_srd = 2 * k * (k + 1)
    return max_srd


def scaling_srd(n_data, d_rank):
    """Scale the SRD values normalized with respect to the maximum possible value.

    Arguments
    ---------
    n_data : int
        Number of comparisons.
    d_rank : dict
        Dictionary with the absolute value of the difference of the ranks of the indices and the
        reference.

    Returns
    -------
    srd_scaled : dict
        Dictionary with the scaled SRD values.
    """
    srd = gen_srd(d_rank)
    max_srd = _srd_maximum(n_data)
    srd_scaled = {}
    for key in srd:
        srd_scaled[key] = 100 * srd[key]/max_srd
    return srd_scaled


def ranking_data_str(rk_data, d_rank, ref, ref_rank):
    """Generate string containing a (verbose) exposition of the SRD analysis.

    Arguments
    ---------
    rk_data : dict
        Dictionary with the ranks for each index.
    d_rank : dict
        Dictionary with the absolute value of the difference of the ranks of the indices and the
        reference.
    ref : np.ndarray
        Reference values for the SRD analysis.
    ref_rank: dict
        Dictionary with the ranks of the reference.

    Returns
    -------
    s : str
        String containing a verbose exposition of the SRD analysis.
    """
    s = "No         " 
    
    indices = []
    for key in sorted(rk_data):
        indices.append(key)
    
    for index in indices:
        s += "{:^{}}  ".format(index, 23)
    s += "     Reference\n            "
    for index in range(len(indices)):
        s += "Ranking        Diff      "
    s += "Reference    Ranking"
    s += "\n"
    
    for i in range(len(ref)):
        s += "{:<13d}".format(i+1)
        for index in indices:
            s += "{:^{}.1f}     ".format(rk_data[index][i], 7)
            s += "{:^{}.1f}     ".format(d_rank[index][i], 8)
        s += "{:^{}.3f}        ".format(ref[i], 5)
        s += "{:^{}.3f}     ".format(ref_rank[i], 5)
        s += "\n"
    return s


def srd_str(n_data, d_rank, srd):
    """Generate string containing a (resumed) version of the SRD analysis.

    Arguments
    ---------
    n_data : int
        Number of comparisons.
    d_rank : dict
        Dictionary with the absolute value of the difference of the ranks of the indices and the
        reference.
    srd : dict
        Dictionary with the SRD.

    Returns
    -------
    s : str
        String containing a resumed version of the SRD analysis.
    """
    
    indices = []
    for key in sorted(d_rank):
        indices.append(key)
    
    s = "{:12}".format(name.split(".")[0])
    for index in indices:
        s += "{:^{}}  ".format(index, 23)
       
    s += "\nSRD        "
    for index in indices:
        srd_value = srd[index]
        l = len(index)
        s += "{:^{}}  ".format(srd_value, 23)
    srd_scaled = scaling_srd(n_data, d_rank)
    s += "\nSRDnorm    "
    for index in indices:
        srd_scaled_value = srd_scaled[index]
        l = len(index)
        s += "{:^{}.3f}  ".format(srd_scaled_value, 23)
    return s


def cv_individual_str(n_data, d_rank, cv_type, i):
    """Generate a str with the results of a single CV study.

    Arguments
    ---------
    n_data : int
        Number of comparisons.
    d_rank : dict
        Dictionary with the absolute value of the difference of the ranks of the indices and the
        reference.
    cv_type : str
        Type of CV performed: 'sequential' or 'random'.
    i : int
        Number of the CV analysis.

    Raises
    ------
    TypeError
        If cv_type is not one of sequential or random.

    Returns
    -------
    s : str
        String containing the results of the i-th CV analysis of the given cv_type.
    """
    if cv_type not in ["random", "sequential"]:
        raise TypeError("cv_type can only be one of 'sequential' or 'random'.")
    s = ""
    if cv_type == "sequential":
        s_help = "Gr{}_A".format(i+1)
    elif cv_type == "random":
        s_help = "Gr{}_B".format(i+1)
    else:
        raise TypeError("cv_type can only be one of 'sequential' or 'random'.")
    s += "{:<11}".format(s_help)
    
    indices = []
    for key in sorted(d_rank):
        indices.append(key)
       
    srd_scaled = scaling_srd(n_data, d_rank)
    for index in indices:
        srd_scaled_value = srd_scaled[index]
        l = len(index)
        s += "{:^{}.3f}  ".format(srd_scaled_value, 23)
    s += "\n"
    return s


def gen_cross_indices(n_data, fraction=5, cv_type="random", repetitions=7):
    """Generate the indices of the rows that will be used as test set in the CV analysis.

    Arguments
    ---------
    n_data : int
        Number of comparisons.
    fraction : int
        Fraction of the data that will be ignored (used as test set).
    cv_type : str
        Type of CV performed: 'sequential' or 'random'.
    repetitions : int
        Number of times that a random CV will be performed.

    Raises
    ------
    TypeError
        If cv_type is not one of sequential or random.

    Returns
    -------
    c_indices : list
        List containing the indices of the rows that will be ignored (used as test set)
        in each of the CV analysis.
    """
    if cv_type not in ["random", "sequential"]:
        raise TypeError("cv_type can only be one of 'sequential' or 'random'.")
    ignore_size = n_data//fraction
    c_indices = []
    if cv_type == "random":
        for i in range(repetitions):
            c_inds = []
            while len(c_inds) < ignore_size:
                index = random.randint(0, n_data - 1)
                if index in c_inds:
                    pass
                else:
                    c_inds.append(index)
            c_indices.append(c_inds)
    elif cv_type == "sequential":
        for i in range(fraction):
            c_indices.append(list(range(n_data)[i * ignore_size:(i + 1) * ignore_size]))
    return c_indices


def cv_total_str(headers, content, n_data, fraction=5, cv_type="random", repetitions=7):
    """Generate string with the results of the CV studies.

    Arguments
    ---------
    headers : list
        List containing the similarity indices considered.
    content : list
        List with the values of the comparisons.
    n_data : int
        Number of comparisons.
    fraction : int
        Fraction of the data that will be ignored (used as test set).
    cv_type : str
        Type of CV performed: 'sequential' or 'random'.
    repetitions : {int, None}
        Number of times that a random CV will be performed.

    Raises
    ------
    TypeError
        If cv_type is not one of sequential or random.

    Returns
    -------
    s : str
        String containing the results of all the CV studies.
    """
    if cv_type not in ["random", "sequential"]:
        raise TypeError("cv_type can only be one of 'sequential' or 'random'.")
    indices = []
    for header in headers:
        indices.append(header)
    
    s = ""
    c_indices = gen_cross_indices(n_data, fraction, cv_type, repetitions)
    for i, ignore in enumerate(c_indices):
        data = get_total_data(headers, content, ignore)
        ref = gen_ref(content, ignore)
        n_data = len(ref)
        ref_rank = _rank_array(ref)
        rk_data = ranked_data(data)
        d_rank = diff_rank(rk_data, ref_rank)
        s += cv_individual_str(n_data, d_rank, cv_type, i)
    return s


def output_header(fp_total, m, n):
    """Generate the header of the (verbose) output.

    Arguments
    ---------
    fp_total : int
        Total number of fingerprints.
    m : int
        Length of the fingerprints.
    n : int
        Number of fingerprints compared simultaneously.

    Returns
    -------
    s : str
        Header of the verbose output.
    """
    s = "Sum of Ranking Differences Analysis\n\n"
    s += "Fingerprint size (m)\n" + str(m) + "\n\n"
    s += "Total number of fingerprints\n" + str(fp_total) + "\n\n"
    s += "Fingerprints compared simultaneously (n)\n" + str(n) + "\n\n"
    return s


def extract_name_data(filename):
    """Extract the calculation information from the name of the .sim file.

    Arguments
    ---------
    filename : str
        Name of the .sim file.

    Returns
    -------
    weight_info, fp_total, m, n
    weight_info : str
        'w' or 'nw'.
    fp_total : int
        Total number of fingerprints.
    m : int
        Length of the fingerprints.
    n : int
        Number of fingerprints compared simultaneously.
    """
    weight_info = filename.split("FP")[0]
    data = filename.split("FP")[-1]
    fp_total = data.split("m")[0]
    m = data.split(".")[0].split("m")[-1].split("n")[0]
    n = data.split("n")[-1].split(".")[0]
    return weight_info, fp_total, m, n


def process_data(sim_file, fraction=7, repetitions=50):
    """Perform SRD analysis.

    Arguments
    ---------
    sim_file : str
        Name of the .sim file
    fraction : int
        Fraction of the data that will be ignored (used as test set).
    repetitions : int
        Number of times that a random CV will be performed.

    Creates the files with the results of the SRD analysis.
    """
    ignore = []
    headers, content = read_data(sim_file)
    headers = header_modifier(headers)
    data = get_total_data(headers, content, ignore)
    ref = gen_ref(content, ignore)
    n_data = len(ref)
    ref_rank = _rank_array(ref)
    rk_data = ranked_data(data)
    d_rank = diff_rank(rk_data, ref_rank)
    srd = gen_srd(d_rank)
    srd_results = srd_str(n_data, d_rank, srd)

    cv_results_seq = cv_total_str(headers, content, n_data, fraction,
                                  cv_type="sequential", repetitions=None)
    cv_results_rand = cv_total_str(headers, content, n_data, fraction,
                                   cv_type="random", repetitions=repetitions)

    s = srd_results + "\n" + cv_results_seq + cv_results_rand
    
    name = "SRD" + sim_file.split(".")[0] + "CV" + str(fraction) + ".txt"
    with open(name, "w") as outfile:
        outfile.write(s)


if __name__ == "__main__":
    # Performs the srd analysis for all the .sim files in the current folder.
    import glob
    fraction = 7
    repetitions = 50
    sim_list = glob.glob("*.sim")
    for name in sim_list:
        process_data(name, fraction, repetitions)
