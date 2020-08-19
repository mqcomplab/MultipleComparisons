import numpy as np
from itertools import combinations
from scipy.special import binom
from indices.indices_info import Indices


def _file_name_gen(key):
    """Generate the name of the file containing a given index.

    Arguments
    ---------
    key : str
        Key of the Indices dict that corresponds to a given index.

    Returns
    -------
    name : str
        Name of the file containing the class corresponding to the given index.
    """
    parts = key.split("-")
    if len(parts) == 1:
        name = key.lower()
    else:
        name = ""
        for part in parts[:-1]:
            name += part.lower()+"_"
        name += parts[-1].lower()
    return name


def generate_bitstring(size):
    """Generate a random fingerprint.

    Arguments
    ---------
    size : int
        Size (or length) of the fingerprint.
    
    Returns
    -------
    fingerprint : np.ndarray
        Fingerprint as a numpy array.
    """
    if not isinstance(size, int):
        raise TypeError("size can only be an integer.")
    return np.random.choice([0, 1], size=size)


def gen_fingerprints(fp_total, fp_size):
    """Generates random fingerprints.

    Arguments
    ---------
    fp_total : int
        Number of fingerprints.
    fp_size : int
        Size (or length) of the fingerprints.

    Returns
    -------
    total_fingerprints : np.array
        Numpy array containing the fingerprints.
    """
    total_fingerprints = []
    for i in range(fp_total):
        total_fingerprints.append(generate_bitstring(fp_size))
    total_fingerprints = np.array(total_fingerprints)
    return total_fingerprints


def calc_indices(indices=Indices, fp_total=2,
                 total_fingerprints=np.array([np.array([1]), np.array([1])]),
                 n=2, c_threshold=None, w_factor="fraction"):
    """Calculate the indices and generates the output.

    Arguments
    ---------
    indices : dict
        Dictionary with the indices that will be calculated.
    fp_total : int
        Total number of fingerprints.
    total_fingerprints : np.ndarray
        Numpy array containing the fingerprints that will be compared.
    n : int
        Number of fingerprints that will be compared simultaneously.
    c_threshold : {None, 'dissimilar', int}
        Coincidence threshold.
    w_factor : {"fraction", "power_n"}
        Type of weight function that will be used.

    Raises
    ------
    TypeError
        If n is not an integer.
    ValueError
        If the number of fingerprints in total_fingerprints is not equal to fp_total.
        If n is less than 2.
        If n is greater than f_total.

    Returns
    -------
    Results : dict
        Dictionary with the results of the comparisons.
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer.")
    if len(total_fingerprints) != fp_total:
        raise ValueError("The number of fingerprints in total_fingerprints must be equal"
                         "to fp_total.")
    if n < 2:
        raise ValueError("n cannot be less than 2.")
    if n > fp_total:
        raise ValueError("n cannot be greater than fp_total.")

    # Dictionary that will contain the results of all the comparisons.
    # Its structure is: Results[index] = (class_name, [index_values])
    Results = {}
    
    for s_index in indices:
        for variant in Indices[s_index][2]:
            Results[indices[s_index][1] + "_" + variant] = (Indices[s_index][0], [])

    # Sets of n numbers that indicate which fingerprints will be compared at a given time.
    index_list = list(combinations(range(fp_total), n))
    
    # Populating the Results dict with the results of the comparisons.
    for inds in index_list:
        fingerprints = total_fingerprints[list(inds)]
        for s_index in sorted(Results):
            h = "fingerprints=np.array([np.array([1, 0]), np.array([0, 1])]), "
            h += "c_threshold = 1, "
            h += "w_factor=None"
            h = "(" + h + ")"
            exec("index = " + Results[s_index][0] + h, None, globals())
            index.__init__(fingerprints=fingerprints, c_threshold=c_threshold, w_factor=w_factor)
            exec("result = index." + s_index, None, globals())
            Results[s_index][1].append(result)
    return Results


def _indices_values(results, fp_total=2, n=2, methods=[]):
    """Generate output with the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    fp_total : int
        Total number of fingerprints.
    n : int
        Number of fingerprints that will be compared simultaneously.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the results of the comparisons of the given indices.
    """
    s = ""
    for method in methods:
        s += method + "      "
    s += "\n"
    for i in range(int(binom(fp_total, n))):
        s += "{:<13d}".format(i + 1)
        for method in methods:
            l = len(method)
            s += "{:^{}.6f}     ".format(results[method][1][i], l + 1)
        s += "\n"
    s += "\n             "
    for method in methods:
        s += method + "      "
    return s


def _max(results, methods=[]):
    """Generate the maxima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the maxima of the comparisons of the given indices.
    """
    s = ""
    s += "\nMax          "
    for method in methods:
        max = np.amax(np.array(results[method][1]))
        l = len(method)
        s += "{:^{}.6f}     ".format(max, l + 1)
    return s


def _min(results, methods=[]):
    """Generate the minima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the minima of the comparisons of the given indices.
    """
    s = ""
    s += "\nMin          "
    for method in methods:
        min = np.amin(np.array(results[method][1]))
        l = len(method)
        s += "{:^{}.6f}     ".format(min, l + 1)
    return s


def _abs_max(results, methods=[]):
    """Generate the abs maxima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the abs maxima of the comparisons of the given indices.
    """
    s = ""
    s += "\nAbsMax       "
    for method in methods:
        absmax = np.amax(np.abs(np.array(results[method][1])))
        l = len(method)
        s += "{:^{}.6f}     ".format(absmax, l + 1)
    return s


def _abs_min(results, methods=[]):
    """Generate the abs minima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the abs minima of the comparisons of the given indices.
    """
    s = ""
    s += "\nAbsMin       "
    for method in methods:
        absmin = np.amin(np.abs(np.array(results[method][1])))
        l = len(method)
        s += "{:^{}.6f}     ".format(absmin, l + 1)
    return s


def _average(results, methods=[]):
    """Generate the averages of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the averages of the comparisons of the given indices.
    """
    s = ""
    s += "\nAverage      "
    for method in methods:
        average = np.average(np.array(results[method][1]))
        l = len(method)
        s += "{:^{}.6f}     ".format(average, l + 1)
    return s


def _abs_average(results, methods=[]):
    """Generate the abs averages of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    String with the abs averages of the comparisons of the given indices.
    """
    s = ""
    s += "\nAbsAverage   "
    for method in methods:
        absaverage = np.average(np.abs(np.array(results[method][1])))
        l = len(method)
        s += "{:^{}.6f}     ".format(absaverage, l + 1)
    return s


def indices_output(results, fp_total=2, fp_size=1, n=2):
    """Generate output file with the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    fp_total : int
        Total number of fingerprints.
    fp_size : int
        Size of the fingerprints.
    n : int
        Number of fingerprints that will be compared simultaneously.
    """
    # Generic header for the output file.
    s = "Similarity analysis\n\n"
    s += "Fingerprint size (m)\n" + str(fp_size) + "\n\n"
    s += "Total number of fingerprints\n" + str(fp_total) + "\n\n"
    s += "Fingerprints compared simultaneously (n)\n" + str(n) + "\n\n"
    s += "#            " 
    
    # Weighted indices.
    w = []
    # Unweighted indices.
    no_w = []
    for s_index in sorted(results):
        if "w" in s_index:
            w.append(s_index)
        else:
            no_w.append(s_index)
    r_w = s
    r_no_w = s

    r_w += _indices_values(results=results, fp_total=fp_total, n=n, methods=w)
    r_no_w += _indices_values(results=results, fp_total=fp_total, n=n, methods=no_w)
    
    r_w += _max(results=results, methods=w)
    r_no_w += _max(results=results, methods=no_w)

    r_w += _abs_max(results=results, methods=w)
    r_no_w += _abs_max(results=results, methods=no_w)

    r_w += _min(results=results, methods=w)
    r_no_w += _min(results=results, methods=no_w)

    r_w += _abs_min(results=results, methods=w)
    r_no_w += _abs_min(results=results, methods=no_w)

    r_w += _average(results=results, methods=w)
    r_no_w += _average(results=results, methods=no_w)

    r_w += _abs_average(results=results, methods=w)
    r_no_w += _abs_average(results=results, methods=no_w)

    with open("wFP"+str(fp_total)+"m"+str(fp_size)+"n"+str(n)+".sim", "w") as outfile:
        outfile.write(r_w)

    with open("nwFP"+str(fp_total)+"m"+str(fp_size)+"n"+str(n)+".sim", "w") as outfile:
        outfile.write(r_no_w)


if __name__ == "__main__":
    # Imports the classes corresponding to the similarity indices.
    for key in Indices:
        s = "from indices." + _file_name_gen(key) + " import " + Indices[key][0]
        exec(s)

    # Sample run with randomly generated fingerprints.

    # Coincidence threshold.
    c_threshold = None

    # Weight factor.
    w_factor = "fraction"

    # Fingerprint sizes. The analysis will be ran for every given size.
    fp_sizes = [1000]

    # Total fingerprint numbers. The analysis will be ran for every set of fingerprints.
    fp_totals = [100]

    for fp_size in fp_sizes:
        for fp_total in fp_totals:
            total_fingerprints = gen_fingerprints(fp_total, fp_size)
            for n in [2, fp_totals[0]]:
                # n values. Possible values are 2 <= n <= fp_total.
                results = calc_indices(indices=Indices, fp_total=fp_total,
                                       total_fingerprints=total_fingerprints,
                                       n=n, c_threshold=c_threshold, w_factor=w_factor)
                indices_output(results=results, fp_total=fp_total, fp_size=fp_size, n=n)
