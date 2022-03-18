import glob
from matplotlib import pyplot as plt
from srd import header_modifier, extract_name_data


def read_stat_data(sim_file):
    """Read the summary of the statistical results of the comparisons,
    and the similarity indices used.

    Arguments
    ---------
    sim_file : str
        Name of the .sim file with the results of the comparisons.

    Returns
    -------
    headers, content
    headers: similarity indices.
    content: summary of the statistical results of the comparisons.
    """
    with open(sim_file, "r") as infile:
        lines = infile.readlines()
    for i in range(len(lines)):
        if len(lines[i].strip()) == 0:
            start = i + 2
        elif "AbsAverage" in lines[i]:
            end = i
        elif "#" in lines[i]:
            headers_line = i
    headers = lines[headers_line].strip().split()[1:]
    content = lines[start: end + 1]
    return headers, content


def process_stat_content(content):
    """Process the statistical results of the comparisons.


    Arguments
    ---------
    content : list
        List containing the raw statistical results of the comparisons.

    Returns
    -------
    max_data, abs_max_data, min_data, abs_min_data, average_data, abs_average_data
    max_data: maximum value of the index.
    abs_max_data: maximum absolute value of the index.
    min_data: minimum value of the index.
    abs_min_data: minimum absolute value of the index.
    average_data: average value of the index.
    abs_average_data: average of the absolute values of the index.
    """
    max_data = content[0].strip().split()[1:]
    abs_max_data = content[1].strip().split()[1:]
    min_data = content[2].strip().split()[1:]
    abs_min_data = content[3].strip().split()[1:]
    average_data = content[4].strip().split()[1:]
    abs_average_data = content[5].strip().split()[1:]
    return max_data, abs_max_data, min_data, abs_min_data, average_data, abs_average_data


def _gen_ind_order(new_headers):
    """Generate auxiliary dictionary with the used indices.

    Arguments
    ---------
    new_headers : list
        List containing the modified headers corresponding to each similarity index.

    Returns
    -------
    ind_order : dict
        Dictionary used to order the similarity indices.
    """
    ind_order = {}
    for i, s_index in enumerate(new_headers):
        ind_order[i + 1] = s_index
    return ind_order


def gen_output_file(new_headers, stat_dict, n_values, name):
    """Generate the summary file for a given statistical measure.

    Arguments
    ---------
    new_headers : list
        List containing the modified headers corresponding to each similarity index.
    stat_dict : dict
        Dictionary with the values of the statistical values for different values of n.
    n_values : list
        List with the studied values of n (number of fingerprints compared simultaneously).
    name : str
        Name of the output .dat file.

    Generate the .dat file for the given statistical measure.
    """
    s = "       "
    for header in new_headers:
        s += "{:^22}      ".format(header)
    s += "\n     "
    for header in new_headers:
        s += "       w         nw         "
    s += "\n"
    for n in n_values:
        s += "n{:<4}".format(n)
        for header in new_headers:
            s += "{:>10}  {:>10}      ".format(stat_dict[header]["w"][n],
                                               stat_dict[header]["nw"][n])
        s += "\n"
    with open(name, "w") as outfile:
        outfile.write(s)
    

def process_sim_files():
    """Process the .sim files in the current folder.

    Returns
    -------
    n_values, max, abs_max, min, abs_min, average, abs_average
    n_values: List with the n values of the .sim files in the current folder.
    max: Dictionary with the maximum values of the indices.
    abs_max: Dictionary with the maximum of the absolute values of the indices.
    min: Dictionary with the minimum values of the indices.
    abs_min: Dictionary with the minimum of the absolute values of the indices.
    average: Dictionary with the average values of the indices.
    abs_average: Dictionary with the average of the absolute values of the indices.
    """
    file_list = glob.glob("*.sim")
    max = {}
    abs_max = {}
    min = {}
    abs_min = {}
    average = {}
    abs_average = {}
    dict_list = [max, min, abs_max, abs_min, average, abs_average]
    n_values = []

    headers, _ = read_stat_data(file_list[0])
    new_headers = header_modifier(headers)
    ind_order = _gen_ind_order(new_headers)

    for header in new_headers:
        for dict in dict_list:
            dict[header] = {"w": {}, "nw": {}}

    for name in file_list:
        weight_info, fp_total, m, n = extract_name_data(name)
        if int(n) not in n_values:
            n_values.append(int(n))

    n_values = sorted(n_values)

    for name in file_list:
        weight_info, fp_total, m, n = extract_name_data(name)
        _, content = read_stat_data(name)
        max_data, abs_max_data, min_data, abs_min_data, average_data, abs_average_data = \
            process_stat_content(content)
        for i, data in enumerate(max_data, start=1):
            max[ind_order[i]][weight_info][int(n)] = data
        for i, data in enumerate(abs_max_data, start=1):
            abs_max[ind_order[i]][weight_info][int(n)] = data
        for i, data in enumerate(min_data, start=1):
            min[ind_order[i]][weight_info][int(n)] = data
        for i, data in enumerate(abs_min_data, start=1):
            abs_min[ind_order[i]][weight_info][int(n)] = data
        for i, data in enumerate(average_data, start=1):
            average[ind_order[i]][weight_info][int(n)] = data
        for i, data in enumerate(abs_average_data, start=1):
            abs_average[ind_order[i]][weight_info][int(n)] = data
    return n_values, max, abs_max, min, abs_min, average, abs_average


def plot_index(index_stat_dict, n_values, name):
    """Generate the plot for the given statistical measure for a similarity index.

    Arguments
    ---------
    index_stat_dict : dict
        Dictionary with the value of the given statistical measure for a given index.
    n_values : list
        List with the studied values of n (number of fingerprints compared simultaneously).
    name : str
        Name of the output file.

    Generate the .png file plotting the given statistical measure for a given index
    with respect to the different values of n, for the weighted and unweighted cases.
    """
    w_data = []
    nw_data = []
    for i in sorted(index_stat_dict["w"]):
        w_data.append(float(index_stat_dict["w"][i]))
        nw_data.append(float(index_stat_dict["nw"][i]))
    plt.style.use('seaborn-poster')
    plt.plot(n_values, w_data, label="w")
    plt.plot(n_values, nw_data, label="nw")
    plt.xlabel("n_values")
    plt.ylabel("values")
    plt.xticks(n_values)
    plt.title("{}".format(name))
    plt.legend()
    plt.tight_layout()
    plt.savefig("{}.png".format(name))


def plot_dict(full_stat_dict, n_values, name):
    """Generate the plot for the given statistical measure for all the similarity indices.

    Arguments
    ---------
    full_stat_dict : dict
        Dictionary with the value of the given statistical measure for all the indices.
    n_values : list
        List with the studied values of n (number of fingerprints compared simultaneously).
    name : str
        Root of the name of the output file.

    Generate the .png files plotting the given statistical measure for all the indices
    with respect to the different values of n, for the weighted and unweighted cases.
    """
    for s_index in full_stat_dict:
        plt.close()
        plot_index(full_stat_dict[s_index], n_values, name=s_index + name)


if __name__ == "__main__":
    # Generates the .dat files with the summary of the statistical measures for all the indices.
    # It assumes that all the .sim files in the folder correspond to the same values of:
    # FP (fp_total): total number of fingerprints.
    # m : fingerprint length.
    # The .sim files will have different values of n (fingerprints compared simultaneously).
    file_list = glob.glob("*.sim")
    if len(file_list) == 0:
        raise TypeError("There are no .sim files in the current folder.")
    weight_info, fp_total, m, n = extract_name_data(file_list[0])
    headers, _ = read_stat_data(file_list[0])
    new_headers = header_modifier(headers)
    n_values, max, abs_max, min, abs_min, average, abs_average = process_sim_files()
    stat_measures = {"max": (max, "Max"), "abs_max": (abs_max, "AbsMax"),
                     "min": (min, "Min"), "abs_min": (abs_min, "AbsMin"),
                     "average": (average, "Average"), "abs_average": (abs_average, "AbsAverage")}
    for stat_measure in stat_measures:
        gen_output_file(new_headers, stat_measures[stat_measure][0], n_values,
                        stat_measures[stat_measure][1]+"FP{}m{}.dat".format(fp_total, m))

    # Plotting.
    # In this particular example we plot the average of the absolute values of the indices.
    # We could also plot the maximum (max) or minimum (min)
    # of the indices, their respective absolute values (abs_max, abs_min),
    # as well as the normal average (average).
    plot_dict(full_stat_dict=abs_average, n_values=n_values,
              name="_FP{}m{}_AbsAverage".format(fp_total, m))
