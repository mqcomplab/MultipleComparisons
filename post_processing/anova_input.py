import random
import re


def gen_headers(cv_file):
    """Generate the headers that will be used in the ANOVA.

    Arguments
    ---------
    cv_file : str
        File containing the results of the resumed CV analysis.

    Raises
    ------
    TypeError
        If cv_vile is not a resumed CV file.

    Returns
    -------
    headers : list
        List with the headers that will be used.
    No: number of the data entry.
    F1: one of 'Norm': normalized srd value, 'A': sequential CV, 'B': random CV.
    F2: nN, with N being the number of fingerprints compared simultaneously.
    F3: one of 'w': weighted index, 'nw': unweighted index.
    The rest of the headers are the similarity indices.
    """
    if not re.match(r'SRD.*FP.*m.*n.*CV.*txt', cv_file):
        raise TypeError("cv_file must be a CV file.")
    with open(cv_file, "r") as infile:
        lines = infile.readlines()
    headers = lines[0].strip().split()[1:]
    headers = ["No", "F1", "F2", "F3"] + headers
    return headers


def read_data(cv_file):
    """Read the data in the resumed cv file.

    Arguments
    ---------
    cv_file : str
        File containing the results of the resumed CV analysis.

    Raises
    ------
    TypeError
        If cv_vile is not a resumed CV file.

    Returns
    -------
    content, na, nb
    content: Raw content of the cv analysis.
    na: Line where the sequential CV results end.
    nb: Line where the random CV results end.
    """
    if not re.match(r'SRD.*FP.*m.*n.*CV.*txt', cv_file):
        raise TypeError("cv_file must be a CV file.")
    with open(cv_file, "r") as infile:
        lines = infile.readlines()
    content = lines[2:]
    na = 0
    nb = 0
    for line in content:
        if "_A" in line:
            na += 1
        elif "_B" in line:
            nb += 1
    infile.close()
    return content, na, nb


def extract_cv_name_data(cv_file):
    """Extract the calculation information from the name of the cv_file.

    Arguments
    ---------
    cv_file : str
        File containing the results of the resumed CV analysis.

    Raises
    ------
    TypeError
        If cv_vile is not a resumed CV file.

    Returns
    -------
    weight_info, n, fp_total, m, CV
    weight_info: 'w' or 'nw'.
    n: Number of fingerprints compared simultaneously.
    fp_total: Total number of fingerprints.
    m: Length of the fingerprints.
    CV: Fraction of the data left out as test set.
    """
    if not re.match(r'SRD.*FP.*m.*n.*CV.*txt', cv_file):
        raise TypeError("cv_file must be a CV file.")
    data = cv_file.split("SRD")[-1]
    weight_info = data.split("FP")[0]
    n = data.split("n")[-1].split("CV")[0]
    fp_total = data.split("FP")[-1].split("m")[0]
    m = data.split("m")[-1].split("n")[0]
    CV = data.split("CV")[-1].split("_")[0]
    return weight_info, n, fp_total, m, CV


def gen_b_indices(na, nb):
    """Select the indices of the random CV results that will be used in the ANOVA.

    Arguments
    ---------
    na : int
        Line where the sequential CV results end.
    nb : int
        Line where the random CV results end.

    Returns
    -------
    List of indices of the random CV results that will be used in the ANOVA.
    """
    return random.sample(range(na + 1, nb + 1), na)


def process_cv_file(cv_file, b_indices):
    """Extract the CV results that will be used in the ANOVA from a single cv file.

    Arguments
    ---------
    cv_file : str
        File containing the results of the resumed CV analysis.
    b_indices : list
        List of indices of the random CV results that will be used in the ANOVA.

    Raises
    ------
    TypeError
        If cv_vile is not a resumed CV file.

    Returns
    -------
    file_lists : list
        List that contains two sub-lists. The first has the results of the sequential CV.
        The second has the results of the selected random CV.
    """
    if not re.match(r'SRD.*FP.*m.*n.*CV.*txt', cv_file):
        raise TypeError("cv_file must be a CV file.")
    content, na, nb = read_data(cv_file)
    weight_info, n, FP, m, CV = extract_cv_name_data(cv_file)
    norm_line = content[0]
    a_lines = content[1: na + 1]
    b_lines = []
    for index in b_indices:
        b_lines.append(content[index])
    file_lists = []
    norm_line = ["Norm", "n" + n, weight_info] + norm_line.strip().split()[1:]
    file_lists.append(norm_line)
    for a_line in a_lines:
        a_line = ["A", "n" + n, weight_info] + a_line.strip().split()[1:]
        file_lists.append(a_line)
    for b_line in b_lines:
        b_line = ["B", "n" + n, weight_info] + b_line.strip().split()[1:]
        file_lists.append(b_line)
    return file_lists


def gen_total_data(cv_file_list, b_indices):
    """Extract the CV results that will be used in the ANOVA from a list of cv files.

        Arguments
        ---------
        cv_file_list : list
            List with the files containing the results of the resumed CV analysis.
        b_indices : list
            List of indices of the random CV results that will be used in the ANOVA.

        Returns
        -------
        total_data : list
            List that contains all the CV results that will be used in the ANOVA.
        """
    total_data = []
    cv_files = {"w": {}, "nw": {}}
    for cv_file in cv_file_list:
        weight_info, n, FP, m, CV = extract_cv_name_data(cv_file)
        cv_files[weight_info][int(n)] = cv_file
    
    for w in sorted(cv_files):
        for n in sorted(cv_files[w]):
            cv_file = cv_files[w][n]
            total_data += process_cv_file(cv_file, b_indices)
    total_data = [[str(i), str(i)] + line for i, line in enumerate(total_data, 1)]
    return total_data


def gen_outfile_name(cv_file_list):
    """Generate the name of the input file for the ANOVA.

    Arguments
    ---------
    cv_file_list : list
        List with the files containing the results of the resumed CV analysis.

    Returns
    -------
    s : str
        Name of the file that will serve as input to the ANOVA.
    """
    weights = []
    n_min = 2
    n_max = 2
    for cv_file in cv_file_list:
        weight_info, n, FP, m, CV = extract_cv_name_data(cv_file)
        if weight_info not in weights:
            weights.append(weight_info)
        if int(n) < n_min:
            n_min = int(n)
        elif int(n) > n_max:
            n_max = int(n)
    s = "ARSRD" + sorted(weights)[0] + sorted(weights)[1] + "FP" + FP + "m" + m + "n" + str(n_min)
    s += "n" + str(n_max) + "CV" + CV + ".dat"
    return s
    

def output_file(headers, total_data, outfile_name):
    """Generate the file that will be used as input in the ANOVA.

    Arguments
    ---------
    headers : list
        List with the headers that will be used.
    total_data : list
        List that contains all the CV results that will be used in the ANOVA.
    outfile_name : str
        Name of the file that will serve as input to the ANOVA.

    Creates the .dat file that will be used as input in the ANOVA.
    """
    s = "              "
    for header in headers:
        s += "{:14}".format(header)
    s += "\n"
    for line in total_data:
        for item in line:
            s += "{:14}".format(item)
        s += "\n"
    with open(outfile_name, "w") as outfile:
        outfile.write(s)


if __name__ == "__main__":
    # Generates the ANOVA input .dat file for all the resumed CV files in the given folder.
    import glob
    cv_file_list = [cv_file for cv_file in glob.glob("*.txt")
                    if re.match(r'SRD.*FP.*m.*n.*CV.*txt', cv_file)]
    headers = gen_headers(cv_file_list[0])
    _, na, nb = read_data(cv_file_list[0])
    b_indices = gen_b_indices(na, nb)
    total_data = gen_total_data(cv_file_list, b_indices)
    outfile_name = gen_outfile_name(cv_file_list)
    output_file(headers, total_data, outfile_name)
