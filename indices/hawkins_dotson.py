from indices.base import BaseComparisons


class HawkinsDotson(BaseComparisons):
    """Class to calculate the Hawkins-Dotson index.

    n=2 formula:
        0.5 * (a/(a + b + c) + d/(d + b + c))

    Attributes
    ----------
    fingerprints : np.ndarray
        Numpy array with the fingerprints that will be compared.
        The fingerprints must be also given as Numpy arrays.
    c_threshold : {None, 'dissimilar', int}
        Coincidence threshold.

    Properties
    ----------
    n_fingerprints : int
        Number of fingerprints that will be compared.

    Methods
    -------
    __init__(self, fingerprints, c_threshold=None, w_factor="fraction")
        Initialize the object.
    assign_fingerprints(fingerprints)
        Assign fingerprints.
    assign_c_threshold(c_threshold)
        Assign coincidence threshold.
    matches()
        Calculate the matches between the fingerprints.
    set_d_vector()
        Calculate the d vector.
    set_w_factor(w_factor)
        Calculate weight factors.
    set_weighted_matches()
        Calculate weighted matches.
    set_a()
        Calculate the (unweighted) 1-similarity counter.
    set_d()
        Calculate the (unweighted) 0-similarity counter.
    set_weighted_a()
        Calculate the (weighted) 1-similarity counter.
    set_weighted_d()
        Calculate the (weighted) 0-similarity counter.
    set_dis_counters()
        Calculate the (unweighted) dissimilarity counters.
    set_weighted_dis_counters()
        Calculate the (weighted) dissimilarity counters.
    set_total_sim_counter()
        Calculate the total number of (unweighted) similarity counters.
    set_total_weighted_sim_counter()
        Calculate the total number of (unweighted) similarity counters.
    total_dis_counters()
        Calculate total number of (unweighted) dissimilarity counters.
    total_weighted_dis_counters()
        Calculate total number of (weighted) dissimilarity counters.
    set_p()
        Calculate p.
    set_weighted_p()
        Calculate weighted p.
    hd_1sim_wdis()
        Calculate the index with 1-sim-counters and with weighted denominator.
    hd_1sim_dis()
         Calculate the index with 1-sim-counters and with unweighted denominator.
    """

    def __init__(self, fingerprints, c_threshold=None, w_factor="fraction"):
        """Initialize the object.

        Parameters
        ----------
        fingerprints : np.ndrarray
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
        c_threshold : {None, 'dissimilar', int}
            Coincidence threshold.
        w_factor : {"fraction", "power_n"}
            Type of weight function that will be used.
        """
        super().__init__(fingerprints, c_threshold, w_factor)
        self.hd_1sim_wdis()
        self.hd_1sim_dis()

    def hd_1sim_wdis(self):
        """Calculate the index with 1-sim-counters and with weighted denominator.

        Note
        ----
        0.5 * (w_a/(w_a + w_b + w_c) + w_d/(w_d + w_b + w_c))
        """
        numerator1 = self.w_a
        numerator2 = self.w_d
        denominator1 = self.w_a + self.total_w_dis
        denominator2 = self.w_d + self.total_w_dis
        self.HD_1sim_wdis = 0.5 * (numerator1/denominator1 + numerator2/denominator2)

    def hd_1sim_dis(self):
        """Calculate the index with 1-sim-counters and with unweighted denominator.

        Note
        ----
        0.5 * (w_a/(a + b + c) + w_d/(d + b + c))
        """
        numerator1 = self.w_a
        numerator2 = self.w_d
        denominator1 = self.a + self.total_dis
        denominator2 = self.d + self.total_dis
        self.HD_1sim_dis = 0.5 * (numerator1/denominator1 + numerator2/denominator2)
