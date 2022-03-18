from indices.base import BaseComparisons


class RussellRao(BaseComparisons):
    """Class to calculate the Russell-Rao index.

    n=2 formula:
        a/p

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
    rr_1sim_wdis()
        Calculate the index with 1-sim-counters and with weighted denominator.
    rr_1sim_dis()
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
        self.rr_sim_wdis()
        self.rr_1sim_wdis()
        self.rr_sim_dis()
        self.rr_1sim_dis()

    def rr_1sim_wdis(self):
        """Calculate the index with 1-sim-counters and with weighted denominator.

        Note
        ----
        w_a/w_p
        """
        numerator = self.w_a
        denominator = self.w_p
        self.RR_1sim_wdis = numerator/denominator

    def rr_1sim_dis(self):
        """Calculate the index with 1-sim-counters and with unweighted denominator.

        Note
        ----
        w_a/p
        """
        numerator = self.w_a
        denominator = self.p
        self.RR_1sim_dis = numerator/denominator
