import numpy as np
from math import ceil
from scipy.special import binom


class BaseComparisons(object):
    """Base class to compare arbitrary numbers of fingerprints.

    Attributes
    ----------
    fingerprints : {np.ndrarray, int}
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
            If an int is given it assumes that one is comparing
            n random fingerprints of infinite length.
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
    """

    def __init__(self, fingerprints, c_threshold=None, w_factor="fraction"):
        """Initialize the object.

        Parameters
        ----------
        fingerprints : {np.ndrarray, int}
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
            If an int is given it assumes that one is comparing
            n random fingerprints of infinite length.
        c_threshold : {None, 'dissimilar', int}
            Coincidence threshold.
        w_factor : {"fraction", "power_n"}
            Type of weight function that will be used.
        """
        self.assign_fingerprints(fingerprints)
        self.assign_c_threshold(c_threshold)
        self.set_matches()
        self.set_d_vector()
        self.set_w_factor(w_factor)
        self.set_weighted_matches()
        self.set_a()
        self.set_d()
        self.set_weighted_a()
        self.set_weighted_d()
        self.set_dis_counters()
        self.set_weighted_dis_counters()
        self.set_total_sim_counter()
        self.set_total_weighted_sim_counter()
        self.total_dis_counters()
        self.total_weighted_dis_counters()
        self.set_p()
        self.set_weighted_p()

    @property
    def n_fingerprints(self):
        """Return number of fingerprints.

        Returns
        -------
        n_fingerprints : int
            Number of fingerprints that will be compared.

        Note: If fingerprints is an int this is taken as the number of fingerprints
              that will be compared.
        """
        if isinstance(self.fingerprints, int):
            return self.fingerprints
        else:
            return len(self.fingerprints)

    def assign_fingerprints(self, fingerprints):
        """Assign fingerprints.

        Parameters
        ----------
        fingerprints : {np.ndrarray, int}
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
            If an int is given it assumes that one is comparing
            n random fingerprints of infinite length.

        Raises
        ------
        TypeError
            If fingerprints is not a numpy array.
            If the elements of fingerprints are not numpy arrays.
        ValueError
            If fingerprints is not a positive integer.
            If less than two fingerprints are provided.
            If not all the fingerprints have the same length.
        """
        if isinstance(fingerprints, int):
            if fingerprints <= 0:
                raise ValueError("If fingerprints is given as an integer,"
                                 "it should be positive integer")
            self.fingerprints = fingerprints
        else:
            if not isinstance(fingerprints, np.ndarray):
                raise TypeError("Fingerprints must be a numpy array or an int.")
            if not all(isinstance(fingerprint, np.ndarray) for fingerprint in fingerprints):
                raise TypeError("The elements of fingerprints must be a numpy array.")
            if len(fingerprints) < 2:
                raise ValueError("A minimum of 2 fingerprints must be provided.")
            if not all([len(fingerprint) == len(fingerprints[0]) for fingerprint in fingerprints]):
                raise ValueError("All the fingerprints must have the same length.")
            self.fingerprints = fingerprints

    def assign_c_threshold(self, c_threshold):
        """Assign coincidence threshold.

        Parameters
        ----------
        c_threshold : {None, 'dissimilar', int}
            Coincidence threshold.
            None : Default, c_threshold = n_fingerprints % 2
            'dissimilar' : c_threshold = ceil(n_fingerprints / 2)
            int : Integer number < n_fingerprints

        Raises
        ------
        TypeError
            If c_threshold is not None, 'dissimilar', or an integer.
        ValueError
            If c_threshold is an integer equal or greater than n_fingerprints
        """
        if not c_threshold:
            self.c_threshold = self.n_fingerprints % 2
        if isinstance(c_threshold, str):
            if c_threshold != 'dissimilar':
                raise TypeError("c_threshold must be None, 'dissimilar', or an integer.")
            else:
                self.c_threshold = ceil(self.n_fingerprints / 2)
        if isinstance(c_threshold, int):
            if c_threshold >= self.n_fingerprints:
                raise ValueError("c_threshold cannot be equal or greater than n_fingerprints.")
            self.c_threshold = c_threshold

    def set_matches(self):
        """Calculate the matches between the fingerprints."""
        if isinstance(self.fingerprints, int):
            matches = [int(binom(self.n_fingerprints, k)) for k in range(self.n_fingerprints + 1)]
        else:
            c_total = np.sum(self.fingerprints, axis=0)
            matches = (self.n_fingerprints + 1) * [0]
            for i in range(self.n_fingerprints + 1):
                matches[i] = np.count_nonzero(c_total == i)
        self.matches = np.array(matches)

    def set_d_vector(self):
        """Calculate the d vector.

        Notes
        -----
        The entries of this vector are the numbers |2k - n_fingerprints|,
        which measure the degree of coincidence between the given fingerprints.
        """
        self.d_vector = np.array([abs(2 * k - self.n_fingerprints) for k in range(
                                 self.n_fingerprints + 1)])

    def set_w_factor(self, w_factor):
        """Calculate weight factors.

        Parameters
        ----------
        w_factor : {"fraction", "power_n"}
            Type of weight function that will be used.
            'fraction' : similarity = d[k]/n
                         dissimilarity = 1 - (d[k] - n_fingerprints % 2)/n_fingerprints
            'power_n' : similarity = n**-(n_fingerprints - d[k])
                        dissimilarity = n**-(d[k] - n_fingerprints % 2)
            other values : similarity = dissimilarity = 1
        """
        if w_factor == "power_n":
            power = int(w_factor.split("_")[-1])

            def f_s(d):
                return power**-(self.n_fingerprints - d)

            def f_d(d):
                return power**-(d - self.n_fingerprints % 2)
        elif w_factor == "fraction":
            def f_s(d):
                return d/self.n_fingerprints

            def f_d(d):
                return 1 - (d - self.n_fingerprints % 2)/self.n_fingerprints
        else:
            def f_s(d):
                return 1

            def f_d(d):
                return 1
        weights = (self.n_fingerprints + 1) * [0]
        for k in range(self.n_fingerprints + 1):
            if self.d_vector[k] > self.c_threshold:
                weights[k] = f_s(self.d_vector[k])
            else:
                weights[k] = f_d(self.d_vector[k])
        self.weights = np.array(weights)

    def set_weighted_matches(self):
        """Calculate weighted matches."""
        self.weighted_matches = self.matches * self.weights

    def set_a(self):
        """Calculate the (unweighted) 1-similarity counter."""
        a = 0
        for k in range(self.n_fingerprints + 1):
            if 2 * k - self.n_fingerprints > self.c_threshold:
                a += self.matches[k]
        self.a = a

    def set_d(self):
        """Calculate the (unweighted) 0-similarity counter."""
        d = 0
        for k in range(self.n_fingerprints + 1):
            if self.n_fingerprints - 2 * k > self.c_threshold:
                d += self.matches[k]
        self.d = d

    def set_weighted_a(self):
        """Calculate the (weighted) 1-similarity counter."""
        w_a = 0
        for k in range(self.n_fingerprints + 1):
            if 2 * k - self.n_fingerprints > self.c_threshold:
                w_a += self.weighted_matches[k]
        self.w_a = w_a

    def set_weighted_d(self):
        """Calculate the (weighted) 0-similarity counter."""
        w_d = 0
        for k in range(self.n_fingerprints + 1):
            if self.n_fingerprints - 2 * k > self.c_threshold:
                w_d += self.weighted_matches[k]
        self.w_d = w_d

    def set_dis_counters(self):
        """Calculate the (unweighted) dissimilarity counters."""
        d_counters = []
        for k in range(self.n_fingerprints + 1):
            if self.d_vector[k] <= self.c_threshold:
                d_counters.append(self.matches[k])
        self.dis_counters = np.array(d_counters)

    def set_weighted_dis_counters(self):
        """Calculate the (weighted) dissimilarity counters."""
        w_d_counters = []
        for k in range(self.n_fingerprints + 1):
            if self.d_vector[k] <= self.c_threshold:
                w_d_counters.append(self.weighted_matches[k])
        self.w_dis_counters = np.array(w_d_counters)

    def set_total_sim_counter(self):
        """Calculate the total number of (unweighted) similarity counters."""
        self.total_sim = self.a + self.d

    def set_total_weighted_sim_counter(self):
        """Calculate the total number of (weighted) similarity counters."""
        self.total_w_sim = self.w_a + self.w_d

    def total_dis_counters(self):
        """Calculate total number of (unweighted) dissimilarity counters."""
        self.total_dis = np.sum(self.dis_counters)

    def total_weighted_dis_counters(self):
        """Calculate total number of (weighted) dissimilarity counters."""
        self.total_w_dis = np.sum(self.w_dis_counters)

    def set_p(self):
        """Calculate p."""
        self.p = self.total_sim + self.total_dis

    def set_weighted_p(self):
        """Calculate weighted p."""
        self.w_p = self.total_w_sim + self.total_w_dis
