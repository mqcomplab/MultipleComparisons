"""Contains miscellaneous information about the implemented indices
[index-name] : (ClassName, Abbreviation, (ways of calculating the index,))
"""

Indices = {"Austin-Colwell": ("AustinColwell", "AC", ("1sim_wdis", "1sim_dis")),
           "Baroni-Urbani-Buser": ("BaroniUrbaniBuser", "BUB", ("1sim_wdis", "1sim_dis")),
           "Consoni-Todeschini1": ("ConsoniTodeschini1", "CT1", ("1sim_wdis", "1sim_dis")),
           "Consoni-Todeschini2": ("ConsoniTodeschini2", "CT2", ("1sim_wdis", "1sim_dis")),
           "Consoni-Todeschini3": ("ConsoniTodeschini3", "CT3", ("1sim_wdis", "1sim_dis")),
           "Consoni-Todeschini4": ("ConsoniTodeschini4", "CT4", ("1sim_wdis", "1sim_dis")),
           "Faith": ("Faith", "Fai", ("1sim_wdis", "1sim_dis")),
           "Gleason": ("Gleason", "Gle", ("1sim_wdis", "1sim_dis")),
           "Goodman-Kruskal": ("GoodmanKruskal", "GK", ("1sim_wdis", "1sim_dis")),
           "Hawkins-Dotson": ("HawkinsDotson", "HD", ("1sim_wdis", "1sim_dis")),
           "Jaccard": ("Jaccard", "Ja", ("1sim_wdis", "1sim_dis", "sim_wdis", "sim_dis")),
           "Jaccard-Tanimoto": ("JaccardTanimoto", "JT", ("1sim_wdis", "1sim_dis")),
           "Rogers-Tanimoto": ("RogersTanimoto", "RT", ("1sim_wdis", "1sim_dis")),
           "Rogot-Goldberd": ("RogotGoldberd", "RG", ("1sim_wdis", "1sim_dis")),
           "Russell-Rao": ("RussellRao", "RR", ("1sim_wdis", "1sim_dis")),
           "Sokal-Michner": ("SokalMichner", "SM", ("1sim_wdis", "1sim_dis")),
           "Sokal-Sneath1": ("SokalSneath1", "SS1", ("1sim_wdis", "1sim_dis")),
           "Sokal-Sneath2": ("SokalSneath2", "SS2", ("1sim_wdis", "1sim_dis"))}
