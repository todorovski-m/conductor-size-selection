import numpy as np

tg = np.tan(np.arccos(0.9))
CASEDATA = {
    'Vs': 10,   # voltage at node 1 (kV)
    'Vmin': 9.5,    # minimum voltage (kV)
    'alpha': 0.6,   # load factor
    'g': 0.1,   # capital recovery rate
    'cl': 0.1,  # cost of energy loss ($/kWh)
    'delta': 50,    # precentage deviation in load Smax*(1'+/-' delta)
    'branch': np.array([
        # from to length(km) P(kW) Q(kVA)
        [0,  1, 1.2,  225.0, 225.0*tg],
        [1,  2, 1.05, 144.0, 144.0*tg],
        [2,  3, 0.75, 90.0, 90.0*tg],
        [1,  4, 1.00, 90.0, 90.0*tg],
        [1,  5, 0.75, 45.0, 45.0*tg],
        [2,  6, 1.25, 90.0, 90.0*tg],
        [3,  7, 1.65, 90.0, 90.0*tg],
        [5,  8, 1.75, 225.0, 225.0*tg],
        [5,  9, 2.00, 144.0, 144.0*tg],
        [6, 10, 1.75, 90.0, 90.0*tg],
        [3, 11, 1.50, 144.0, 144.0*tg],
        [10, 12, 0.45, 90.0, 90.0*tg],
        [9, 13, 0.50, 90.0, 90.0*tg],
        [8, 14, 2.20, 90.0, 90.0*tg],
        [1, 15, 0.20, 135.0, 135.0*tg],
        [15, 16, 1.00, 72.0, 72.0*tg],
        [16, 17, 1.05, 36.0, 36.0*tg],
        [3, 18, 0.50, 90.0, 90.0*tg],
        [18, 19, 1.75, 36.0, 36.0*tg],
        [11, 20, 1.55, 54.0, 54.0*tg],
        [20, 21, 0.75, 36.0, 36.0*tg],
        [7, 22, 0.75, 72.0, 72.0*tg],
        [3, 23, 0.50, 90.0, 90.0*tg],
        [8, 24, 0.50, 27.0, 27.0*tg],
    ]),

    'cable': np.array([
        # area(mm^2) r(ohm/km) x(ohm/km) Imax(A) instalation($/km)
        [16,   2.003,  0.136,   85,  12200],
        [25,   1.282,  0.131,  130,  16000],
        [35,   0.866,  0.126,  155,  17200],
        [50,   0.641,  0.121,  180,  20300],
        [70,   0.443,  0.116,  225,  23300],
        [95,   0.320,  0.112,  270,  27600],
        [120,  0.253,  0.108,  305,  34500],
        [150,  0.206,  0.105,  340,  39000],
        [185,  0.130,  0.076,  380,  43500],
    ]),

    'xy': np.array([
        # x(mm) y(mm)
        [3.26, 0.59],
        [4.56, 1.71],
        [3.51, 2.53],
        [2.32, 3.01],
        [6.92, 2.33],
        [5.56, 3.35],
        [3.57, 4.63],
        [1.70, 4.02],
        [7.09, 3.82],
        [5.77, 5.36],
        [1.48, 6.01],
        [1.01, 3.38],
        [3.13, 7.05],
        [7.78, 7.05],
        [8.45, 2.94],
        [5.12, 0.9],
        [5.87, 1.53],
        [7.64, 1.68],
        [1.70, 1.80],
        [1.36, 2.56],
        [0.53, 3.89],
        [1.06, 4.31],
        [1.50, 5.17],
        [3.09, 3.51],
        [7.42, 4.38],
    ])
}