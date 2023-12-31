import numpy as np

tg = np.tan(np.arccos(0.7))
l=0.5
CASEDATA = {
    'Vs': 10,   # voltage at node 1 (kV)
    'Vmin': 9.5,    # minimum voltage (kV)
    'alpha': 0.6,   # load factor
    'g': 0.1,   # capital recovery rate
    'cl': 0.1,  # cost of energy loss ($/kWh)
    'delta': 50,    # precentage deviation in load Smax*(1'+/-' delta)
    'branch': np.array([
        # from to length(km) P(kW) Q(kVA)
        [0,  1, l,  0.01, 0.01*tg],
        [1,  2, l,  0.01, 0.01*tg],
        [2,  3, l,  56.0, 56.0*tg],
        [3,  4, l,  0.01, 0.01*tg],
        [4,  5, l,  35.28, 35.28*tg],
        [5,  6, l,  0.01, 0.01*tg],
        [6,  7, l,  35.28, 35.28*tg],
        [7,  8, l,  0.01, 0.01*tg],
        [8,  9, l,  0.01, 0.01*tg],
        [9,  10, l,  56, 56*tg],
        [10,  11, l,  0.01, 0.01*tg],
        [11,  12, l,  0.01, 0.01*tg],
        [12,  13, l,  35.28, 35.28*tg],
        [13,  14, l,  35.28, 35.28*tg],
        [1,  15, l,  35.28, 35.28*tg],
        [2,  16, l,  112, 112*tg],
        [4,  17, l,  56, 56*tg],
        [17,  18, l,  56, 56*tg],
        [18,  19, l,  35.28, 35.28*tg],
        [19,  20, l,  35.28, 35.28*tg],
        [20,  21, l,  35.28, 35.28*tg],
        [18,  22, l,  56, 56*tg],
        [6,  23, l,  35.28, 35.28*tg],
        [7,  24, l,  35.28, 35.28*tg],
        [24,  25, l,  56, 56*tg],
        [25,  26, l,  0.01, 0.01*tg],
        [26,  27, l,  56, 56*tg],

        [25,  28, l,  56, 56*tg],
        [26,  29, l,  56, 56*tg],

        [8,  30, l,  56, 56*tg],
        [30,  31, l,  0.01, 0.01*tg],
        [31,  32, l,  56, 56*tg],
        [31,  33, l,  56, 56*tg],
        [33,  34, l,  56, 56*tg],
        [34,  35, l, 56, 56*tg],
        [34,  36, l,  14, 14*tg],
        [36,  37, l,  0.01, 0.01*tg],
        [37, 38, l,  0.01, 0.01*tg],
        [38,  39, l,  56, 56*tg],
        [37,  40, l,  0.01, 0.01*tg],
        [40,  41, l,  0.01, 0.01*tg],
        [41,  42, l,  56, 56*tg],
        [42,  43, l,  0.01, 0.01*tg],
        [43,  44, l,  35.28, 35.28*tg],
        [40,  45, l,  56, 56*tg],
        [41,  46, l,  0.01, 0.01*tg],
        [46,  47, l,  56, 56*tg],
        [46,  48, l,  35.28, 35.28*tg],
        [43,  49, l,  56, 56*tg],
        [38,  50, l,  14, 14*tg],
        [9,  51, l,  56, 56*tg],
        [40,  52, l,  35.28, 35.28*tg],
        [11,  53, l,  56, 56*tg],
        [53,  54, l,  0.01, 0.01*tg],
        [54,  55, l,  56, 56*tg],
        [54,  56, l,  35.28, 35.28*tg],
        [56,  57, l,  14, 14*tg],
        [12,  58, l,  35.28, 35.28*tg],

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
        [7, 0], #1
        [7, -1], #2
        [7, -2], #3
        [6, -3], #4
        [6, -4], #5
        [6, -5], #6
        [6, -6], #7
        [5, -7], #8
        [4, -8], #9
        [3, -9], #10
        [2, -10], #11
        [2, -11], #12
        [1, -12], #13
        [0, -13], #14
        [0, -14], #15
        [8, -2], #16
        [7, -3], #17
        [7, -5], #18
        [7, -6], #19
        [7, -7], #20
        [7, -8], #21
        [7, -9], #22
        [8, -7], #23
        [6, -7], #24
        [5, -8], #25
        [5, -9], #26
        [5, -10], #27
        [5, -11], #28
        
        [6, -10], #37
        [6, -11], #38
        
        [4, -9], #57
        [4, -10], #58
        [3, -11], #59
        [4, -11], #60
        [4, -12], #61
        [4, -13], #62
        [5, -12], #63
        [5, -13], #64
        [3, -14], #65
        [2,-15], #66
        [5, -14], #67
        [4, -15], #68
        [4,-16], #69
        [4, -17], #70
        [4, -18], #71
        [5, -15], #72
        [5, -16], #73
        [5, -17], #74
        [6, -17], #75
        [5, -18], #76
        [3, -15], #77
        [3, -10], #78
        [6, -15], #79
        [2, -12], #80
        [2, -13], #81
        [1, -14], #82
        [2, -14], #83
        [1, -15], #84
        [1, -13], #85
    ])
}
