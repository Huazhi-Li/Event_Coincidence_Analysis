# -*- coding: utf-8 -*-
"""
Function to calculate precursor and trigger coincidence rates for a pair of event series for individual cases.
Translated from the R package CoinCalc based on Donges et al. (2016).

@author: Huazhi Li (huazhi.li@vu.nl)
"""
import numpy as np

def coincidence_rate(COUNTRY, seriesA, seriesB, spanA, spanB, locA, locB, delT, tau=0):
    """
    COUNTRY:            The analysed country.
    seriesA,seriesB:    Discrete time stamp of eventA and eventB.
    spanA, spanB:       (start, end) time span of seriesA and seriesB.
    locA, locB:         Location (admin1) of Event A and B.
    delT:               Coincidence window.
    tau:                Time lag between event series; can be incorporated to delT, therefore set to 0。
    """
    
    # ---- ERROR #1: negative parameters
    if tau < 0 or delT < 0:
        raise ValueError(
            "ERROR #1: The time lag (tau) or delta T (delT) is negative."
        )
    
    # ---- ERROR #2 & #3: not vector-like
    if not isinstance(seriesA, (list, np.ndarray)):
        raise TypeError(
            "ERROR #2: seriesA is not a vector."
        )
    
    if not isinstance(seriesB, (list, np.ndarray)):
        raise TypeError(
            "ERROR #3: seriesB is not a vector."
        )
    
    # ---- ERROR #4: no overlapping time span
    
    if spanA[0] > spanB[1] or spanB[0] > spanA[1]:
        raise ValueError(
            "ERROR #4: No common time span found for spanA and spanB."
        )
        
    # EVENT COINCIDENCE ANALYSIS
    seriesA = np.asarray(seriesA)
    seriesB = np.asarray(seriesB)

    N_A = len(seriesA)
    N_B = len(seriesB)

    # # common time span
    span = np.arange(
        max(spanA[0], spanB[0]),
        min(spanA[1], spanB[1]) + 1
    )

    # ---- Precursor coincidence rate
    K_prec = 0

    for a, loca in zip(seriesA, locA):
        if a not in span:
            continue

        for b, locb in zip(seriesB, locB):
            if b not in span or b>a:
                continue
            
            dt  = (a - tau) - b # time window

            if 0 <= dt <= delT and (loca == [COUNTRY] or locb == [COUNTRY] or set(loca).intersection(locb)): # check if the location matches
                    K_prec += 1
                    break


    CRprec = K_prec / N_A if N_A > 0 else np.nan

    # ---- Trigger coincidence rate
    K_trigg = 0

    for b, locb in zip(seriesB, locB):
        if b not in span:
            continue
        
        for a, loca in zip(seriesA, locA):
            if a not in span or a>b:
                continue

            dt = (b - tau) - a # time window

            if 0 <= dt <= delT and (loca == [COUNTRY] or locb == [COUNTRY] or set(loca).intersection(locb)): # check if the location matches
                K_trigg += 1
                break


    CRtrigg = K_trigg / N_B if N_B > 0 else np.nan

    return CRprec, CRtrigg
