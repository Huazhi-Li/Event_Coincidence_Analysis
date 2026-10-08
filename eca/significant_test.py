# -*- coding: utf-8 -*-
"""
Function to test the significance level of the two event series.
Translated from the R package CoinCalc based on Donges et al. (2016).

@author: Huazhi Li (huazhi.li@vu.nl)
"""

import numpy as np
from math import comb
from .coincidence_rate import coincidence_rate

def poisson_significance_test(CRprec, CRtrigg, N_A, N_B, delT, Tlen, alpha=0.05, tau=0):
    """
    CRprec:     Precursor coincidence rate
    CRtrigg:    Trigger coincidence rate
    N_A:        Number of events in series A
    N_B:        Number of events in series B
    delT:       Tolerance window
    Tlen:       Length of observation span
    tau :       Time lag
    """
   
    if N_A==1 or N_B==1:
        raise ValueError(
            "ERROR: The poisson test is only valid when the number of events is way larger than 1."
        )
    
    # Precursor coincidence probability
    K_prec = int(CRprec*N_A)
    K_trigg = int(CRtrigg*N_B)
    Pprec = 0.0
    base_prob_prec = (1 - (delT / (Tlen - tau))) ** N_B
    success_prob_prec = 1 - base_prob_prec

    for Ktmp in range(K_prec, N_A + 1):
        Ptmp = (
            comb(N_A, Ktmp)
            * (success_prob_prec ** Ktmp)
            * (base_prob_prec ** (N_A - Ktmp))
        )
        Pprec += Ptmp

    # Trigger coincidence probability
    Ptrigg = 0.0
    base_prob_trig = (1 - (delT / (Tlen - tau))) ** N_A
    success_prob_trig = 1 - base_prob_trig

    for Ktmp in range(K_trigg, N_B + 1):
        Ptmp = (
            comb(N_B, Ktmp)
            * (success_prob_trig ** Ktmp)
            * (base_prob_trig ** (N_B - Ktmp))
        )
        Ptrigg += Ptmp

    return {
        "p_value_precursor": Pprec,
        "p_value_trigger": Ptrigg
    }

def MC_sim(COUNTRY, seriesA, seriesB, spanA, spanB, locA, locB, delT, tau=0, sigtest="wt.surrogate", reps=1000):
    """
    COUNTRY:            The analysed country.    
    seriesA,seriesB:    Discrete time stamp of eventA and eventB.
    spanA, spanB:       (start, end) time span of seriesA and seriesB.
    locA, locB:         Location (admin1) of Event A and B.
    delT:               Coincidence window.
    tau:                Time lag between event series; can be incorporated to delT, therefore set to 0
    sigtest:            Surrogate event generation method.
    reps:               Repetition runs.
    """

    N_A = len(seriesA)
    N_B = len(seriesB)

    span = np.arange(max(spanA[0], spanB[0]),
                     min(spanA[1], spanB[1]) + 1)
    
    Tlen = len(span)
    surdist = np.full((reps,2), np.nan)

    # ---- waiting time surrogate ----
    # create 'reps' (standard = 1000) surrogate sime series for seriesA and seriesB, having the same average waiting times between two events. Then performe
    # coincidence analysis as given above for all 1000 pairs of event series.
    
    if sigtest == "wt.surrogate":

        def waiting_times(series, span_start):
            if len(series) <= 1:
                return np.array([np.nan])
            s = np.sort(series)
            wt = np.zeros(len(series))
            wt[0] = s[0] - span_start
            for i in range(1, len(series)):
                wt[i] = s[i] - s[i-1]
            return wt
        
        # both event A and B have more than one case
        if (N_A>1 and N_B>1):
            wtA = waiting_times(seriesA, spanA[0]).astype(int)
            wtB = waiting_times(seriesB, spanB[0]).astype(int)
    
            for surno in range(reps):
    
                # surrogate A
                surA = [np.random.choice(wtA)]
                for i in range(1, N_A+1):
                    tmp = surA[i-1] + np.random.choice(wtA)
                    if tmp > span[-1]:
                        break
                    surA.append(tmp)
                surA = np.array(surA)
    
                # surrogate B
                surB = [np.random.choice(wtB)]
                for i in range(1, N_B+1):
                    tmp = surB[i-1] + np.random.choice(wtB)
                    if tmp > span[-1]:
                        break
                    surB.append(tmp)
                surB = np.array(surB)
    
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
        
        # event A has more than 1 but envet B has less than 1 case, surrogate A uses the waiting time method while surrogate B uses the shuflle method
        if (N_A>1 and N_B<=1):
            wtA = waiting_times(seriesA, spanA[0]).astype(int)
            
            for surno in range(reps):
    
                # surrogate A
                surA = [np.random.choice(wtA)]
                for i in range(1, N_A+1):
                    tmp = surA[i-1] + np.random.choice(wtA)
                    if tmp > span[-1]:
                        break
                    surA.append(tmp)
                surA = np.array(surA)
    
                # surrogate B
                surB = np.random.choice(span, size=N_B, replace=False)
    
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
            
        # event A has less than 1 but envet B has more than 1 case, surrogate A uses the shuflle method while surrogate B uses the waiting time method
        if (N_A<=1 and N_B>1):
            wtB = waiting_times(seriesB, spanB[0]).astype(int)
            
            for surno in range(reps):
    
                # surrogate A
                surA = np.random.choice(span, size=N_A, replace=False)
                
                # surrogate B
                surB = [np.random.choice(wtB)]
                for i in range(1, N_B+1):
                    tmp = surB[i-1] + np.random.choice(wtB)
                    if tmp > span[-1]:
                        break
                    surB.append(tmp)
                surB = np.array(surB)
    
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
                    
        # event A and B both have less than 1 but envet B has more than 1 case, use the shuflle method 
        if (N_A<=1 and N_B<=1):   
            for surno in range(reps):

                surA = np.random.choice(span, size=N_A, replace=False)
                surB = np.random.choice(span, size=N_B, replace=False)
                
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
            
    # ---- shuffle surrogate ----
    # create 'reps' (standard = 1000) shuffled time series for seriesA and seriesB, having the same number of events. Then performe
    # coincidence analysis as given above for all 1000 pairs of event series. 
    if sigtest == "shuffle.surrogate":

        span = np.arange(1, Tlen + 1)

        for surno in range(reps):

            surA = np.random.choice(span, size=N_A, replace=False)
            surB = np.random.choice(span, size=N_B, replace=False)
            
            # precursor and trigger coincidence rates
            surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
            
    # # ---- mix surrogate ----
    # # in some cases, only one series such as floods needs to capture seasonality, i.e. having the averaged waiting times; the other one shall have use the shuffle 
    # # create 'reps' (standard = 1000) time series for seriesA and seriesB, with seriesA having the same number of events and B having the same waiting times. Then perform
    # # coincidence analysis as given above for all 1000 pairs of event series.
    # if sigtest == "mixed":

    #     span = np.arange(1, Tlen + 1)

    #     # waiting times for B
    #     if N_B > 1:
    #         sB = np.sort(seriesB)
    #         wtB = np.zeros(N_B)
    #         wtB[0] = sB[0] - spanB[0]
    #         for i in range(1, N_B):
    #             wtB[i] = sB[i] - sB[i - 1]
    #     else:
    #         wtB = np.array([np.nan])
        
    #     for surno in range(reps):
            
    #         # shuffle samples for A
    #         surA = np.random.choice(span, size=N_A, replace=False)
    #         surB = [np.random.choice(wtB)]
            
    #         for i in range(1, N_B+1):
    #             tmp = surB[i-1] + np.random.choice(wtB)
    #             if tmp > span[-1]:
    #                 break
    #             surB.append(tmp)
    #         surB = np.array(surB)
    #         surB = np.array(surB)
            
    #         # precursor and trigger coincidence rates
    #         surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
        

    # ---- significance level ----
    # the percentile of the empirically found coincidence rate using the ecdf of the 1000 surrogate coincidence rates.
    if reps == 1:
        return surdist[0,0]*len(surA), surdist[0,1]*len(surB),len(surA),len(surB)
    else: 
        return np.nanquantile(surdist[:,0], 0.95), np.nanquantile(surdist[:,0], 0.99), np.nanquantile(surdist[:,1], 0.95), np.nanquantile(surdist[:,1], 0.99)

def MC_sim_agg(COUNTRY, G, seriesA, seriesB, spanA, spanB, locA, locB, delT, tau=0, sigtest="wt.surrogate", reps=1000):
    """
    COUNTRY:            The analysed country.    
    seriesA,seriesB:    Discrete time stamp of eventA and eventB.
    spanA, spanB:       (start, end) time span of seriesA and seriesB.
    locA, locB:         Location (admin1) of Event A and B.
    delT:               Coincidence window.
    tau:                Time lag between event series; can be incorporated to delT, therefore set to 0
    sigtest:            Surrogate event generation method.
    reps:               Repetition runs.
    """
    
    for surno in range(reps):
        
        
        N_A = len(seriesA)
        N_B = len(seriesB)
    
        span = np.arange(max(spanA[0], spanB[0]),
                         min(spanA[1], spanB[1]) + 1)
        
        Tlen = len(span)
        surdist = np.full((reps,2), np.nan)
    
        # ---- waiting time surrogate ----
        # create 'reps' (standard = 1000) surrogate sime series for seriesA and seriesB, having the same average waiting times between two events. Then performe
        # coincidence analysis as given above for all 1000 pairs of event series.
        
        if sigtest == "wt.surrogate":
    
            def waiting_times(series, span_start):
                if len(series) <= 1:
                    return np.array([np.nan])
                s = np.sort(series)
                wt = np.zeros(len(series))
                wt[0] = s[0] - span_start
                for i in range(1, len(series)):
                    wt[i] = s[i] - s[i-1]
                return wt
    
            wtA = waiting_times(seriesA, spanA[0]).astype(int)
            wtB = waiting_times(seriesB, spanB[0]).astype(int)
    
            for surno in range(reps):
    
                # surrogate A
                surA = [np.random.choice(wtA)]
                for i in range(1, N_A+1):
                    tmp = surA[i-1] + np.random.choice(wtA)
                    if tmp > span[-1]:
                        break
                    surA.append(tmp)
                surA = np.array(surA)
    
                # surrogate B
                surB = [np.random.choice(wtB)]
                for i in range(1, N_B+1):
                    tmp = surB[i-1] + np.random.choice(wtB)
                    if tmp > span[-1]:
                        break
                    surB.append(tmp)
                surB = np.array(surB)
    
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
         
        # ---- shuffle surrogate ----
        # create 'reps' (standard = 1000) shuffled time series for seriesA and seriesB, having the same number of events. Then performe
        # coincidence analysis as given above for all 1000 pairs of event series. 
        if sigtest == "shuffle.surrogate":
    
            span = np.arange(1, Tlen + 1)
    
            for surno in range(reps):
    
                surA = np.random.choice(span, size=N_A, replace=False)
                surB = np.random.choice(span, size=N_B, replace=False)
                
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
                
        # ---- mix surrogate ----
        # in some cases, only one series such as floods needs to capture seasonality, i.e. having the averaged waiting times; the other one shall have use the shuffle 
        # create 'reps' (standard = 1000) time series for seriesA and seriesB, with seriesA having the same number of events and B having the same waiting times. Then perform
        # coincidence analysis as given above for all 1000 pairs of event series.
        if sigtest == "mixed":
    
            span = np.arange(1, Tlen + 1)
    
            # waiting times for B
            if N_B > 1:
                sB = np.sort(seriesB)
                wtB = np.zeros(N_B)
                wtB[0] = sB[0] - spanB[0]
                for i in range(1, N_B):
                    wtB[i] = sB[i] - sB[i - 1]
            else:
                wtB = np.array([np.nan])
            
            for surno in range(reps):
                
                # shuffle samples for A
                surA = np.random.choice(span, size=N_A, replace=False)
                surB = [np.random.choice(wtB)]
                
                for i in range(1, N_B+1):
                    tmp = surB[i-1] + np.random.choice(wtB)
                    if tmp > span[-1]:
                        break
                    surB.append(tmp)
                surB = np.array(surB)
                surB = np.array(surB)
                
                # precursor and trigger coincidence rates
                surdist[surno,:] = coincidence_rate(COUNTRY, surA, surB, spanA, spanB, locA, locB, delT, tau=0)
        

    # ---- significance level ----
    # the percentile of the empirically found coincidence rate using the ecdf of the 1000 surrogate coincidence rates.
    return np.nanquantile(surdist[:,0], 0.95), np.nanquantile(surdist[:,0])