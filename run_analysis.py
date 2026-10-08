import numpy as np
import pandas as pd
import os, sys
from pathlib import Path
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

logger = logging.getLogger(__name__)

# Set work directory
path = '/gpfs/work4/0/FWC2/Spatial_dependence/Disease_outbreak/global'

from eca import (
    ts,
    ts2es,
    coincidence_rate,
    MC_sim,
    poisson_significance_test
)

def run_eca_analysis(
    COUNTRY,
    seriesA,
    seriesB,
    locA,
    locB,
    spanA,
    spanB,
    reps=10000,
    delT_all=None,
    tau=0,
    sigtest="wt.surrogate",
    all_events=False
):
    if delT_all is None:
        delT_all = np.arange(0, 13)  # 0 to 12 months

    # Get all unique locations where events were reported
    loc_all = list(dict.fromkeys(
        item
        for sub in locA + locB
        for item in sub
    ))

    res = pd.DataFrame(
        index=range(len(delT_all)),
        columns=[
            "delT",
            "precursor_rate",
            "trigger_rate",
            "precursor_95",
            "precursor_99",
            "trigger_95",
            "trigger_99",
        ],
    )

    res["delT"] = delT_all

    if all_events:
        df_sur = pd.DataFrame(index=range(reps),columns=delT_all) # document all the surrogate event results

    for row_idx, delT in enumerate(delT_all):

        # Calculating aggregated observed rates
        prec_agg = 0
        trig_agg = 0
        N_A = 0
        N_B = 0

        for loc in loc_all:
            idxA = [i for i, sub in enumerate(locA) if loc in sub]
            seriesA_k = seriesA[idxA]

            idxB = [i for i, sub in enumerate(locB) if loc in sub]
            seriesB_k = seriesB[idxB]

            if len(seriesA_k) > 0 and len(seriesB_k) > 0:
                prec_k, trig_k = coincidence_rate(
                    COUNTRY,
                    seriesA_k,
                    seriesB_k,
                    spanA,
                    spanB,
                    [loc] * len(seriesA_k),
                    [loc] * len(seriesB_k),
                    delT,
                    tau=tau,
                )

                prec_k = (0 if pd.isna(prec_k) else prec_k) * len(seriesA_k)
                trig_k = (0 if pd.isna(trig_k) else trig_k) * len(seriesB_k)

                prec_agg += prec_k
                trig_agg += trig_k
                N_A += len(seriesA_k)
                N_B += len(seriesB_k)

        res.loc[row_idx, "precursor_rate"] = prec_agg / N_A if N_A != 0 else 0
        res.loc[row_idx, "trigger_rate"] = trig_agg / N_B if N_B != 0 else 0

        if N_A == 0:
            print(COUNTRY)

        # Calculating significance levels using surrogate events
        surdist = np.full((reps, 2), np.nan)

        for n in range(reps):
            n_prec_agg = 0
            n_trig_agg = 0
            sur_N_A = 0
            sur_N_B = 0

            for loc in loc_all:
                idxA = [i for i, sub in enumerate(locA) if loc in sub]
                seriesA_k = seriesA[idxA]

                idxB = [i for i, sub in enumerate(locB) if loc in sub]
                seriesB_k = seriesB[idxB]

                if len(seriesA_k) > 0 and len(seriesB_k) > 0:
                    n_prec_k, n_trig_k, NA_k, NB_k = MC_sim(
                        COUNTRY,
                        seriesA_k,
                        seriesB_k,
                        spanA,
                        spanB,
                        [loc] * len(seriesA_k),
                        [loc] * len(seriesB_k),
                        delT,
                        tau=tau,
                        sigtest=sigtest,
                        reps=1,
                    )

                    n_prec_agg += n_prec_k
                    n_trig_agg += n_trig_k
                    sur_N_A += NA_k
                    sur_N_B += NB_k

            surdist[n, 0] = n_prec_agg / sur_N_A if sur_N_A != 0 else 0
            surdist[n, 1] = n_trig_agg / sur_N_B if sur_N_B != 0 else 0
            if all_events:
                df_sur.iloc[:,delT] = surdist[:, 0]

        res.loc[row_idx, [
            "precursor_95",
            "precursor_99",
            "trigger_95",
            "trigger_99",
        ]] = [
            np.nanquantile(surdist[:, 0], 0.95),
            np.nanquantile(surdist[:, 0], 0.99),
            np.nanquantile(surdist[:, 1], 0.95),
            np.nanquantile(surdist[:, 1], 0.99),
        ]

    if all_events:
        return res, df_sur
    else:
        return res

def str_to_bool(value):
    return value.strip().lower() in ["true", "1", "yes", "y"]

## Main
args = sys.argv[1:]

for i in range(0, len(args), 2):
    key = args[i]
    value = args[i + 1]

    if key == "--scenario":
        scenario = value

    if key == "--waterborne":
        waterborne = str_to_bool(value)

    if key == "--hazard_grouping":
        hazard_grouping = str_to_bool(value)

    if key == "--natural_all":
        natural_all = str_to_bool(value)

if waterborne:
    sub_folder = 'waterborne'
else:
    sub_folder = 'epidemic'

logger.info("Starting ECA analysis")
logger.info(f"Scenario: {scenario}. "+f"Disease type set to {sub_folder} "+f"Hazard grouping set to {hazard_grouping}")

if hazard_grouping:
    hazard_type_all = ['Climatological', 'Geophysical', 'Hydrological', 'Meteorological']
else:
    hazard_type_all = ['Flood', 'Drought', 'Storm', 'Mass movement','Volcanic activity', 'Wildfire', 'Extreme temperature','Earthquake']

if natural_all:
    hazard_type_all = ['Natural Hazard']

# run simulations
if scenario == "all":
    for hazard_type in hazard_type_all:
        data_folder = Path(os.path.join(path,f'data/eca_binary/{sub_folder}/country',hazard_type))
        country_list = [
            file.stem.split("_")[-1]
            for file in data_folder.glob("df_eca_Epidemic_*.csv")
        ]
        logger.info(f"Hazard type: {hazard_type}, "+"found {:d} countries".format(len(country_list)))
        
        for COUNTRY in country_list:
            logger.info(f"Calculating ECA rates for {hazard_type} for {COUNTRY}")

            df_hazards = pd.read_csv(os.path.join(path,f'data/eca_binary/{sub_folder}/country', 
                                                hazard_type,
                                                "df_eca_{:s}_".format(hazard_type) + COUNTRY + '.csv'))
            df_epidemics = pd.read_csv(os.path.join(path,f'data/eca_binary/{sub_folder}/country',
                                                    hazard_type,
                                                "df_eca_Epidemic_" + COUNTRY + '.csv'))

            # Preparing data for ECA analysis
            date_range = (min(df_hazards["Start date"].min(), df_epidemics["Start date"].min()), 
                        max(df_hazards["Start date"].max(), df_epidemics["Start date"].max()))

            # Epidemics
            event_dates = df_epidemics.loc[df_epidemics["Epidemic"].astype(bool), "Start date"]
            epidemic_binary, locations_epidemic = ts(date_range, event_dates, locations=df_epidemics["Location"], pad=182, freq="monthly")
            seriesA, spanA = ts2es(epidemic_binary)
            locA = [list(dict.fromkeys(loc)) for loc in locations_epidemic if loc]

            # Hazards
            event_dates = df_hazards.loc[df_hazards[hazard_type].astype(bool), "Start date"]
            hazard_binary, locations_hazard = ts(date_range, event_dates, locations=df_hazards["Location"], pad=182, freq="monthly")
            seriesB, spanB = ts2es(hazard_binary)
            locB = [list(dict.fromkeys(loc)) for loc in locations_hazard if loc]

            # Run ECA
            res = run_eca_analysis(
                COUNTRY=COUNTRY,
                seriesA=seriesA,
                seriesB=seriesB,
                locA=locA,
                locB=locB,
                spanA=spanA,
                spanB=spanB,
                reps=1000,
            )

            # Save
            outdir = os.path.join(path,'output', f'{sub_folder}/{scenario}', 'country', hazard_type)
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            res.to_csv(os.path.join(outdir, "ECA_result_" + COUNTRY + '_{:s}.csv'.format(hazard_type)), index=False)

else:

    data_folder = Path(os.path.join(path,f'data/eca_binary/multi_hazard/{sub_folder}',f'{scenario}_hazard'))
    country_list = [
        file.stem.split("_")[-1]
        for file in data_folder.glob(f"df_{scenario}_hazard_*.csv")
    ]
    logger.info(f"Hazard type: {scenario} hazards, "+"found {:d} countries".format(len(country_list)))
    hazard_type = scenario

    for COUNTRY in country_list:
        logger.info(f"Calculating ECA rates for {scenario} hazards for {COUNTRY}")

        df_hazards = pd.read_csv(os.path.join(data_folder, 
                                            f"df_{scenario}_hazard_{COUNTRY}.csv"))
        df_hazards[hazard_type] = True
        df_epidemics = pd.read_csv(os.path.join(path,f'data/eca_binary/multi_hazard/{sub_folder}/Natural Hazard',
                                            "df_eca_Epidemic_" + COUNTRY + '.csv'))

        # Preparing data for ECA analysis
        date_range = (min(df_hazards["Start date"].min(), df_epidemics["Start date"].min()), 
                    max(df_hazards["Start date"].max(), df_epidemics["Start date"].max()))

        # Epidemics
        event_dates = df_epidemics.loc[df_epidemics["Epidemic"].astype(bool), "Start date"]
        epidemic_binary, locations_epidemic = ts(date_range, event_dates, locations=df_epidemics["Location"], pad=182, freq="monthly")
        seriesA, spanA = ts2es(epidemic_binary)
        locA = [list(dict.fromkeys(loc)) for loc in locations_epidemic if loc]

        # Hazards
        event_dates = df_hazards.loc[df_hazards[hazard_type].astype(bool), "Start date"]
        hazard_binary, locations_hazard = ts(date_range, event_dates, locations=df_hazards["Location"], pad=182, freq="monthly")
        seriesB, spanB = ts2es(hazard_binary)
        locB = [list(dict.fromkeys(loc)) for loc in locations_hazard if loc]

        # Run ECA
        res, df_sur = run_eca_analysis(
            COUNTRY=COUNTRY,
            seriesA=seriesA,
            seriesB=seriesB,
            locA=locA,
            locB=locB,
            spanA=spanA,
            spanB=spanB,
            reps=1000,
            all_events=True  # for documenting all surrogate events for multi and single ######## will set to False for normal runs
        )

        # Save
        outdir = os.path.join(path,'output',  f'{sub_folder}/{scenario}')
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        res.to_csv(os.path.join(outdir, "ECA_result_" + COUNTRY + '_{:s}.csv'.format(scenario)), index=False)
        df_sur.to_csv(os.path.join(outdir, "Surrogate_result_" + COUNTRY + '_{:s}.csv'.format(scenario)), index=False)