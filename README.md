# Event_Coincidence_Analysis
This package provides a python version of the Event Coincidence Analysis, based on the R package CoinCalc (Donges et al., 2016). The generated scripts are used to assess the statistical correlations (i.e. precursor and trigger rates) between natural hazards and the subsequent (waterborne) disease outbreaks. For details, see the manuscript de Ruiter, M. C.\*, **Li, H.**\*, & Jäger, W. S. (2026). Global hotspots and temporal dynamics of post-hazard waterborne disease outbreaks. *To be submitted*. A preprint doi will be provided soon.

**Input dataset**

Historic natural hazard and disease outbreak records during 1950 and 2024 are retrieved from EM-DAT. An addtional EM-DAT dataset is used to provide the administrative 1 information of these events from 2000 onwards. These data can be accessed via https://www.emdat.be/.

**Scripts**

*preprocess*

combining_emdat_admin1_datasets.py-->combining the EM-DAT (1950-2024) with the additonal dataset containing detailed administrative 1 information for events from 2000 onwards
hazard_epidemic_extracation.py-->extracting natural hazard and disease outbreak records per country from the combined EM-DAT dataset
multi_hazard_identification.py-->identifying single-hazard and multi-hazard events per country

*eca*

coincidence_rate.py-->calculating the precursor and triggering rates
event_series.py-->converting discrete event records to binary event time series at a monthly time step
significant_test.py-->testing the significance of calculated rate, including 1) poisson significance test; and 2) surrogate event significance test
run_analysis.py-->running the event coincidence analysis

*postprocess*

agg_global_results-->showing global results at a country scale for each natural hazard - disease outbreak combination

*plot (for producing figures in the manuscript)*

global_lag_figure.py-->plotting global aggregated precursor rates for each natural hazard - disease outbreak combination
identify_most_sig_lag.py-->identifying the first significant lag and plotting a global figure
multi_versus_single.py-->plotting the relative difference in single and multi precursor rates per country

**References**

Donges, J. F., Schleussner, C. F., Siegmund, J. F., & Donner, R. V. (2016). Event coincidence analysis for quantifying statistical interrelationships between event time series: on the role of flood events as triggers of epidemic outbreaks. The European Physical Journal Special Topics, 225(3), 471-487.
