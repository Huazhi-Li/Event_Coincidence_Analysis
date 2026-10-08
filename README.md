# Event_Coincidence_Analysis
This package provides a python version of the Event Coincidence Analysis, based on the R package CoinCalc (Donges et al., 2016). The generated scripts are used to assess the statistical correlations (i.e. precursor and trigger rates) between natural hazards and the subsequent (waterborne) disease outbreaks. For details, see the manuscript de Ruiter, M. C.\*, **Li, H.**\*, & Jäger, W. S. (2026). Global hotspots and temporal dynamics of post-hazard waterborne disease outbreaks. *To be submitted*. A preprint doi will be provided soon.

## **Input dataset**

Historic natural hazard and disease outbreak records during 1950 and 2024 are retrieved from EM-DAT. An addtional EM-DAT dataset is used to provide the administrative 1 information of these events from 2000 onwards. These data can be accessed via https://www.emdat.be/.

## **Scripts**

### Preprocess

- `combining_emdat_admin1_datasets.py`  
  Combining the EM-DAT dataset (1950–2024) with an additional dataset containing detailed administrative level-1 information for events from 2000 onwards.

- `hazard_epidemic_extracation.py`  
  Extracting natural hazard and disease outbreak records for each country from the combined EM-DAT dataset.

- `multi_hazard_identification.py`  
  Identifying single-hazard and multi-hazard events for each country.

### ECA (Event Coincidence Analysis)

- `coincidence_rate.py`  
  Calculating precursor and trigger coincidence rates.

- `event_series.py`  
  Converting discrete event records into binary event time series at a monthly resolution.

- `significant_test.py`  
  Testing the statistical significance of the calculated coincidence rates using (1) a Poisson significance test and (2) a surrogate event significance test.

- `run_analysis.py`  
  Running the event coincidence analysis.

### Postprocess

- `agg_global_results`  
  Aggregating and presenting global results at the country level for each natural hazard–disease outbreak combination.

### Plot (Figures in the manuscript)

- `global_lag_figure.py`  
  Plotting globally aggregated precursor coincidence rates for each natural hazard–disease outbreak combination.

- `identify_most_sig_lag.py`  
  Identifying the first statistically significant lag and visualizing the results on a global map.

- `multi_versus_single.py`  
  Plotting the relative differences in precursor coincidence rates between single-hazard and multi-hazard events for each country.

## **References**

Donges, J. F., Schleussner, C. F., Siegmund, J. F., & Donner, R. V. (2016). Event coincidence analysis for quantifying statistical interrelationships between event time series: on the role of flood events as triggers of epidemic outbreaks. The European Physical Journal Special Topics, 225(3), 471-487.
