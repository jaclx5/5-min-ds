---
title: COVID Vaccines How Much 1% Matters?
layout: post
author: jaclx5
---

_1% increase on vaccinated population correlates with 352 saved lives each week in Europe._


> Check [the companion notebook](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/covid_vaccination/covid_vaccination.ipynb) for the complete data analysis.



A few days ago someone posted the following chart on Twitter:

<div align="center">
	<img src="/images/covid_vaccination/covid_vaccination-ext_1.png" alt="COVID-19 Vaccination" width="400"/>
</div>

The chart compares the percentages of fully vaccinated adult population with the COVID-19 death rates in EU countries.

The negative correlation between vaccination and COVID-19 death rates seems pretty obvious. However, the actual decrease in number of deaths, for each each percent point increase in vaccinated population, is not so easy to estimate just by looking at the chart.

In this post we try to obtain an estimate (rough) for this quantity from the data.


First we went for the original data sources at the [European Centre for Disease Prevention and Control](https://www.ecdc.europa.eu/), and reproduced the initial chart:

<div align="center">
	<img src="/images/covid_vaccination/covid_vaccination-fig_1.png" alt="COVID-19 Vaccination" width="700"/>
</div>

As expected results are pretty close except for: (i) small changes in the order of  countries (we don't really know the exact week of the initial chart data); (ii) in our analysis we used the whole population instead of the adult population, so percentages are smaller.

Comming back to the original question, we performed a simple linear regression in order to "predict" the 14-days average number of deaths per 1 million inhabitants given the percentage of fully vaccinated population:

<div align="center">
	<img src="/images/covid_vaccination/covid_vaccination-fig_2.png" alt="COVID-19 Vaccination" width="400"/>
</div>

Although with some caveats (see bellow) we get a regression line with a slope of around -570, which means that, in average, an increase of a single 1% in fully vaccinated population correlates with a decrease of 352 deaths per week:

<div align="center">
	<img src="/images/covid_vaccination/covid_vaccination-fig_3.png" alt="COVID-19 Vaccination" width="400"/>
</div>


__352 people is not small deal! Get vaccinated, please!__


## Caveats and Sources

The goal on this post was to obtain a rough estimation of the correlation between vaccination and deaths rates in Europe. Rather than a precise number, which would be impossible to compute from this data alone, we want to have a sense of the order of magnitude of the effect of vaccination on the death rates.




