
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ncov-nowcast

<!-- badges: start -->

<!-- badges: end -->

Contact: Eric Marty (<emarty@uga.edu>)

Contributors:  
Tim Wildauer (<twildauer@blc.edu>) - deconvolution John Drake
(<jdrake@uga.edu>)  
Eric Marty (<emarty@uga.edu>) - visualization, nowcast architecture
Rachel Mercaldo (<mercaldo@uga.edu>) - forecasting

## Contributing

  - Code for functions should go in scripts in the R/ folder.
  - High level scripts or Rmd documents should go at the top level.
  - Temporary data products are in the /data folder.
  - Source data is pulled directly from wikipedia or from
    <https://github.com/CEIDatUGA/COVID-19-DATA>

## Objective

Estimate the current size of the epidemic.

## Rationale

A key problem in making management decisions is estimating the size of
the epidemic. This project aims to estimate the size of the unknown
epidemic. Case notifications are a poor indicator of epidemic size for
several reasons.

  - Case notifications are incomplete (there is under-reporting).
  - Case notifications are primarily associated with patient isolation
    and therefore are not contributing greatly to transmission.
  - 2019-nCov has a significant incubation period. These presymptomatic
    cases are also part of the epidemic.

## Strategy

This project will use a non-parametric approach to deconvolving the case
notification record to construct the actual size of the epidemic as it
existed at past times (backcasting) and then use the backcasted
estimates to feed a forecasting model that “predicts” the present time
(nowcasting). Our nonparametric approach derives from Tim’s REU project.
Backcasting will proceed in two steps using individual-level
observations for the wait time distributions: (1) construct estimated
curve of patients with symptom onset; (2) construct estimated curve of
patients with active infection. From these curves we will use a
statistical model (perhaps time-varying autoregressive model) to predict
the epidemic size at the current time.
