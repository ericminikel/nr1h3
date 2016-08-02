options(stringsAsFactors=FALSE)
require(binom)

# Penetrance formulae adopted from Minikel 2016 (PMID: 26791950) / Kirov 2014 (PMID: 23992924)
# See original code in https://github.com/ericminikel/prnp_penetrance/blob/master/src/generate_figures.r
# And explanation under https://github.com/ericminikel/prnp_penetrance/blob/master/manuscript.md#lifetime-risk-estimation
penetrance = function(af_case, af_control, baseline_risk) {
  calculated_penetrance = af_case * baseline_risk / af_control
  estimated_penetrance = pmin(1,pmax(0,calculated_penetrance)) # trim to [0,1] support
  return (estimated_penetrance)
}
penetrance_confint = function (ac_case, n_case, ac_control, n_control, baseline_risk) {
  # for a genotypic model, use 1*n_case; for allelic, use 2*n_case
  # here, results are virtually identical.
  case_confint = binom.confint(x=ac_case,n=2*n_case,method='wilson')
  control_confint = binom.confint(x=ac_control,n=2*n_control,method='wilson')
  lower_bound = penetrance(case_confint$lower,control_confint$upper,baseline_risk)
  best_estimate = penetrance(case_confint$mean,control_confint$mean,baseline_risk)
  upper_bound = penetrance(case_confint$upper,control_confint$lower,baseline_risk)
  return ( c(lower_bound, best_estimate, upper_bound) )
}

# Incidence, prevalence, and lifetime risk of MS are reviewed in Alonso & Hernan 2008 (PMID: 18606967)
# See review of lifetime risk in Table 1 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109189/table/T111/
# Best estimates of lifetime risk are 2.5% for women and 1.4% in men.
baseline_risk_women = .0025
baseline_risk_men = .0014

# Details of case series are given in paper: http://www.cell.com/neuron/pdfExtended/S0896-6273(16)30126-X
case_series_n = 2053

# Details of ExAC variant calls are at: http://exac.broadinstitute.org/variant/11-47290147-G-A
exac_ac = 21
exac_n = 33369

# AF in cases if only 1 case is considered part of the case series
1 / (case_series_n * 2)
# AF if 2 cases are considered as part of the case series.
2 / (case_series_n * 2)

# AF in ExAC non-Finnish European population controls
exac_ac / (exac_n * 2)

# Is the variant significantly enriched in (predominantly European?) cases over (European) controls?
# If only 1 case:
fisher.test(matrix(c(1,case_series_n*2-1,exac_ac,exac_n*2-exac_ac),nrow=2),alternative='two.sided')
# If 2 cases:
fisher.test(matrix(c(2,case_series_n*2-2,exac_ac,exac_n*2-exac_ac),nrow=2),alternative='two.sided')

# Supposing that the MS case series in the NR1H3 paper is of predominantly European ancestry and can be
# (very roughly) compared to ExAC non-Finnish Europeans, and using 2 as the AC in cases,
# what is the 95% CI of possible penetrance?
penetrance_confint(2, case_series_n, exac_ac, exac_n, baseline_risk = baseline_risk_women)
penetrance_confint(2, case_series_n, exac_ac, exac_n, baseline_risk = baseline_risk_men)

# For MS, to what odds ratio would 50% penetrance correspond?
odds_hets = 50/50 # if a dominant variant were 50% penetrant, that would mean heterozygotes have 50/50 risk
odds_baseline_women = baseline_risk_women / (1-baseline_risk_women) # odds of developing MS for women in the general population
odds_ratio_for_50percent_penetrance = odds_hets / odds_baseline_women # ratio of the odds
odds_ratio_for_50percent_penetrance

# Here are some additional calculations to support the claim that the inclusion of some inflammatory bowel disease
# exomes in ExAC does not affect interpretation.
exac_males = 33644
exac_females = 27062
# If ExAC ascertainment is neutral with respect to MS, then just by chance, how many individuals in ExAC would be 
# expected to eventually develop MS?
exac_males * baseline_risk_men + exac_females * baseline_risk_women # ~115 individuals
# Now, suppose a worst case scenario, that all 1,675 of the IBD exomes are cases as opposed to controls, and that all are
# women. 
exac_ibd_exomes = 1675
# Suppose these people have 2X the normal level of MS risk - then how many additional individuals who will eventually develop MS
# would be added to ExAC?
exac_ibd_exomes * baseline_risk_women * 2 # ~8 additional individuals
# Or suppose they have 10X the normal level of MS risk - then how many additional individuals who will eventually develop MS
# would be added to ExAC?
exac_ibd_exomes * baseline_risk_women * 10 # ~42 additional individuals
# Conclusion: even if people with IBD had 10X increased risk of MS, the total number of MS cases in ExAC would be increased
# by only about 36% rather than by a factor of hundreds, as would be required to allow NR1H3 R415Q to be Mendelian.

