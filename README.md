# SPP-Study-1

Data analysis for the Swimmer's Phenomics Project - Study 1

<strong> Study Overview </strong>

This study analysed the metabolite profile before, immediately after, 8 and 22 hours after two 7 x 200 m standardised swimming trials and one control (non-exercise) trial. We recruited 14 high-performance male and female swimmers for this study.

<strong> Project Team </strong>

<li> Dr Andrew Govus (La Trobe University, Principal Investigator) </li>
<li> Ms Laine Heidenreich (La Trobe University; Victorian Institute of Sport, PhD Candidate) </li>
<li> Dr Nathan Lawler (Australian National Phenome Centre, Murdoch Univeristy, Co-investigator: Metabolomic Phenotyping) </li>
<li> Dr Chloe Goldsmith (University of Western Sydney, Co-investigator: Epigenetics) </li>
<li> Dr Lachlan Mitchell (Victorian Institute of Sport, Co-investigator: Swimming Physiology) </li>
<li> Ms Louise Cato (Victorian Institute of Sport, Co-investigator: Dietetics) </li>
<li> E/Prof David Pyne (University of Canberra, Co-Investigator: Swimming Physiology) </li>

<br>

<strong> Sample & Data Analysis </strong>

<li><strong> Metabolomic data analysis: </strong> Blood plasma was analysed by untargeted LC-MS (reverse phase and HILIC positive and negative ionisation modes) by Dr Nathan Lawler (Mudroch University) at the Australian National Phenome Centre.
<li> <strong> Bioinformatics </strong>: Bioinformatics for metabolomic data was performed by Laine Heindenreich (La Trobe University), Dr Andrew Govus (La Trobe University) and Dr Nathan Lawler. </li>

<br>

<strong> Project Funding </strong>

<li> Project funding was provided by an internal research support grant from La Trobe University and a High Performance Sport Research Grant from the Australian Institute of Sport. </li>
<br>

<strong> Bioinformatics Data Analysis: </strong>

<li> Data analysis is peformed using the R statistical programming language. </li>
<br>

<strong> Bioinformatics Approach - Metabolomics </strong>

<li> Data cleaning: Data > RSD 30% or >X% zeros was excluded </li>
<li> Inputation of missing data & outlier detection: Random forest inputation of missing data using missRanger </li>
<li> Supervised  multivariate data analysis: Linear Mixed Model ANOVA Simultaneous Components Analysis (LiMM-ASCA) of targeted and untargeted metabolomics using the ALASCA R package. </li>
