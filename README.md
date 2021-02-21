# circadian
Circadian Rhythm in Critically Ill Patients: Insights from the eICU Database

## Overview

1. Data required to run the analysis:
* patient.csv
* apachePatientResult.csv
* vitalaperiodic.csv
* vitalperiodic.csv
* medication.csv
* infusionDrug.csv
* angus_sepsis.csv

We do not provide the eICU data. You must acquire the data yourself from https://eicu-crd.mit.edu/ after completion of the CITI "Data or Specimens Only Research" course. You will also need to create a materialized view of Angus sepsis criteria following the query provided [here](https://github.com/kseverso/eicu-code/blob/5875ea8e400519d62a7a9d52e3ab94550dc00b41/concepts/angus_sepsis.sql), and save it into a csv file titled ‘angus_sepsis.csv’.

2. Code required to run the analysis:
* dataPreparation.ipynb
* cv_analysis_submission.do
* cv_analysis__mi_submission.do

3. Files generated during analysis:
* patient_data_extraction.csv
* map_data_extraction.csv
* map_median.csv

## First step: dataPreparation.ipynb
This notebook is used to extract and preprocess data. You will need the eICU tables mentioned before in order to run this notebook. 

Notebook outline:
* Import tables
* Inclusion criteria
* Preprocess mean arterial blood pressure (MAP)
* Add medications and infusions
* Calculate the median
* Save the preprocessed data

You can use this notebook if you need to change the inclusion/exclusion criteria (e.g. ICU LOS). Before using the notebook make sure you change the path to the eICU data in 1. Import tables and 3. Preprocess mean arterial blood pressure (MAP).

## Second step: cv_analysis_submission.do

## Third step: cv_analysis__mi_submission.do
