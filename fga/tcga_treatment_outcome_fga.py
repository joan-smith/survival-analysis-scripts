#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 09:34:32 2021

@author: Joan Smith
"""

#%%

import pandas as pd
import biomarker_survival as surv

#%%
dropbox_dir = '/media/joan/data/Dropbox/comprehensive-tcga-survival/'
clinical = dropbox_dir + 'raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'
fga_path = dropbox_dir + 'FGA_response/TCGA_Fraction_Genome_Altered.xlsx'

#%%
clin = surv.TCGA_CDR_util(clinical)
df = clin.cancer_type_data('*', extra_cols=['treatment_outcome_first_course'])

relevant_outcomes = ['Complete Remission/Response',
 'Progressive Disease',
 'Stable Disease',
 'Partial Remission/Response']
#%%

ctype_counts = df.groupby('type').apply(lambda x:
                                        x[x['treatment_outcome_first_course'].isin(relevant_outcomes)]['treatment_outcome_first_course'].value_counts())
ctype_counts.T.to_csv(dropbox_dir + 'FGA_response/cancer_type_counts.csv')

#%%

fga = pd.read_excel(fga_path).set_index('Patient ID')

ctype_counts_w_fga = df.join(fga, how='inner').groupby('type').apply(lambda x:
                                        x[x['treatment_outcome_first_course'].isin(relevant_outcomes)]['treatment_outcome_first_course'].value_counts())
ctype_counts_w_fga.T.to_csv(dropbox_dir + 'FGA_response/cancer_type_counts_w_fga.csv')
#%%%

def percentile_fga_outcome(ctype_df):
  df = ctype_df.join(fga, how='inner')
  df = df[df['treatment_outcome_first_course'].isin(relevant_outcomes)]
  q25, q75 = df['Fraction Genome Altered'].quantile([.25, .75])
  print(ctype_df.name)
  print(q25, q75)
  lte_q25 = df.groupby('treatment_outcome_first_course').apply(lambda x:
                                                               x[x['Fraction Genome Altered'] <= q25]['treatment_outcome_first_course'].value_counts())
  lte_q25.name = ctype_df.name + '_' + 'lte_q50'
  if lte_q25.shape[0] > 0:
    lte_q25.index = lte_q25.index.droplevel()
  print(lte_q25)

  gte_q75 = df.groupby('treatment_outcome_first_course').apply(lambda x:
                                                               x[x['Fraction Genome Altered'] >= q75]['treatment_outcome_first_course'].value_counts())

  gte_q75.name = ctype_df.name + '_' + 'gte_q50'
  if gte_q75.shape[0] > 0:
    gte_q75.index = gte_q75.index.droplevel()
  print(gte_q75)
  if gte_q75.shape[0] > 0:
    return lte_q25.to_frame().join(gte_q75, lsuffix='_' + lte_q25.name, rsuffix='_' + gte_q75.name)


fga_outcome = df.groupby('type').apply(percentile_fga_outcome)
fga_outcome.index = fga_outcome.index.droplevel(0)


def cleanup_data(outcome_type):
  return outcome_type.melt().dropna().set_index('variable')


out_df = fga_outcome.groupby(level=0).apply(cleanup_data).unstack()
out_df.to_csv(dropbox_dir + 'FGA_response/FGA_outcome_quantiles.csv')



#%%

df = clin.cancer_type_data('*', extra_cols=['tumor_status'])

relevant_status = ['TUMOR FREE', 'WITH TUMOR']

ctype_counts_w_fga = df.join(fga, how='inner').groupby('type').apply(lambda x:
                                        x[x['tumor_status'].isin(relevant_status)]['tumor_status'].value_counts())
ctype_counts_w_fga.T.to_csv(dropbox_dir + 'FGA_response/tumor_status_cancer_type_counts_w_fga.csv')
#%%%

def percentile_fga_outcome_tumor_status(ctype_df):
  df = ctype_df.join(fga, how='inner')
  df = df[df['tumor_status'].isin(relevant_status)]
  q25, q75 = df['Fraction Genome Altered'].quantile([.25, .75])
  print(ctype_df.name)
  print(q25, q75)
  lte_q25 = df.groupby('tumor_status').apply(lambda x:
                                                               x[x['Fraction Genome Altered'] <= q25]['tumor_status'].value_counts())
  lte_q25.name = ctype_df.name + '_' + 'lte_q25'
  if lte_q25.shape[0] > 0:
    lte_q25.index = lte_q25.index.droplevel()
  print(lte_q25)

  gte_q75 = df.groupby('tumor_status').apply(lambda x:
                                                               x[x['Fraction Genome Altered'] >= q75]['tumor_status'].value_counts())

  gte_q75.name = ctype_df.name + '_' + 'gte_q75'
  if gte_q75.shape[0] > 0:
    gte_q75.index = gte_q75.index.droplevel()
  print(gte_q75)
  if gte_q75.shape[0] > 0:
    return lte_q25.to_frame().join(gte_q75, lsuffix='_' + lte_q25.name, rsuffix='_' + gte_q75.name)


fga_outcome = df.groupby('type').apply(percentile_fga_outcome_tumor_status)
fga_outcome.index = fga_outcome.index.droplevel(0)

def cleanup_data(outcome_type):
  return outcome_type.melt().dropna().set_index('variable')


out_df = fga_outcome.groupby(level=0).apply(cleanup_data).unstack()
out_df.to_csv(dropbox_dir + 'FGA_response/FGA_outcome_quantiles_tumor_status.csv')



#%%%
#%%%
#%%%
#%%

tcga_cdr_local = surv.TCGA_CDR_util(clinical)
cancer_types = tcga_cdr_local.cancer_types()

outdir = dropbox_dir + 'FGA_response/'

univariate = {}
for t in cancer_types:
  df = clin.cancer_type_data(t)
  df = df.join(fga, how='inner')
  univariate[t] = (surv.do_cox(df['time'], df['censor'], df['Fraction Genome Altered']))

df = pd.concat((clin.cancer_type_data('COAD'), clin.cancer_type_data('READ')))
df = df.join(fga, how='inner')
univariate['COADREAD'] = (surv.do_cox(df['time'], df['censor'], df['Fraction Genome Altered']))

pd.DataFrame(univariate).to_csv(outdir + '/univariate.csv')


