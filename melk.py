#!/usr/bin/env python
# encoding: utf-8
"""
melk.py

Created by Joan Smith
on 2016-8-15.

Copyright (c) 2016 . All rights reserved.

Given processed breast cancer microarray files, perform one of two survival analyses:
  1. univariate cox analysis of MCM2, PCNA, MKI67, TOP2A, and CCNB1 normalized and averaged
  2. multivariate cox analysis of MCM2, PCNA, MKI67, TOP2A, and CCNB1 normalized and averaged
     and MELK as a second variable.
Note that in both analyses, a GSE file is only included if all six genes appear.

Required inputs:
  1. Processed breast cancer microarray files with excel-safe gene names in the first column,
     organized in folders by GPL. (e.g. GPL96/GSE*.csv)
  2. A file containing four columns: GPL, GSE, Time, and Censor. This allows automatic
    processing of each microarray file.

"""
import getopt
import sys
import re
import os
import glob

import pandas as pd
import numpy as np

from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
from rpy2.robjects import r

rpy2.robjects.numpy2ri.activate()

MICROARRAY_NOTES = '/Users/joans/Dropbox/survival/microarray_notes.csv'
BASE_MICROARRAY_PATH = 'microarrays/'
GENE_SET = ['MELK', 'MCM2', 'PCNA', 'MKI67', 'TOP2A', 'CCNB1']
GENE_SET = set(["'" + gene for gene in GENE_SET])
COX_GENE_SET = ['MCM2', 'PCNA', 'MKI67', 'TOP2A', 'CCNB1']
COX_GENE_SET = set(["'" + gene for gene in COX_GENE_SET])

def do_cox(time, censor, variable, melk=[]):
  surv = importr("survival")

  time = np.array(time, dtype=np.float)
  censor = np.array(censor, dtype=np.float)
  variable = np.array(variable, dtype=np.float)
  if len(melk):
    melk = np.array(melk, dtype=np.float)
  # remove missing data
  skip_cols = []
  for i in range(len(variable)):
    if np.isnan(variable[i]):
      skip_cols.append(i)
    elif np.isnan(time[i]):
      skip_cols.append(i)
    elif np.isnan(censor[i]):
      skip_cols.append(i)
    elif len(melk) and np.isnan(melk[i]):
      skip_cols.append(i)

  variable = np.delete(variable, skip_cols)
  time = np.delete(time, skip_cols)
  censor = np.delete(censor, skip_cols)
  if len(melk):
    melk = np.delete(melk, skip_cols)
  r.assign('time',robjects.FloatVector(time))
  r.assign('censor', robjects.IntVector(censor))
  r.assign('variable', robjects.FloatVector(variable))
  if len(melk):
    r.assign('melk', robjects.FloatVector(melk))

  if len(melk):
    coxuh_output = r('summary(coxph(formula = Surv(time, censor) ~ variable + melk))')
  else:
    coxuh_output = r('summary(coxph(formula = Surv(time, censor) ~ variable))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  patient_count_ind = list(coxuh_output.names).index('n')
  patient_count = coxuh_output[patient_count_ind][0]

  cox_dict = {
      'n': patient_count,
      'z': coeffs.rx('variable', 'z')[0],
      'p': coeffs.rx('variable', 'Pr(>|z|)')[0]
      }
  if len(melk):
    cox_dict['melk-z'] = coeffs.rx('melk', 'z')[0]
    cox_dict['melk-p'] = coeffs.rx('melk', 'Pr(>|z|)')[0]
  return cox_dict

def get_gpl_gse_from_path(f):
  sp = f.split('/')
  gpl = sp[1]
  gse_path = sp[2].split('.')[0]
  gse = gse_path.split('_')[0]
  return gpl, gse

def read_and_check_gse(gse_path, notes_data, multivariate=False):
  gse_data = pd.read_csv(gse_path, index_col=0, low_memory=False, header=None)
  index_set = set(gse_data.index)
  if index_set.issuperset(GENE_SET):
    gpl, gse = get_gpl_gse_from_path(gse_path)
    relevant_genes = gse_data.loc[[unicode(i) for i in COX_GENE_SET]].astype(float)
    relevant_gene_avg = relevant_genes.mean(axis=0)
    columns = notes_data.loc[gpl, gse]
    time = gse_data.loc[columns[u'Time']]
    censor = gse_data.loc[columns[u'Censor']]
    if not multivariate:
      cox = do_cox(time, censor, relevant_gene_avg)
    else:
      melk = gse_data.loc[u'\'MELK']
      cox = do_cox(time, censor, relevant_gene_avg, melk)
    print cox
    return cox, gpl, gse
  else:
    return None, None, None

def do_work(multivariate=False):
  microarray_notes_data = pd.read_csv(MICROARRAY_NOTES)
  microarray_notes_data = microarray_notes_data.set_index(['GPL', 'GSE'])

  formatstr = ''
  outname = ''
  if not multivariate:
    formatstr = '%s,%s,%f,%f,%d\n'
    outname = 'melk_breast.csv'
  else:
    formatstr = '%s,%s,%f,%f,%f,%f,%d\n'
    outname = 'melk_breast_multivariate.csv'

  with open(outname, 'w') as out_file:
    if not multivariate:
      out_file.write('GPL,GSE,p,z,n\n')
    else:
      out_file.write('GPL,GSE,p,z,melk-p,melk-z,n\n')

    breast_cancer_files = glob.glob(BASE_MICROARRAY_PATH + 'GPL*/GSE*.csv')
    for gse_path in breast_cancer_files:
      print gse_path
      cox_dict, gpl, gse = read_and_check_gse(gse_path, microarray_notes_data, multivariate)
      if cox_dict:
        if not multivariate:
          out_file.write(formatstr % (gpl, gse, cox_dict['p'], cox_dict['z'], cox_dict['n']))
        else:
          out_file.write(formatstr % (gpl, gse, cox_dict['p'], cox_dict['z'],
            cox_dict['melk-p'], cox_dict['melk-z'], cox_dict['n']))

def main():
  do_work(multivariate=False)

if  __name__ == '__main__':
  main()
