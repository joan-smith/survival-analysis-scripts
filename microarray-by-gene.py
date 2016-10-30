#!/usr/bin/env python
# encoding: utf-8
"""
microarray-by-gene.py

Created by Joan Smith
on 2016-5-30.

Copyright (c) 2016 . All rights reserved.

Convert processed microarrays with probe names in the first column to microarrays with standard (excel-safe)
gene names in the first column.

Inputs:
  1. Microarrays to process, organized in folders by GPL (e.g. GPL96/GSE*). Each GPL folder must
     also contain an annotation file named e.g. GPL96.xlsx
  2. Conversion files for Ensembl, ENTREZ, and Genbank that convert from probe to gene.
  3. Excel file with 'Formats' sheet that provides information about which platform was used to
      produce each GPL data set.

Output:
  An identical directory structure with each microarray converted to use gene names.
  Every probe is normalized to zero by its mean. For genes with multiple probes, the value for the
  gene is obtained by taking the normalized mean of the relevant probes.

"""
import sys
import glob
import os
import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import islice
import unicodedata

SURVIVAL_NOTES = 'survival/Survival-notes.xlsx'

def un_unicode(string):
  if string == None:
    return None
  if isinstance(string, str):
    ascii_str =  string
  else:
    ascii_str = unicodedata.normalize('NFKD', string).encode('ascii','ignore')
  ascii_str = ascii_str.strip().lower()
  if ascii_str == 'nan':
    return None
  return ascii_str.lower()

def fix_floats(contents):
  sp = contents.split('.')
  if len(sp) == 2 and sp[1] == '0':
    return sp[0]
  else:
    return contents

def process_gene_ids(col_contents):
  if col_contents == 'nan':
    return None
  col_contents = fix_floats(col_contents)
  col_contents = col_contents.strip()
  if col_contents and (\
    '/' in col_contents or \
    '|' in col_contents or \
    ';' in col_contents):
    vals = col_contents.split('///')
    vals = [val.split('|') for val in vals]
    vals = [x for y in vals for x in y]
    vals = [val.split(';') for val in vals]
    vals = [x for y in vals for x in y]
    vals = [id for id in vals if id]
  else:
    vals = [col_contents]

  vals = [val.upper() for val in vals if val != None]
  return vals

def make_annotation_df(gpl, annotation_type, id_column):
  print gpl, annotation_type, id_column

  annotation_path = os.path.join('data', gpl, gpl + '.xlsx')
  annotations = pd.read_excel(annotation_path, convert_float=True)
  annotations[id_column] = annotations[id_column].astype('str')
  probe_col = annotations.columns[0]

  gene_cardinality = defaultdict(int)
  probes = []
  all_gene_ids = []
  for i,row in annotations.iterrows():
    probe = str(row[probe_col])
    gene_vals = []
    gene_ids = process_gene_ids(str(row[id_column]))
    if probe == None or gene_ids == None:
      continue
    if len(gene_ids) > 0 and probe:
      probes.extend([probe]*len(gene_ids))
      all_gene_ids.extend(gene_ids)
  probe_to_gene_ids_df = pd.DataFrame(data={'ID_REF': probes, 'gene_ids': all_gene_ids})
  probe_to_gene_ids_df = probe_to_gene_ids_df.set_index('ID_REF')

  return probe_to_gene_ids_df

def get_patient_row(microarray_file):
  with open(microarray_file, 'rU') as f:
    row_titles = []
    for i,line in enumerate(f.readlines()):
      row_title = line.split(',')[0]
      if row_title == 'patient' or row_title == 'ID_REF':
        return len(row_titles), row_title
      elif i > 50:
        print 'Failed to find patient row for microarray file: ', microarray_file
        print 'Exiting...'
        sys.exit(1)
      else:
        row_titles.append(row_title)
  return len(row_titles)

def get_new_file_name(gpl, microarray_file):
  basename = os.path.basename(microarray_file)
  split = basename.split('.')
  name = '.'.join(split[:-1]) + '.genes.csv'
  path = os.path.join('by-gene-name', gpl, name)

  return path

def rename_gene_id_to_symbol(df):
    new_columns = list(df.columns)
    gene_id_index = new_columns.index('gene_ids')
    new_columns[gene_id_index] = u'Symbol'
    df.columns = new_columns
    return df

def write_metadata(new_file, patient_row, microarray_file):
  with open(microarray_file, 'rU') as microarray:
    metadata = list(islice(microarray, patient_row))
    metadata = ''.join(metadata)
    new_file.write(metadata)

def write_genes_file(gpl, microarray_file, annotations, conversion_file):
  patient_row, patient_row_title = get_patient_row(microarray_file)
  microarray = pd.read_csv(microarray_file, header=patient_row,  low_memory=False)
  microarray[patient_row_title] = microarray[patient_row_title].astype('str')
  microarray = microarray.set_index(patient_row_title)
  normed_microarray = microarray.sub(microarray.mean(axis=1), axis=0)

  if conversion_file == None:
    joined_with_gene_id = normed_microarray.join(annotations, how='right')
    joined_with_gene_id = joined_with_gene_id.reset_index()
    joined_with_gene_id = joined_with_gene_id.dropna(how='any', subset=['gene_ids'])
    joined_with_gene_id['gene_ids'] = '\'' + joined_with_gene_id['gene_ids']
    joined_with_gene_id = rename_gene_id_to_symbol(joined_with_gene_id)
    gene_groups = joined_with_gene_id.groupby(u'Symbol').agg(np.mean)
  else:
    joined_with_gene_id = normed_microarray.join(annotations)
    joined_with_gene_id = joined_with_gene_id.reset_index()
    joined_with_gene_id = joined_with_gene_id.dropna(how='any', subset=['gene_ids'])
    joined_with_gene_id = joined_with_gene_id.set_index('gene_ids')

    conversion = pd.read_csv(conversion_file, dtype='str', low_memory=False)
    conversion = conversion.dropna(how='any')
    conversion = conversion.set_index(u'Annotation')

    joined = joined_with_gene_id.join(conversion[u'Symbol'])
    joined = joined.dropna(subset=[u'Symbol'])
    gene_groups = joined.groupby(u'Symbol').agg(np.mean)

  new_file_name = get_new_file_name(gpl, microarray_file)
  print new_file_name
  if gene_groups.size == 0:
    print '   Error: no data'
  with open(new_file_name, 'w') as new_file:
    write_metadata(new_file, patient_row, microarray_file)
    gene_groups.to_csv(new_file, index_label='ID_REF')

def do_work():
  conversion_files = {
      'Ensembl': 'conversions/ensembl to HGNC.csv',
      'Genbank': 'conversions/genbank to HGNC.csv',
      'ENTREZ': 'conversions/Entrez to HGNC.csv'
  }
  # Note: ORF has no conversion file.

  survival_notes = pd.read_excel(SURVIVAL_NOTES,
      sheetname='Formats',
      index_col=u'Array format')

  for gpl, gpl_row in survival_notes.iterrows():
    try:
      os.mkdir(os.path.join('by-gene-name', gpl))
    except OSError as err:
      pass
    annotation_dict = make_annotation_df(gpl,
        gpl_row[u'Merge into this file'],
        gpl_row[u'Column name'])
    microarray_files =  glob.glob(os.path.join('data', gpl, 'GSE*.csv'))
    if gpl_row[u'Merge into this file'] == 'ORF':
      conversion_file = None
    else:
      conversion_file = conversion_files[gpl_row[u'Merge into this file']]
    for microarray_file in microarray_files:
      write_genes_file(gpl, microarray_file, annotation_dict, conversion_file)

def main():
  do_work()

if __name__ == '__main__':
  main()

