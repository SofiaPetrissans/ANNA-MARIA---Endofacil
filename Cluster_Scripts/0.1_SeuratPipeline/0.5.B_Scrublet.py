#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# PYTHON
# SCRIPT: 5.2º Step for Scrublet
# AUTHOR: Sofía Petrissans Moll
# DATE: 10 junio 2025


# OBJETIVOS:
# Identificar dobletes con Scublet
  # Correr Scrublet en cada muestra
# Guardar objeto resultante 

"""

# LIBRERÍAS .....................................................................

import scrublet as scr
import scipy.io
import pandas as pd
import numpy as np
import os
import sys


# ...............................................................................
# ARGUMENTOS 
# ...............................................................................

if len(sys.argv) < 2:
    print("Uso: python 02_Scrublet_Run.py <input_folder>")
    sys.exit(1)

input_folder = sys.argv[1]

# ...............................................................................
# Cargar Datos 
# ...............................................................................

samples = [f.split(".mtx")[0] for f in os.listdir(input_folder) if f.endswith(".mtx")]

# ...............................................................................
# EJECTUAR SCRUBLET
# ...............................................................................

results = []

for sample in samples:
    counts_matrix = scipy.io.mmread(os.path.join(input_folder, f"{sample}.mtx")).T.tocsr()
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.075)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    df = pd.DataFrame({
        'Barcode': pd.read_csv(os.path.join(input_folder, f"{sample}_barcodes.tsv"), header=None)[0],
        'Scrublet_Score': doublet_scores,
        'Scrublet_Classification': ['Doublet' if x else 'Singlet' for x in predicted_doublets]
    })
    
    df.to_csv(os.path.join(input_folder, f"{sample}_Scrublet.csv"), index=False)
    results.append(df.assign(Sample=sample))

pd.concat(results).to_csv(os.path.join(input_folder, "All_Scrublet_Results.csv"), index=False)
