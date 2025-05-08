# -*- coding: utf-8 -*-
# %% ============================================================================
#
# Extraction and Calculation of Protein Properties from Uniprot Database
# Author: Alexa Canchola
# Advisor: Wei-Chun Chou
# Date: April 22, 2025
# ==============================================================================
"""
This Python script provides a comprehensive workflow for the extraction of molecular weight (MW) and amino acid sequence information from the UniProt Protein Knowledgebase. 
Additional functionality to calculate protein properties including:
    - Isoelectric point (pI)
    - Grand Average of Hydropathicity (GRAVY)
using the BioPython (v1.85) package
""" 
# %% ============================================================================
# Install Required Libraries
# ==============================================================================
!pip install biopython

# %% ============================================================================
# Import Libraries
# ==============================================================================
# Standard Libraries for Data Handling
import pandas as pd
import numpy as np


# Libraries for Uniprot Query
import requests as r
import re
from urllib.request import urlopen
from ast import AnnAssign
import warnings
warnings.filterwarnings("ignore")

# Bio.SeqUtils Package
# Used to calculate protein information from amino acid sequence
from Bio.SeqUtils.ProtParam import ProteinAnalysis



# %% ============================================================================
#Import & Preprocess Data
# ============================================================================

# Load data; data should only include valid UniProt Entry IDs
proteins=pd.read_csv('/content/proteinIDs.csv')

proteinids = proteins['ID'].to_list() #requires list format
print(proteisnids)

# %% ============================================================================
# Define Data Retrieval Functions
# ============================================================================

#Function to query MW and AA sequence for proteins in list
def get_mw_and_seq_from_uniprot_id(ID):
    ID = ID.strip()
    url = f"http://www.uniprot.org/uniprot/{ID}.txt"
    try:
        response = urlopen(url)
        data = response.read().decode()

        # Optional: Print the raw data for inspection
        # print(f"Raw data for {ID}:\n{data}\n{'-'*40}")

        # Ensure the response contains the sequence information
        if 'SQ ' in data:
            result = data.split('SQ ')[-1]

            # Use regex to find molecular weight (MW)
            mw_match = re.search(r"(\d+)\s*MW", result)
            if mw_match:
                mw = int(mw_match.group(1))
            else:
                print(f"Molecular weight not found for {ID}")
                mw = None

            #use regex to find amino acid length (AA)
            AA_match = re.search(r"(\d+)\s*AA", result)
            if AA_match:
                AA = int(AA_match.group(1))
            else:
                print(f"Amino acid length not found for {ID}")

            # Extract sequence
            sequence_lines = result.split('\n')[1:]
            sequence = ''.join([line.strip().replace(' ', '') for line in sequence_lines if line and not line.startswith('//')])

            print(f"Extracted sequence for ID {ID}: {sequence}")
            return mw, sequence, AA
        else:
            print(f"Sequence data not found for UniProt ID {ID}")
            return None, None, None

    except Exception as e:
        print(f"Error fetching data for UniProt ID {ID}: {e}")
        return None, None, None

#Function to clean formatting of extracted AA sequences
def clean_sequence(sequence):
  clean_sequence = ''.join(sequence.split()).upper()
  valid_amino_acids = set( "ARNDCEQGHILKMFPSTWYV")

  #check for invalid characters
  if not all(char in valid_amino_acids for char in clean_sequence):
      raise ValueError("Sequence contains invalid characters.")

  return clean_sequence
# %% ============================================================================
# Define Functions for pI and GRAVY calculation
# ============================================================================

#define function for calculating
def calculate_isoelectric_point (sequence):
  try:
    sequence = clean_sequence(sequence)
    analysis = ProteinAnalysis(sequence)
    pI = analysis.isoelectric_point()
    return pI
  except Exception as e:
    print(f"Error calculating pI for sequence: {e}")
    return None

def calculate_gravy(sequence):
    try:
        sequence = clean_sequence(sequence)
        analysis = ProteinAnalysis(sequence) 
        gravy = analysis.gravy() # calculate the GRAVY according to Kyte and Doolitle, 1982
        return gravy
    except Exception as e:
        print(f"Error calculating GRAVY for sequence: {e}")
        return None

#create dataframe for calculated properties
ProteinInfo = pd.DataFrame(proteinids)
ProteinInfo.rename(columns={ProteinInfo.columns[0]: 'ID'}, inplace=True)
ProteinInfo['MW'] = None
ProteinInfo['AA'] = None
ProteinInfo['pI'] = None
ProteinInfo['gravy'] = None

# %% ============================================================================
# Define Functions for pI and GRAVY calculation
# ============================================================================

# Loop through protein ID list to extract relevant data; if data is available, calculate specified properties
for index, row in ProteinInfo.iterrows():
    mw, seq, AA = get_mw_and_seq_from_uniprot_id(row['ID'])
    if mw and seq:
       ProteinInfo.at[index, 'MW'] = mw
       ProteinInfo.at[index,'AA'] = AA
       ProteinInfo.at[index, 'pI'] = calculate_isoelectric_point(seq)
       ProteinInfo.at[index, 'gravy'] = calculate_gravy(seq)

# Check dataframe
print(ProteinInfo.head(30))