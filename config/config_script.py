import sys, os
from pathlib import Path

ROOT_DIR = Path(__file__).parents[1].resolve()
print(f'ROOT_DIR: {ROOT_DIR}')

SUPPORTING_DIR = Path(ROOT_DIR / 'Supporting_Files')
PATIENT_01_DIR = Path(SUPPORTING_DIR / 'patient_01.fasta')

REFERENCE_DIR = Path(SUPPORTING_DIR / 'reference.fasta')

