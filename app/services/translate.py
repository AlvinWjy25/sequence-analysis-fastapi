from typing import Dict, List, Tuple, Any

import sys, os
from pathlib import Path

ROOT_DIR = Path(__file__).parents[2].resolve()
print(f'ROOT_DIR: {ROOT_DIR}')

sys.path.insert(0, str(ROOT_DIR))

from config.config_script import FASTAValidationError, parse_and_validate_fasta

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: AMINO[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in BASES for b in BASES for c in BASES
    )
}

def translate_sequence(sequence: str, frame: int):
    seq = sequence[frame - 1:]
    protein = []
    has_stop = False
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i + 3]
        if "N" in codon:
            protein.append("X")
            continue
        aa = CODON_TABLE[codon]
        if aa == "*":
            has_stop = True
            break
        protein.append(aa)
    return "".join(protein), has_stop

def translate(
    fasta_content: str, reference_id: str, reference_seq: str, frame: int,
) -> Dict[str, Any]:
    sequence_id, sequence = parse_and_validate_fasta(fasta_content)
    protein, has_stop = translate_sequence(sequence, frame)

    return {
        "sequence_id": sequence_id,
        "frame": frame,
        "protein": protein,
        "length": len(protein),
        "has_stop_codon": has_stop,
    }