import sys, os
from pathlib import Path

ROOT_DIR = Path(__file__).parents[1].resolve()
print(f'ROOT_DIR: {ROOT_DIR}')

SUPPORTING_DIR = Path(ROOT_DIR / 'Supporting_Files')
PATIENT_01_DIR = Path(SUPPORTING_DIR / 'patient_01.fasta')

REFERENCE_DIR = Path(SUPPORTING_DIR / 'reference.fasta')

class FASTAValidationError(Exception):
    """Custom exception for FASTA validation failures."""
    pass

def parse_and_validate_fasta(fasta_content: str) -> tuple[str, str]:
    """Parse a FASTA string and perform strict validation according to the specifications.

    Returns (sequence_id, sequence).
    """
    if not fasta_content or not fasta_content.strip():
        raise FASTAValidationError("FASTA content is empty.")

    lines = [line.strip() for line in fasta_content.splitlines() if line.strip()]

    # Check for header presence and single record restriction
    header_lines = [line for line in lines if line.startswith(">")]
    if len(header_lines) == 0:
        raise FASTAValidationError("Malformed FASTA: Missing header line starting with '>'.")
    if len(header_lines) > 1:
        raise FASTAValidationError(
            f"Multiple records detected ({len(header_lines)}). Only exactly 1 record is allowed."
        )

    # Extract sequence_id (first whitespace-delimited token after '>')
    header = header_lines[0]
    header_body = header[1:].strip()
    if not header_body:
        raise FASTAValidationError("FASTA header is empty after '>' symbol.")

    sequence_id = header_body.split()[0]

    # Reconstruct sequence lines
    seq_lines = [line for line in lines if not line.startswith(">")]

    # Reject embedded whitespace within sequence lines
    for line in seq_lines:
        if any(char.isspace() for char in line):
            raise FASTAValidationError("Embedded whitespace detected within sequence.")

    sequence = "".join(seq_lines).upper()

    # Validate sequence length (1 - 1000 bp)
    if not (1 <= len(sequence) <= 1000):
        raise FASTAValidationError(
            f"Sequence length ({len(sequence)} bp) out of allowed bounds (1 to 1000 bp)."
        )

    # Validate allowed characters (A, C, G, T, N only)
    allowed_bases = set("ACGTN")
    invalid_chars = set(sequence) - allowed_bases
    if invalid_chars:
        raise FASTAValidationError(
            f"Sequence contains invalid character(s): {sorted(list(invalid_chars))}"
        )

    return sequence_id, sequence
