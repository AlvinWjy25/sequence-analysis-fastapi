from typing import Dict, List, Tuple, Any
from pathlib import Path
import sys, os

ROOT_DIR = Path(__file__).parents[2].resolve()
print(f'ROOT_DIR: {ROOT_DIR}')

sys.path.insert(0, str(ROOT_DIR))

from config.config_script import FASTAValidationError, parse_and_validate_fasta

def needleman_wunsch(
    ref_seq: str,
    query_seq: str,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_penalty: float = -2.0,
) -> tuple[float, str, str]:
    """Perform Global Alignment between reference and query sequence using Dynamic Programming.

    Returns (score, ref_aligned, query_aligned).
    """
    n = len(ref_seq)
    m = len(query_seq)

    # Initialize DP scoring matrix (N+1) x (M+1)
    dp = [[0.0] * (m + 1) for _ in range(n + 1)]

    # Fill boundaries with cumulative gap penalties
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap_penalty
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap_penalty

    # Fill DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            ref_base = ref_seq[i - 1]
            query_base = query_seq[j - 1]

            # Calculate match/mismatch score
            if ref_base == query_base:
                diag_score = dp[i - 1][j - 1] + match_score
            else:
                diag_score = dp[i - 1][j - 1] + mismatch_score

            up_score = dp[i - 1][j] + gap_penalty
            left_score = dp[i][j - 1] + gap_penalty

            # Select max score with consistent tie-breaking preference: Diagonal > Up > Left
            dp[i][j] = max(diag_score, up_score, left_score)

    # Traceback phase to construct aligned strings
    ref_aligned = []
    query_aligned = []
    i, j = n, m

    while i > 0 or j > 0:
        current_score = dp[i][j]

        # Check diagonal move
        if i > 0 and j > 0:
            ref_base = ref_seq[i - 1]
            query_base = query_seq[j - 1]
            step_score = match_score if ref_base == query_base else mismatch_score

            if abs(current_score - (dp[i - 1][j - 1] + step_score)) < 1e-6:
                ref_aligned.append(ref_base)
                query_aligned.append(query_base)
                i -= 1
                j -= 1
                continue

        # Check up move (deletion in query / gap in query)
        if i > 0 and abs(current_score - (dp[i - 1][j] + gap_penalty)) < 1e-6:
            ref_aligned.append(ref_seq[i - 1])
            query_aligned.append("-")
            i -= 1
            continue

        # Check left move (insertion in query / gap in reference)
        if j > 0 and abs(current_score - (dp[i][j - 1] + gap_penalty)) < 1e-6:
            ref_aligned.append("-")
            query_aligned.append(query_seq[j - 1])
            j -= 1
            continue

    total_score = dp[n][m]
    return total_score, "".join(reversed(ref_aligned)), "".join(reversed(query_aligned))


def calculate_identity(ref_aligned: str, query_aligned: str) -> float:
    """Calculate sequence identity according to rules.

    identity = identical A/C/G/T columns / total alignment columns.
    N/N is NOT counted as a match.
    """
    total_columns = len(ref_aligned)
    if total_columns == 0:
        return 0.0

    identical_matches = 0
    valid_bases = set("ACGT")

    for r_base, q_base in zip(ref_aligned, query_aligned):
        if r_base in valid_bases and r_base == q_base:
            identical_matches += 1

    return round(identical_matches / total_columns, 2)


def detect_variants(ref_aligned: str, query_aligned: str) -> List[Dict[str, Any]]:
    """Detect SNVs, insertions, and deletions from aligned sequence strings.

    Groups consecutive gaps into single indel events and formats coordinates.
    """
    variants = []
    ref_coord = 0
    i = 0
    n_cols = len(ref_aligned)

    while i < n_cols:
        r_char = ref_aligned[i]
        q_char = query_aligned[i]

        # Case 1: Base match or mismatch (SNV or N)
        if r_char != "-" and q_char != "-":
            ref_coord += 1
            # Call SNV only between differing A/C/G/T bases (exclude N and matches)
            if r_char != q_char and r_char in "ACGT" and q_char in "ACGT":
                variants.append(
                    {
                        "position": ref_coord,
                        "ref": r_char,
                        "alt": q_char,
                        "type": "SNV",
                    }
                )
            i += 1

        # Case 2: Deletion event (gap in query)
        elif q_char == "-":
            start_ref_pos = ref_coord + 1
            deleted_ref_bases = []

            while i < n_cols and query_aligned[i] == "-":
                ref_coord += 1
                deleted_ref_bases.append(ref_aligned[i])
                i += 1

            variants.append(
                {
                    "position": start_ref_pos,
                    "ref": "".join(deleted_ref_bases),
                    "alt": "-",
                    "type": "deletion",
                }
            )

        # Case 3: Insertion event (gap in reference)
        elif r_char == "-":
            insertion_ref_pos = ref_coord  # Preceding reference base position
            inserted_query_bases = []

            while i < n_cols and ref_aligned[i] == "-":
                inserted_query_bases.append(query_aligned[i])
                i += 1

            variants.append(
                {
                    "position": insertion_ref_pos,
                    "ref": "-",
                    "alt": "".join(inserted_query_bases),
                    "type": "insertion",
                }
            )

    # Sort variant events by reference position
    variants.sort(key=lambda x: x["position"])
    return variants


def align_sequences(
    fasta_content: str, reference_id: str, reference_seq: str
) -> Dict[str, Any]:
    """Main wrapper function to execute full pipeline for alignment endpoint."""
    # Step 1: Validate input
    sequence_id, query_seq = parse_and_validate_fasta(fasta_content)

    # Step 2: Global alignment
    score, ref_aligned, query_aligned = needleman_wunsch(reference_seq, query_seq)

    # Step 3: Compute identity
    identity = calculate_identity(ref_aligned, query_aligned)

    # Step 4: Detect variants
    variants = detect_variants(ref_aligned, query_aligned)

    return {
        "sequence_id": sequence_id,
        "reference_id": reference_id,
        "alignment": {
            "score": score,
            "identity": identity,
            "reference_aligned": ref_aligned,
            "query_aligned": query_aligned,
        },
        "variants": variants,
    }