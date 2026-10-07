import os, sys
from pathlib import Path
ROOT_DIR = Path(__file__).resolve().parents[1]
print(ROOT_DIR)

sys.path.insert(0, str(ROOT_DIR))

from fastapi import (FastAPI, UploadFile,  
                     File, HTTPException, status)

from app.services.alignment import align_sequences, FASTAValidationError, parse_and_validate_fasta
from config.config_script import REFERENCE_DIR


from pydantic import BaseModel
from contextlib import asynccontextmanager

reference_store = {}

async def lifespan(app: FastAPI):
    ref_path = REFERENCE_DIR
    if not ref_path.exists():
        raise RuntimeError(f"Reference file not found at {ref_path}")

    with open(ref_path, "r", encoding="utf-8") as f:
        ref_content = f.read()

    # Parse and validate reference sequence using your existing FASTA parser
    ref_id, ref_seq = parse_and_validate_fasta(ref_content)

    # Store in memory
    reference_store["id"] = ref_id
    reference_store["sequence"] = ref_seq
    print(f"Loaded reference sequence '{ref_id}' ({len(ref_seq)} bp) into memory.")

    yield  # Application is running and serving requests

    # --- SHUTDOWN LOGIC ---
    # Cleanup tasks if any (optional)
    reference_store.clear()

app = FastAPI (
    title = "DNA Sequence Analysis RestAPI",
    description = "Bioinformatics Internship Kable Technical Test",
    version = "1.0.0",
    lifespan=lifespan
)

@app.get("/")
def read_root():
    return {"message": "DNA Analysis API is up and running!"}

@app.post("/api/test-upload")
async def test_upload(file: UploadFile = File(...)):
    content = await file.read()
    fasta_str = content.decode("utf-8")

    return {
        "filename": file.filename,
        "content_length": len(fasta_str),
        "preview": fasta_str[:50] if fasta_str else ""
    }

@app.post("/api/align")
async def align_fasta(file: UploadFile = File(...)):
    if not file or not file.filename:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail="Missing file."
        )

    try:
        content = await file.read()
        fasta_content = content.decode("utf-8")

        # Use the reference sequence stored in memory (No disk read needed!)
        result = align_sequences(
            fasta_content=fasta_content,
            reference_id=reference_store["id"],
            reference_seq=reference_store["sequence"],
        )

        return result

    except FASTAValidationError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    
    except UnicodeDecodeError:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Uploaded file must be a valid text FASTA file.",
        )