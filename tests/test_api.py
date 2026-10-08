from pathlib import Path

import pytest, sys, os
from fastapi.testclient import TestClient

PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))
PATIENT_02 = PROJECT_DIR / "Supporting_Files" / "patient_02.fasta"

from app.main import app
@pytest.fixture
def client():
    with TestClient(app) as test_client:
        yield test_client

@pytest.mark.parametrize("patient_number", range(1, 6))
def test_align_patient(client, patient_number):
    patient_id = f"patient_{patient_number:02}"
    fasta_path = PROJECT_DIR / "Supporting_Files" / f"{patient_id}.fasta"

    response = client.post(
        "/api/align",
        files={
            "file": (
                fasta_path.name,
                fasta_path.read_bytes(),
                "text/plain",
            )
        },
    )

    assert response.status_code == 200

    result = response.json()
    assert result["sequence_id"] == patient_id
    assert result["reference_id"] == "REF_001"

    alignment = result["alignment"]
    assert len(alignment["reference_aligned"]) == len(alignment["query_aligned"])
    assert isinstance(result["variants"], list)


@pytest.mark.parametrize("patient_number", range(1, 6))
def test_translate_patient(client, patient_number):
    patient_id = f"patient_{patient_number:02}"
    fasta_path = PROJECT_DIR / "Supporting_Files" / f"{patient_id}.fasta"

    response = client.post(
        "/api/translate",
        files={
            "file": (
                fasta_path.name,
                fasta_path.read_bytes(),
                "text/plain",
            )
        },
        data={"frame": "1"},
    )

    assert response.status_code == 200

    result = response.json()
    assert result["sequence_id"] == patient_id
    assert result["frame"] == 1
    assert isinstance(result["protein"], str)
    assert result["length"] == len(result["protein"])
    assert isinstance(result["has_stop_codon"], bool)
@pytest.mark.parametrize(
    "fasta",
    [
        b"ATGAAATAA\n",                       # Tidak ada header
        b">\nATGAAATAA\n",                    # Header kosong
        b">sample_1\nATG\n>sample_2\nCCC\n",  # Lebih dari satu record
    ],
)
@pytest.mark.parametrize("endpoint", ["/api/align", "/api/translate"])
def test_invalid_fasta_header_is_rejected(client, endpoint, fasta):
    request_data = {
        "files": {"file": ("invalid.fasta", fasta, "text/plain")}
    }

    if endpoint == "/api/translate":
        request_data["data"] = {"frame": "1"}

    response = client.post(endpoint, **request_data)

    assert response.status_code == 400
    assert "detail" in response.json()

@pytest.mark.parametrize("frame, expected_protein, expected_stop", [
    ("1", "MK", True),
    ("2", "", True),
    ("3", "EI", False),
])
def test_translate_reading_frames(client, frame, expected_protein, expected_stop):
    fasta = b">example\nATGAAATAA\n"

    response = client.post(
        "/api/translate",
        files={"file": ("example.fasta", fasta, "text/plain")},
        data={"frame": frame},
    )

    assert response.status_code == 200
    result = response.json()
    assert result["protein"] == expected_protein
    assert result["has_stop_codon"] is expected_stop


def test_translate_n_codon_as_x(client):
    fasta = b">example\nATGNNNTAA\n"

    response = client.post(
        "/api/translate",
        files={"file": ("example.fasta", fasta, "text/plain")},
        data={"frame": "1"},
    )

    assert response.status_code == 200
    assert response.json()["protein"] == "MX"
    assert response.json()["has_stop_codon"] is True


def test_translate_ignores_incomplete_final_codon(client):
    fasta = b">example\nATGAAAC\n"

    response = client.post(
        "/api/translate",
        files={"file": ("example.fasta", fasta, "text/plain")},
        data={"frame": "1"},
    )

    assert response.status_code == 200
    assert response.json()["protein"] == "MK"
    assert response.json()["has_stop_codon"] is False


def test_invalid_fasta_is_rejected(client):
    invalid_fasta = b"this is not a FASTA record"

    for endpoint in ("/api/align", "/api/translate"):
        response = client.post(
            endpoint,
            files={"file": ("invalid.fasta", invalid_fasta, "text/plain")},
        )

        assert response.status_code == 400


@pytest.mark.parametrize("frame", ["0", "4"])
def test_invalid_frame_is_rejected(client, frame):
    fasta = b">example\nATGAAATAA\n"

    response = client.post(
        "/api/translate",
        files={"file": ("example.fasta", fasta, "text/plain")},
        data={"frame": frame},
    )

    assert response.status_code == 422