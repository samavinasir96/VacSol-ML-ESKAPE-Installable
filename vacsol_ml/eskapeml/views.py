from pathlib import Path
import os
import sys
import time
import csv
import zipfile
import logging
import subprocess

import numpy as np
import pandas as pd
import joblib
import iedb
import biolib

from Bio import SeqIO
from keras.models import load_model

from django.shortcuts import render, redirect
from django.http import (
    HttpResponse, HttpResponseRedirect, JsonResponse
)
from django.http import FileResponse, Http404
from django.views.decorators.csrf import csrf_exempt
from django.core.cache import cache

import shutil

# =========================
# PATH SETUP (CRITICAL)
# =========================
APP_DIR = Path(__file__).resolve().parent            # eskapeml/
BASE_DIR = APP_DIR.parent                            # vacsol_ml/
DATA_DIR = BASE_DIR / "Data"
TEMP_DIR = APP_DIR / "Temp_Results"
RESULTS_DIR = APP_DIR / "Analysis_Results"

TEMP_DIR.mkdir(exist_ok=True, parents=True)
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

SEQUENCE_FILE = APP_DIR / "sequences.fasta"

IFEATURE_SCRIPT = DATA_DIR / "iFeature" / "iFeature.py"
SCALER_PATH = APP_DIR / "eskape_scaler.pkl"
ENCODER_PATH = APP_DIR / "encoder_model.keras"
ENSEMBLE_MODEL_PATH = APP_DIR / "eskape_ensemble_model.pkl"

# =========================
# LOGGING
# =========================
LOG_FILE = APP_DIR / "vacsol_ml.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

# =========================
# PROGRESS UTILITIES
# =========================


from django.http import JsonResponse

def get_progress(request):
    job_id = request.session.get("job_id")
    if not job_id:
        return JsonResponse({"progress": 0, "logs": []})

    return JsonResponse({
        "progress": cache.get(f"progress_{job_id}", 0),
        "logs": cache.get(f"logs_{job_id}", [])
    })


# =========================
# BASIC VIEWS
# =========================
@csrf_exempt
def home(request):
    return render(request, "index.html")

def faqs(request):
    return render(request, "faqs.html")

def glossary(request):
    return render(request, "glossary.html")

def processing_view(request):
    return render(request, "processing.html")

# =========================
# FILE UPLOAD
# =========================
from django.shortcuts import redirect
from django.core.cache import cache
import threading
import uuid

def sanitize_fasta_headers(raw_fasta: str, header_len: int = 10):
    sanitized = []
    current_seq = []
    header_map = {}

    for line in raw_fasta.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_seq:
                sanitized.append("".join(current_seq))
                current_seq = []

            original = line[1:]
            truncated = original.split()[0][:header_len]
            sanitized.append(f">{truncated}")
            header_map[truncated] = original
        else:
            current_seq.append(line)

    if current_seq:
        sanitized.append("".join(current_seq))

    return "\n".join(sanitized), header_map

def upload_sequence(request):
    if request.method == "POST":

        raw_sequence = request.POST.get("sequence")
        if not raw_sequence:
            return JsonResponse({"error": "No sequence provided"}, status=400)

        sanitized_fasta, protein_ids = sanitize_fasta_headers(raw_sequence)

        SEQUENCE_FILE.write_text(sanitized_fasta)

        job_id = str(uuid.uuid4())
        request.session["job_id"] = job_id

        cache.set(f"progress_{job_id}", 0, timeout=3600)
        cache.set(f"logs_{job_id}", ["Upload started"], timeout=3600)

        threading.Thread(
            target=calculate_features,
            args=(request, job_id),
            daemon=True
        ).start()

        return redirect("eskapeml:processing")

    return render(request, "upload_vacsolml.html")


# =========================
# CORE PIPELINE
# =========================
def calculate_features(request, job_id):

    def advance(step, msg):
        percent = int((step / 20) * 100) 
        cache.set(f"progress_{job_id}", percent, 3600)

        logs = cache.get(f"logs_{job_id}", [])
        logs.append(msg)
        cache.set(f"logs_{job_id}", logs, 3600)

        logger.info(f"{percent}% - {msg}")

    advance(1, "Sequence uploaded")

    # -------------------------
    # iFeature
    # -------------------------
    advance(3, "Running iFeature")

    ifeature_outputs = [
        ("Moran", "ifeature3.tsv"),
        ("Geary", "ifeature4.tsv"),
        ("NMBroto", "ifeature5.tsv"),
        ("CTriad", "ifeature9.tsv"),
        ("KSCTriad", "ifeature10.tsv"),
        ("SOCNumber", "ifeature11.tsv"),
        ("QSOrder", "ifeature12.tsv"),
        ("CKSAAP", "ifeature15.tsv"),
        ("CKSAAGP", "ifeature19.tsv"),
    ]

    for ftype, outfile in ifeature_outputs:
        subprocess.run([
            sys.executable,
            str(IFEATURE_SCRIPT),
            "--file", str(SEQUENCE_FILE),
            "--type", ftype,
            "--out", str(TEMP_DIR / outfile)
        ], check=True)

    dfs = [pd.read_table(TEMP_DIR / f[1]) for f in ifeature_outputs]
    df = pd.concat(dfs, axis=1).loc[:, ~pd.concat(dfs, axis=1).columns.duplicated()]
    df.rename(columns={"#": "Protein_ID"}, inplace=True)
    df.to_csv(TEMP_DIR / "ifeatures.csv", index=False)
    df.to_csv(RESULTS_DIR / "Descriptors_by_iFeature.csv", index=False)

    advance(4, "Physicochemical features extracted")

    # -------------------------
    # IEDB EPITOPES
    # -------------------------
    sequences = SeqIO.to_dict(SeqIO.parse(SEQUENCE_FILE, "fasta"))

    mhci_scores, mhci_ranks = [], []
    mhcii_scores, mhcii_ranks = [], []
    bcell_scores, surface_scores, antigenicity_scores = [], [], []

    for record in sequences.values():
        seq = str(record.seq)

        advance(6, "Running NetMHCpan")

        mhci = iedb.query_mhci_binding("recommended", seq, "HLA-A*02:01", 9)
        mhci["percentile_rank"] = pd.to_numeric(mhci["percentile_rank"], errors="coerce")
        mhci_f = mhci[mhci["percentile_rank"] <= 0.5]
        #mhci_scores.append(mhci_f["score"].mean() if not mhci_f.empty else 0)
        scores = pd.to_numeric(mhci_f["score"], errors="coerce")
        mhci_scores.append(scores.mean() if not scores.empty else 0)
        #mhci_ranks.append(mhci_f["percentile_rank"].mean() if not mhci_f.empty else 0)
        ranks = pd.to_numeric(mhci_f["percentile_rank"], errors="coerce")
        mhci_ranks.append(ranks.mean() if not ranks.empty else 0)
        mhci_epitopes = mhci_f['peptide'].values
        
        advance(8, "Running NetMHCIIpan")

        mhcii = iedb.query_mhcii_binding("recommended", seq, "HLA-DRB1*01:01", 15)
        mhcii["rank"] = pd.to_numeric(mhcii["rank"], errors="coerce")
        mhcii_f = mhcii[mhcii["rank"] <= 1]
        #mhcii_scores.append(mhcii_f["score"].mean() if not mhcii_f.empty else 0)
        scoresii = pd.to_numeric(mhcii_f["score"], errors="coerce")
        mhcii_scores.append(scoresii.mean() if not scores.empty else 0)

        #mhcii_ranks.append(mhcii_f["rank"].mean() if not mhcii_f.empty else 0)
        rankii = pd.to_numeric(mhcii_f["rank"], errors="coerce")
        mhcii_ranks.append(rankii.mean() if not rankii.empty else 0)
        mhcii_epitopes = mhcii_f['peptide'].values

        def safe_mean(series):
            numeric = pd.to_numeric(series, errors="coerce")
            return numeric.mean() if not numeric.empty else 0

        advance(10, "Running B-cell Epitopes Query")

        bcell_df = iedb.query_bcell_epitope("Bepipred", seq, 9)
        bcell_scores.append(
        safe_mean(bcell_df["Score"])
        )

        surface_df = iedb.query_bcell_epitope("Emini", seq, 9)
        surface_scores.append(
            safe_mean(surface_df["Score"])
        )

        antigenicity_df = iedb.query_bcell_epitope("Kolaskar-Tongaonkar", seq, 9)
        antigenicity_scores.append(
            safe_mean(antigenicity_df["Score"])
        )


    advance(13, "Epitope prediction completed")

    # -------------------------
    # SIGNALP
    # -------------------------
    advance(15,"Running SignalP")

    signalp = biolib.load("DTU/SignalP_6")
    job = signalp.cli(args=f"--fastafile {SEQUENCE_FILE} --output_dir output")
    job.save_files(TEMP_DIR / "signalP", overwrite=True)

    sp_df = pd.read_table(
        TEMP_DIR / "signalP" / "prediction_results.txt",
        sep="\t",
        skiprows=1
    )
    sp_df.columns = ["ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "TATLIPO(Sec/SPII)", "PILIN(Sec/SPIII)", "CS Position"]
    #sp_df.drop(index=sp_df.index[0], axis=1, inplace=True)

    SP = sp_df['SP(Sec/SPI)'].values
    LIPO = sp_df['LIPO(Sec/SPII)'].values
    TAT = sp_df['TAT(Tat/SPI)'].values
    TATLIPO = sp_df['TATLIPO(Sec/SPII)'].values
    PILIN = sp_df['PILIN(Sec/SPIII)'].values
    OTHER = sp_df['OTHER'].values

    advance(16, "Signal peptide prediction completed")

    # -------------------------
    # MERGE & SCALE
    # -------------------------
    features = pd.read_csv(TEMP_DIR / "ifeatures.csv")
    features["antigenicity_1"] = antigenicity_scores
    features["b_cells_probability_score"] = bcell_scores
    features["mhci_probability_score"] = mhci_scores
    features["mhci_rank"] = mhci_ranks
    features["mhcii_rank"] = mhcii_ranks
    features["mhcii_score"] = mhcii_scores
    features["surface_probability"] = surface_scores

    features["signal_peptide_SP"] = SP
    features["signal_peptide_LIPO"] = LIPO
    features["signal_peptide_TAT"] = TAT
    features["signal_peptide_TATLIPO"] = TATLIPO
    features["signal_peptide_PILIN"] = PILIN
    features["signal_peptide_OTHER"] = OTHER

    if request.method == "POST":
        selected_pathogen = request.POST.get("pathogen")

        pathogen_encoding = {
            "Staphylococcus aureus": 5,
            "Enterococcus faecium": 2,
            "Pseudomonas aeruginosa": 4,
            "Acinetobacter baumannii": 0,
            "Klebsiella pneumoniae": 3,
            "Enterobacter": 1,
        }

        if selected_pathogen not in pathogen_encoding:
            raise ValueError(f"Invalid pathogen selected: {selected_pathogen}")

        encoded_pathogen = pathogen_encoding[selected_pathogen]
        features["Organism.1"] = encoded_pathogen

        features.to_csv(TEMP_DIR / "finalfeatures.csv", index=False)

        protein_ids = features["Protein_ID"]
        dataset_numeric = features.drop(columns=["Protein_ID"])

        with open(SCALER_PATH, "rb") as f:
            scaler = joblib.load(f)

        rescaled = scaler.transform(dataset_numeric)

        rescaled_df = pd.DataFrame(
            rescaled,
            columns=dataset_numeric.columns
        )
        rescaled_df.insert(0, "Protein_ID", protein_ids.values)
        rescaled_df.to_csv(TEMP_DIR / "scaledfinalfeatures.csv", index=False)

        biological_df =  pd.DataFrame({
            'Protein_ID': protein_ids.values,
            'antigenicity_1': antigenicity_scores,
            'b_cells_probability_score': bcell_scores,
            'mhci_probability_score': mhci_scores,
            'mhci_rank': mhci_ranks,
            'mhcii_rank': mhcii_ranks,
            'mhcii_score': mhcii_scores,
            'surface_probability': surface_scores,
        })
        biological_df = biological_df.round(6)


        threshold_adhesion = [
            'Threshold for adhesion probability = 1 Or greater than 1',
        ]
        threshold_adesion_df = pd.DataFrame({'Thresholds': threshold_adhesion})

        
        df_mhci_epitopes = pd.DataFrame({
                                            'Epitope_Type': ['mhci'] * len(mhci_epitopes),
                                            'Epitope_Sequence': mhci_epitopes})
        
        df_mhcii_epitopes = pd.DataFrame({
                                            'Epitope_Type': ['mhcii'] * len(mhcii_epitopes),
                                            'Epitope_Sequence': mhcii_epitopes})
        
        comments = [
            'Threshold for strong "MHC Class I" binder: % Rank = Lower than 0.5',
            'Threshold for strong "MHC Class II" binder: % Rank = Lower than 1',
            'Threshold for strong "B cell epitope": Score = Greater than 0.5',
            'Antigenicity = Value close to 1 or Greater than 1 indicates high antigenicity'
        ]
        comments_df = pd.DataFrame({'Thresholds': comments})

        df_biological = pd.concat([biological_df, sp_df, df_mhci_epitopes, df_mhcii_epitopes, comments_df], ignore_index=True)
            
        df_biological.to_csv(RESULTS_DIR / "Biological_Properties.csv", index=False)


    advance(20, "ML predictions")

# =========================
# RESULTS VIEW
# =========================
def get_results(request):
    features = pd.read_csv(TEMP_DIR / "scaledfinalfeatures.csv")
    encoder = load_model(ENCODER_PATH)
    encoded = encoder.predict(features.drop(columns=["Protein_ID"]))

    with open(ENSEMBLE_MODEL_PATH, "rb") as f:
        model = joblib.load(f)

    # Predictions
    preds = (
        model["Random_Forest"].predict(encoded) +
        model["Gradient_Boosting"].predict(encoded) +
        model["Logistic_Regression"].predict(encoded)
    ) >= 1

    # Probabilities (averaged)
    prob_rf = model["Random_Forest"].predict_proba(encoded)[:,1]
    prob_gb = model["Gradient_Boosting"].predict_proba(encoded)[:,1]
    prob_lr = model["Logistic_Regression"].predict_proba(encoded)[:,1]
    avg_prob = (prob_rf + prob_gb + prob_lr) / 3

    results = features.copy()
    results["Prediction"] = preds.astype(int)
    results["Probability_Class_1"] = avg_prob
    results["Probability_Class_0"] = 1 - avg_prob

    # confidence aligned with prediction
    results["Confidence"] = np.where(
        results["Prediction"] == 1,
        results["Probability_Class_1"],
        results["Probability_Class_0"]
    )

    results["Confidence_percent"] = results["Confidence"] * 100


    return render(
        request,
        "results_vacsolml.html",
        {"eskapeml_results": results.to_dict(orient="records")}
    )

def eskapeml_download_csv(request):

    df = pd.read_csv(TEMP_DIR / "scaledfinalfeatures.csv")
    encoder = load_model(ENCODER_PATH)
    encoded = encoder.predict(df.drop(columns=["Protein_ID"]))

    with open(ENSEMBLE_MODEL_PATH, "rb") as f:
        model = joblib.load(f)
    
    prob = (
        model["Random_Forest"].predict_proba(encoded)[:, 1] +
        model["Gradient_Boosting"].predict_proba(encoded)[:, 1] +
        model["Logistic_Regression"].predict_proba(encoded)[:, 1]
    ) / 3

    results_df = pd.DataFrame({
        "Protein_ID": df["Protein_ID"],
        "Prediction": (prob >= 0.5).astype(int),
        "Probability_Class_1": prob.round(4)
    })

    response = HttpResponse(content_type="text/csv")
    response["Content-Disposition"] = 'attachment; filename="VacSolML_Predictions.csv"'
    results_df.to_csv(response, index=False)
    return response


def download_biological_csv(request):
    path = RESULTS_DIR / "Biological_Properties.csv"
    if not path.exists():
        raise Http404("File not found")

    return FileResponse(path.open("rb"), as_attachment=True, filename="Biological_Properties.csv")

def download_physiochemical_csv(request):
    path = RESULTS_DIR / "Descriptors_by_iFeature.csv"
    if not path.exists():
        raise Http404("File not found")

    return FileResponse(path.open("rb"), as_attachment=True, filename="Descriptors_by_iFeature.csv")

@csrf_exempt
def track_citation(request):
    if request.method == "POST":
        # log to DB / file / analytics system
        return JsonResponse({"status":"ok"})

def clear_temp(request):
    if TEMP_DIR.exists():
        for item in TEMP_DIR.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
    return redirect('eskapeml:home') 