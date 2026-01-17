# VacSol-ML(ESKAPE)

**VacSol-ML(ESKAPE)** is a machine learning–driven framework for **vaccine target discovery**, specifically designed to analyze and prioritize candidate proteins from **ESKAPE pathogens**. The tool integrates bioinformatics feature extraction with trained ML models to assist researchers in identifying promising vaccine candidates in a reproducible and scalable manner.

---

## 🔬 Background

ESKAPE pathogens (*Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, Acinetobacter baumannii, Pseudomonas aeruginosa,* and *Enterobacter spp.*) are responsible for a large proportion of antimicrobial-resistant infections worldwide. Identifying effective vaccine targets against these pathogens is a critical challenge.

VacSol-ML(ESKAPE) addresses this by combining:

* Sequence-derived biological and physicochemical features
* Machine learning–based classification
* Automated prediction and confidence scoring

---

## ✨ Key Features

* 🧬 Protein-level vaccine target prediction
* 🤖 Machine learning–based classification pipeline
* 📊 Confidence scoring for predictions
* 🖥️ Immunological and Physiochemical descriptors extraction of the input sequences.
* 📦 Distributed as a pip-installable package

---

## 📦 Installation

### Requirements

* Python **>=3.10, <3.13**

### Install from PyPI

```bash
pip install vacsol-ml
```

This will install VacSol-ML along with its required dependencies.

⚠️ TensorFlow is a required dependency and will be installed automatically during setup. Depending on your system and internet speed, this step may take some time.

---

## 🚀 Usage

VacSol-ML is designed to be used **after installation via pip** through a simple command-line workflow that launches a local web interface.

### Step 1: Run database migrations

After installing the package, initialize the required database tables by running:

```bash
vacsol-ml migrate
```

This command prepares the internal database required for analysis.

---

### Step 2: Launch the web interface

Start the local web server using:

```bash
vacsol-ml web
```

This will launch a local server at:

```
http://127.0.0.1:8000/
```

Open this URL in your browser to access the VacSol-ML web interface and perform analyses interactively.

---

### Step 3: Run analysis via the web UI

Once the server is running:

* Upload or select protein sequence data
* Configure analysis parameters
* Execute the ML-based vaccine target prediction pipeline
* View and download results directly from the interface

> ⚠️ Note: The web interface is intended for **local or server-based use**. It is **not recommended for Google Colab** or notebook-based environments.

---

## 📁 Input Format

* Input sequences must be provided in **FASTA** format
* Each entry should represent a protein sequence

Example:

```fasta
>Protein_1
MKKLLPTAAAGLLLLAAQPAMA...
```

---

## 📤 Output

VacSol-ML generates tabular outputs containing:

* Protein ID
* Predicted class (e.g., Vaccine Candidate / Non-candidate)
* Model confidence score
* Downloadable datasets containing Immunological and Physiochemical features of the input sequences.

---


## 📚 Citation

If you use VacSol-ML in your research, please cite the associated publication:

> *VacSol-ML: An ML-driven framework for vaccine target discovery against ESKAPE pathogens*
> Vaccine (2024). [https://doi.org/10.1016/j.vaccine.2024.126204](https://doi.org/10.1016/j.vaccine.2024.126204)

---

## 👩‍🔬 Author

**Samavi Nasir**
MSc Industrial Biotechnology (ASAB, NUST)
📧 Email: [samavi.nasir@gmail.com](mailto:samavi.nasir@gmail.com)
🔗 Google Scholar: [https://scholar.google.com/citations?hl=en&user=K_okdSMAAAAJ](https://scholar.google.com/citations?hl=en&user=K_okdSMAAAAJ)

---

## 📄 License

This project is released under the **MIT License**. See the `LICENSE` file for details.

---

## ⚠️ Disclaimer

VacSol-ML(ESKAPE) is intended for **research purposes only**. VacSol-ML(ESKAPE) provides in silico predictions to support vaccine target identification. These predictions necessitate downstream experimental validation before clinical translation.

---

## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request via the project repository.
