# Least Squares Approximation of ECG Waveforms

This repository contains the implementation and experimental code for the paper  
**“Least Squares Approximation of ECG Waveforms for Cardiac Signal Analysis”**,  
developed as part of the **Linear Algebra and Geometry** course.

The project demonstrates how **least squares approximation** can be applied to real ECG signals by modeling them as vectors in a finite-dimensional space and projecting them onto polynomial subspaces.

---

## 📌 Overview

Electrocardiogram (ECG) signals are time-series measurements of cardiac electrical activity. From a mathematical perspective, each ECG segment can be represented as a vector in \( \mathbb{R}^n \).  
This project applies **polynomial least squares approximation** to ECG segments from the **MIT-BIH Arrhythmia Database**, with the goal of:

- illustrating least squares as a **projection operation**
- analyzing the effect of polynomial degree on approximation accuracy
- visualizing approximation errors geometrically
- extracting simple, interpretable features from approximation results

The implementation focuses on **mathematical modeling and analysis**, not medical diagnosis.

---

## 📂 Repository Structure

```
.
├── src/
│   ├── ecganalysis.py   # Main implementation script
│
├── figures/
│   ├── ecg_vs_ls_degree6.png     # Example approximation plot
│   ├── rmse_vs_degree.png        # RMSE vs polynomial degree
│
└── README.md
```

---

## 🧠 Methodology Summary

1. **Dataset**
   - MIT-BIH Arrhythmia Database (PhysioNet)
   - Record 101, MLII lead
   - Sampling rate: 360 Hz

2. **Preprocessing**
   - Baseline removal (mean subtraction)
   - Noise reduction (moving average filter)
   - Z-score normalization
   - Fixed-length segmentation

3. **Mathematical Model**
   - ECG segment represented as a vector \( \mathbf{y} \in \mathbb{R}^n \)
   - Polynomial basis constructed using a Vandermonde matrix
   - Least squares problem:
     \[
     \min_{\mathbf{c}} \|\mathbf{y} - A\mathbf{c}\|_2^2
     \]
   - Solution interpreted as orthogonal projection onto a polynomial subspace

4. **Evaluation**
   - Mean Squared Error (MSE)
   - Root Mean Squared Error (RMSE)
   - Residual \( \ell_2 \) norm
   - Visual comparison of original and approximated signals

---

## 🧪 Experimental Results

- Increasing polynomial degree leads to **monotonic error reduction**
- Error improvement shows **diminishing returns beyond degree 6**
- Low-degree polynomials capture global ECG trends
- Sharp QRS peaks remain in the residual component due to their high-frequency nature

These results align with the geometric interpretation of least squares approximation.

---

## 🛠️ Requirements

Install dependencies using:

```bash
pip install -r requirements.txt
```

### Required packages
- numpy
- matplotlib
- wfdb

---

## ▶️ Running the Code

Execute the main script:

```bash
python src/ecg_ls_approximation.py
```

The script will:
- download ECG data from PhysioNet (via WFDB)
- preprocess the signal
- compute least squares approximations for multiple polynomial degrees
- output error metrics
- generate visualization plots

---

## 📊 Output Examples

- ECG waveform vs polynomial approximation
- RMSE vs polynomial degree
- Residual analysis

These outputs are used directly in the experimental results section of the paper.

---

## ⚠️ Notes and Limitations

- This project **does not perform medical diagnosis**
- ECG annotations are not used for classification
- Polynomial approximation is global and may not capture localized sharp features
- The implementation is intended for **educational and analytical purposes**

---

## 📚 References

- MIT-BIH Arrhythmia Database, PhysioNet  
- G. Strang, *Linear Algebra and Its Applications*  
- L. Xie et al., “Computational diagnostic techniques for electrocardiogram signal analysis,” *Sensors*, 2020  

(Full references are provided in the accompanying paper.)

---

## 🎓 Academic Context

This project was developed for:
- **Linear Algebra and Geometry**
- Undergraduate (Semester 3)
- Informatics / Computer Science
