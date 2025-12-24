import numpy as np
import matplotlib.pyplot as plt


def movingAverage(x: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return x.copy()
    window = int(window)
    kernel = np.ones(window) / window
    x_pad = np.pad(x, (window//2, window-1-window//2), mode='reflect')
    return np.convolve(x_pad, kernel, mode='valid')

def zscoreNormalize(x: np.ndarray) -> np.ndarray:
    mu = np.mean(x) 
    sigma = np.std(x)
    if sigma < 1e-12:
        return x - mu
    return (x - mu) / sigma

def preprocessEcg(seg: np.ndarray, smooth_window: int = 9, normalize: bool = True) -> np.ndarray:
    x = seg.astype(float)
    x = x - np.mean(x)                
    x = movingAverage(x, smooth_window)  
    if normalize:
        x = zscoreNormalize(x)
    return x


def vandermonde(t: np.ndarray, degree: int) -> np.ndarray:
    return np.vander(t, N=degree + 1, increasing=True)

def polynomialLeastsquares(y: np.ndarray, degree: int, solver: str = "lstsq"):
    n = len(y)
    # normalized time in [-1,1] improves numerical stability for polynomial fitting
    t = np.linspace(-1.0, 1.0, n)
    A = vandermonde(t, degree)

    if solver == "normal_eq":
        # Normal equations: (A^T A) c = A^T y
        ATA = A.T @ A
        ATy = A.T @ y
        # Use solve (not inverse)
        coeff = np.linalg.solve(ATA, ATy)
    else:
        # Stable least squares solver (QR/SVD internally)
        coeff, _, _, _ = np.linalg.lstsq(A, y, rcond=None)

    y_hat = A @ coeff
    return y_hat, coeff, A, t



def calcMetrics(y: np.ndarray, y_hat: np.ndarray):
    """Return MSE, RMSE, residual norm."""
    residual = y - y_hat
    mse = np.mean(residual**2)
    rmse = np.sqrt(mse)
    rnorm = np.linalg.norm(residual, ord=2)
    return mse, rmse, rnorm



def extractSegment(signal: np.ndarray, fs: int, start_sec: float, duration_sec: float) -> np.ndarray:
    startidx = int(round(start_sec * fs))
    endidx = int(round((start_sec + duration_sec) * fs))
    endidx = min(endidx, len(signal))
    return signal[startidx:endidx]


def loadMitbih(record_name: str = "101", lead_index: int = 0):
    """
    Load MIT-BIH record using wfdb.
    Requires: pip install wfdb
    """
    import wfdb
    rec = wfdb.rdrecord(record_name, pn_dir="mitdb")
    fs = int(rec.fs)
    # p_signal is physical units (mV-like)
    sig = rec.p_signal[:, lead_index]
    lead_name = rec.sig_name[lead_index] if hasattr(rec, "sig_name") else f"lead_{lead_index}"
    return sig, fs, lead_name

def main():
    signal, fs, lead_name = loadMitbih("101", lead_index=0)  

    start_sec = 60.0     
    duration_sec = 4.0   
    raw_seg = extractSegment(signal, fs, start_sec, duration_sec)

    seg = preprocessEcg(raw_seg, smooth_window=9, normalize=True)

    degrees = [2, 4, 6, 8, 10]
    results = []
    fits = {}

    for d in degrees:
        y_hat, coeff, A, t = polynomialLeastsquares(seg, degree=d, solver="lstsq")
        mse, rmse, rnorm = calcMetrics(seg, y_hat)
        results.append((d, mse, rmse, rnorm))
        fits[d] = (y_hat, coeff)

    print(f"MIT-BIH Record 101 | Lead: {lead_name} | fs={fs} Hz")
    print(f"Segment: start={start_sec}s, duration={duration_sec}s, samples={len(seg)}")
    print("Degree |     MSE     |    RMSE    | Residual L2 Norm")
    print("-----------------------------------------------------")
    for d, mse, rmse, rnorm in results:
        print(f"{d:>6} | {mse:>10.6f} | {rmse:>10.6f} | {rnorm:>16.6f}")


    chosen_degree = 6
    y_hat, coeff = fits[chosen_degree]

    x_axis = np.arange(len(seg)) / fs 

    plt.figure()
    plt.plot(x_axis, seg, label="Preprocessed ECG segment")
    plt.plot(x_axis, y_hat, label=f"LS polynomial approx. (degree={chosen_degree})")
    plt.xlabel("Time (s)")
    plt.ylabel("Normalized amplitude")
    plt.title("ECG Waveform vs Least Squares Approximation")
    plt.legend()
    plt.show()

    degs = [r[0] for r in results]
    rmses = [r[2] for r in results]

    plt.figure()
    plt.plot(degs, rmses, marker="o")
    plt.xlabel("Polynomial degree")
    plt.ylabel("RMSE")
    plt.title("Effect of Polynomial Degree on Approximation Error")
    plt.show()

    coeff_l2 = np.linalg.norm(coeff, ord=2)
    mse, rmse, rnorm = calcMetrics(seg, y_hat)

    print("\nFeature extraction :")
    print(f"- Chosen degree: {chosen_degree}")
    print(f"- RMSE (fit error): {rmse:.6f}")
    print(f"- Residual L2 norm: {rnorm:.6f}")
    print(f"- Coefficient L2 norm: {coeff_l2:.6f}")

if __name__ == "__main__":
    main()
