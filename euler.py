import numpy as np

def simulate_response_euler(a, b, u, t):
    # Obliczenie kroku czasowego (odstęp między próbkami)
    Ts = t[1] - t[0]
    # Liczba próbek czasowych
    N = len(t)
    # Uzupełnienie współczynników zerami (aby listy miały stałą długość)
    a = [0.0] * (4 - len(a)) + a
    b = [0.0] * (5 - len(b)) + b
    # Odwrócenie kolejności współczynników dla zgodności z równaniem różniczkowym
    a = a[::-1]
    b = b[::-1]
    # Normalizacja przez b0
    if b[0] > 1e-8:
        a = [coef / b[0] for coef in a]
        b = [coef / b[0] for coef in b]
    a0, a1, a2, a3 = a
    b0, b1, b2, b3, b4 = b
    y = np.zeros(N)
    dy = np.zeros(N)
    d2y = np.zeros(N)
    d3y = np.zeros(N)
    d4y = np.zeros(N)
    du = np.zeros(N)
    d2u = np.zeros(N)
    d3u = np.zeros(N)
    # Określenie rzędu licznika i mianownika (najwyższy niezerowy współczynnik)
    a_order = max([i for i in range(4) if abs(a[i]) > 1e-8], default=-1)
    b_order = max([i for i in range(5) if abs(b[i]) > 1e-8], default=-1)
    system_order = max(a_order, b_order)
    if system_order == -1:
        raise ValueError("Mianownik i licznik transmitancji są zerowe – układ nieprawidłowy.")
    # Obliczenie pochodnych sygnału wejściowego
    if a1 != 0.0 or a2 != 0.0 or a3 != 0.0:
        for k in range(1, N):
            du[k] = (u[k] - u[k-1]) / Ts
    if a2 != 0.0 or a3 != 0.0:
        for k in range(2, N):
            d2u[k] = (du[k] - du[k-1]) / Ts
    if a3 != 0.0:
        for k in range(3, N):
            d3u[k] = (d2u[k] - d2u[k-1]) / Ts
    # Symulacja odpowiedzi w zależności od rzędu układu
    if system_order == 4:
        for k in range(4, N):
            d4y[k] = (a3 * d3u[k] + a2 * d2u[k] + a1 * du[k] + a0 * u[k] -
                    b3 * d3y[k-1] - b2 * d2y[k-1] - b1 * dy[k-1] - b0 * y[k-1]) / b4
            d3y[k] = d3y[k-1] + Ts * d4y[k]
            d2y[k] = d2y[k-1] + Ts * d3y[k]
            dy[k] = dy[k-1] + Ts * d2y[k]
            y[k] = y[k-1] + Ts * dy[k]
    elif system_order == 3:
        for k in range(3, N):
            d3y[k] = (a3 * d3u[k] + a2 * d2u[k] + a1 * du[k] + a0 * u[k] -
                    b2 * d2y[k-1] - b1 * dy[k-1] - b0 * y[k-1]) / b3
            d2y[k] = d2y[k-1] + Ts * d3y[k]
            dy[k] = dy[k-1] + Ts * d2y[k]
            y[k] = y[k-1] + Ts * dy[k]
    elif system_order == 2:
        for k in range(2, N):
            d2y[k] = (a2 * d2u[k] + a1 * du[k] + a0 * u[k] -
                    b1 * dy[k-1] - b0 * y[k-1]) / b2
            dy[k] = dy[k-1] + Ts * d2y[k]
            y[k] = y[k-1] + Ts * dy[k]
    elif system_order == 1:
        for k in range(1, N):
            dy[k] = (a1 * du[k] + a0 * u[k] - b0 * y[k-1]) / b1
            y[k] = y[k-1] + Ts * dy[k]
    elif system_order == 0:
        for k in range(N):
            y[k] = a0 * u[k] / b0
    return y