import numpy as np

def compute_bode_manual(num, den, w_range):
    G_mag = []
    G_phase = []

    for w in w_range:
        jw = 1j * w

        # Zakładamy, że num/den są od najwyższej potęgi do stałej (jak w MATLAB)
        num_eval = sum(c * jw ** (len(num) - i - 1) for i, c in enumerate(num))
        den_eval = sum(c * jw ** (len(den) - i - 1) for i, c in enumerate(den))

        if abs(den_eval) < 1e-12:
            G_mag.append(np.inf)
            G_phase.append(np.nan)
        else:
            G = num_eval / den_eval
            G_mag.append(20 * np.log10(abs(G)))
            G_phase.append(np.angle(G, deg=True))

    # Odwijanie fazy po pętli
    G_phase = -np.unwrap(np.radians(G_phase)) * (180 / np.pi) #minus na poczatku zeby odwrocic faze

    return G_mag, G_phase
