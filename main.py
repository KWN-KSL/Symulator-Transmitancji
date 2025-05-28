import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from euler import simulate_response_euler
from stability import find_roots, check_stability
from body import compute_bode_manual

class App:
    CONFIG = {
        'coeff_limit': 999,  # Maksymalna/minimalna wartość współczynników
        'amp_limit': 100,    # Maksymalna amplituda
        'freq_limit': 100,   # Maksymalna częstotliwość
        'fig_size': (6, 1.5),  # Rozmiar wykresów transmitancji
        'bode_fig_size': (8, 6),  # Rozmiar wykresów Bodego
        'font_label': ("Helvetica", 16, "bold"),  # Czcionka dla etykiet
        'font_entry': ("Helvetica", 14),  # Czcionka dla pól tekstowych
        'font_button': ("Helvetica", 18, "bold")  # Czcionka dla przycisku
    }

    def __init__(self, root):
        # Inicjalizacja aplikacji z oknem głównym Tkinter.
        self.root = root
        self.root.title("Symulator transmitancji G(s)")
        self.root.geometry("900x1050")
        self.root.configure(bg="white")
        self.entries = {}  # Słownik na pola tekstowe dla współczynników
        self.signal_type = tk.StringVar(value="Prostokątny")
        self.build_gui()

    def build_gui(self):
        # Wyświetlanie G(s) na górze okna
        self.fig_static, self.ax_static = plt.subplots(figsize=self.CONFIG['fig_size'])
        self.canvas_static = FigureCanvasTkAgg(self.fig_static, master=self.root)
        self.canvas_static.get_tk_widget().pack(pady=(10, 0))
        self.display_symbolic_transfer_function()
        # Ramka na pola współczynników
        frame = tk.Frame(self.root, bg="white")
        frame.pack(pady=20)
        # Tworzenie pól dla współczynników licznika (a3, a2, a1, a0)
        tk.Label(frame, text="Współczynniki a:", font=self.CONFIG['font_label'], bg="white").grid(row=0, column=0, columnspan=10, pady=10)
        for i, coeff in enumerate(['a3', 'a2', 'a1', 'a0']):
            tk.Label(frame, text=f"{coeff}:", font=self.CONFIG['font_label'], bg="white").grid(row=1, column=i*2+1, sticky='e', padx=10)
            entry = tk.Entry(frame, width=7, font=self.CONFIG['font_entry'])
            entry.grid(row=1, column=i*2+2, padx=10)
            self.entries[coeff] = entry
        # Tworzenie pól dla współczynników mianownika (b4, b3, b2, b1, b0)
        tk.Label(frame, text="Współczynniki b:", font=self.CONFIG['font_label'], bg="white").grid(row=2, column=0, columnspan=10, pady=20)
        for i, coeff in enumerate(['b4', 'b3', 'b2', 'b1', 'b0']):
            tk.Label(frame, text=f"{coeff}:", font=self.CONFIG['font_label'], bg="white").grid(row=3, column=i*2, sticky='e', padx=10)
            entry = tk.Entry(frame, width=7, font=self.CONFIG['font_entry'])
            entry.grid(row=3, column=i*2+1, padx=10)
            self.entries[coeff] = entry
        # Wybór typu sygnału
        tk.Label(self.root, text="Typ sygnału wejściowego:", font=self.CONFIG['font_label'], bg="white").pack(pady=20)
        self.signal_type.trace_add("write", self.update_freq_visibility)
        # Przyciski dla typów sygnałów
        signal_types = ["Prostokątny", "Trójkątny", "Harmoniczny", "Impuls jednostkowy", "Skok jednostkowy"]
        radio_frame = tk.Frame(self.root, bg="white")
        radio_frame.pack()
        for sig in signal_types:
            tk.Radiobutton(radio_frame, text=sig, variable=self.signal_type, value=sig, font=self.CONFIG['font_entry'], bg="white").pack(anchor='w')
        # Parametry sygnału (amplituda i częstotliwość)
        self.param_frame = tk.Frame(self.root, bg="white")
        self.param_frame.pack(pady=10)
        tk.Label(self.param_frame, text="Amplituda:", font=self.CONFIG['font_label'], bg="white").grid(row=0, column=0, padx=10)
        self.amplitude_entry = tk.Entry(self.param_frame, width=7, font=self.CONFIG['font_entry'])
        self.amplitude_entry.insert(0, "1.0")
        self.amplitude_entry.grid(row=0, column=1, padx=10)
        # Pole dla częstotliwości (pokazywane dla wybranych sygnałów)
        self.freq_label = tk.Label(self.param_frame, text="Częstotliwość [Hz]:", font=self.CONFIG['font_label'], bg="white")
        self.freq_label.grid(row=0, column=2, padx=10)
        self.freq_entry = tk.Entry(self.param_frame, width=7, font=self.CONFIG['font_entry'])
        self.freq_entry.insert(0, "1.0")
        self.freq_entry.grid(row=0, column=3, padx=10)
        self.update_freq_visibility()
        # Przycisk do uruchomienia symulacji
        tk.Button(self.root, text="Symuluj", font=self.CONFIG['font_button'], command=self.simulate, bg="#3e8ef7", fg="white", padx=15, pady=8).pack(pady=30)
        # Wyświetlanie dynamicznej postaci G(s) z wprowadzonymi współczynnikami
        self.fig_dynamic, self.ax_dynamic = plt.subplots(figsize=self.CONFIG['fig_size'])
        self.canvas_dynamic = FigureCanvasTkAgg(self.fig_dynamic, master=self.root)
        self.canvas_dynamic.get_tk_widget().pack(pady=(20, 30))
        self.ax_dynamic.axis('off')
        self.canvas_dynamic.draw()

    def update_freq_visibility(self, *args):
        # Pokazuje lub ukrywa pole częstotliwości w zależności od wybranego typu sygnału.
        sig = self.signal_type.get()
        needs_freq = sig in ["Prostokątny", "Trójkątny", "Harmoniczny"]
        if needs_freq:
            self.freq_label.grid()
            self.freq_entry.grid()
        else:
            self.freq_label.grid_remove()
            self.freq_entry.grid_remove()

    def simulate(self):
        plt.close('all')
        try:
            # Odczyt współczynników licznika i mianownika
            a_coeffs = [self.entries[f'a{i}'].get() for i in reversed(range(4))]
            b_coeffs = [self.entries[f'b{i}'].get() for i in reversed(range(5))]
            a = [float(v) if v.strip() else 0.0 for v in a_coeffs]
            b = [float(v) if v.strip() else 0.0 for v in b_coeffs]
            # Walidacja zakresu współczynników
            for val in a + b:
                if abs(val) > self.CONFIG['coeff_limit']:
                    messagebox.showerror("Błąd", f"Współczynniki muszą być w zakresie od -{self.CONFIG['coeff_limit']} do {self.CONFIG['coeff_limit']}.")
                    return
            # Usuwanie wiodących zer, z zachowaniem przynajmniej jednego współczynnika
            while len(a) > 1 and a[0] == 0:
                a.pop(0)
            while len(b) > 1 and b[0] == 0:
                b.pop(0)
            # Walidacja niezerowego licznika i mianownika
            if not any(a):
                messagebox.showerror("Błąd", "Licznik nie może być zerowy.")
                return
            if not any(b):
                messagebox.showerror("Błąd", "Mianownik nie może być zerowy.")
                return
            if len(a) > len(b):
                messagebox.showerror("Błąd", "Rząd licznika nie może być większy niż rząd mianownika.")
                return
            if len(a) == 1 and len(b) == 1:
                messagebox.showerror("Błąd", "Transmitancja nie może być stałą liczbą.")
                return
            # Aktualizacja dynamicznego wyświetlania transmitancji
            self.update_dynamic_transfer_function(a, b)
        except ValueError:
            messagebox.showerror("Błąd", "Wprowadź poprawne liczby dla współczynników.")
            return
        # Generowanie sygnału wejściowego
        t = np.linspace(0, 10, 5000)
        u = self.generate_input(t)
        if u is None:
            return
        # Symulacja odpowiedzi metodą Eulera – oblicza odpowiedź czasową y(t)
        y = simulate_response_euler(a, b, u, t)
        # Rysowanie odpowiedzi czasowej – przygotowuje wykres sygnału wejściowego i wyjściowego
        self.plot_time_response(t, u, y)
        # Rysowanie charakterystyki Bodego – przygotowuje wykres amplitudy i fazy
        self.plot_bode(a, b)
        # Rysowanie analizy stabilności – przygotowuje wykres biegunów
        self.plot_stability(b)
        # Wyświetlenie wszystkich wykresów jednocześnie
        plt.show()

    def generate_input(self, t):
        # Generuje sygnał wejściowy na podstawie wybranego typu i parametrów.
        try:
            A = float(self.amplitude_entry.get())
            if abs(A) > self.CONFIG['amp_limit']:
                messagebox.showerror("Błąd", f"Amplituda musi być mniejsza niż {self.CONFIG['amp_limit']}.")
                return None
        except ValueError:
            messagebox.showerror("Błąd", "Wprowadź poprawną wartość amplitudy.")
            return None
        signal_type = self.signal_type.get()
        f = None
        if signal_type in ['Prostokątny', 'Trójkątny', 'Harmoniczny']:
            try:
                f = float(self.freq_entry.get())
                if f <= 0 or f > self.CONFIG['freq_limit']:
                    messagebox.showerror("Błąd", f"Częstotliwość musi być dodatnia i mniejsza niż {self.CONFIG['freq_limit']}.")
                    return None
            except ValueError:
                messagebox.showerror("Błąd", "Wprowadź poprawną wartość częstotliwości.")
                return None
        # Generowanie sygnału na podstawie typu
        if signal_type == 'Prostokątny':
            return A * signal.square(2 * np.pi * f * t)
        elif signal_type == 'Trójkątny':
            return A * signal.sawtooth(2 * np.pi * f * t, width=0.5)
        elif signal_type == 'Harmoniczny':
            return A * np.sin(2 * np.pi * f * t)
        elif signal_type == 'Impuls jednostkowy':
            u = np.zeros_like(t)
            idx = np.argmin(np.abs(t - 0.01))
            u[idx] = A
            return u
        elif signal_type == 'Skok jednostkowy':
            return A * np.ones_like(t)
        return np.zeros_like(t)

    def coeffs_to_str(self, coeffs):
        # Konwertuje współczynniki na ciąg LaTeX (wzór transmitancji w interfejsie graficznym) do wyświetlenia G(s).
        terms = []
        power = len(coeffs) - 1
        for c in coeffs:
            if isinstance(c, str):
                term = f"{c}s^{power}" if power > 1 else (f"{c}s" if power == 1 else f"{c}")
            elif abs(c) > 1e-6:
                term = f"{c:.2f}s^{power}" if power > 1 else (f"{c:.2f}s" if power == 1 else f"{c:.2f}")
            else:
                term = None
            if term:
                terms.append(term)
            power -= 1
        return " + ".join(terms) if terms else "0"

    def display_symbolic_transfer_function(self):
        # Wyświetla G(s) z nazwami zmiennych.
        self.ax_static.clear()
        num_str = self.coeffs_to_str(['a3', 'a2', 'a1', 'a0'])
        den_str = self.coeffs_to_str(['b4', 'b3', 'b2', 'b1', 'b0'])
        latex = fr"$G(s) = \frac{{{num_str}}}{{{den_str}}}$"
        self.ax_static.text(0.5, 0.5, latex, fontsize=20, ha='center', va='center')
        self.ax_static.axis('off')
        self.canvas_static.draw()

    def update_dynamic_transfer_function(self, a, b):
        # Aktualizuje G(s) z konkretnymi wartościami współczynników.
        self.ax_dynamic.clear()
        num_str = self.coeffs_to_str(a)
        den_str = self.coeffs_to_str(b)
        latex = fr"$G(s) = \frac{{{num_str}}}{{{den_str}}}$"
        self.ax_dynamic.text(0.5, 0.5, latex, fontsize=20, ha='center', va='center', color='blue')
        self.ax_dynamic.axis('off')
        self.canvas_dynamic.draw()

    def plot_time_response(self, t, u, y):
        plt.figure("Odpowiedź czasowa")
        plt.plot(t, u, label="Wejście")
        plt.plot(t, y, label="Wyjście")
        plt.xlabel("Czas [s]")
        plt.ylabel("Amplituda")
        plt.title("Odpowiedź układu")
        plt.grid()
        plt.legend()

    def plot_bode(self, a, b):
        w = np.logspace(-2, 3, 1000)
        mag, phase = compute_bode_manual(a, b, w)
        fig_bode, (ax1, ax2) = plt.subplots(2, 1, figsize=self.CONFIG['bode_fig_size'], num="Charakterystyka Bodego")
        ax1.semilogx(w, mag)
        ax1.set_title("Bode - amplituda")
        ax1.set_ylabel("Mag [dB]")
        ax1.grid()
        ax2.semilogx(w, phase)
        ax2.set_title("Bode - faza")
        ax2.set_xlabel("Częstotliwość [rad/s]")
        ax2.set_ylabel("Faza [°]")
        ax2.grid()
        plt.tight_layout()

    def plot_stability(self, b):
        plt.figure("Stabilność (bieguny)")
        poles = find_roots(b)
        colors = []
        for p in poles:
            if np.real(p) < 0:
                colors.append('green')
            elif np.real(p) == 0:
                colors.append('orange')
            else:
                colors.append('red')
        plt.scatter(np.real(poles), np.imag(poles), color=colors, s=80)
        plt.axvline(0, color='gray', linestyle='--')
        plt.title("Bieguny układu (Re vs Im)")
        plt.xlabel("Re")
        plt.ylabel("Im")
        plt.grid()
        # Wyświetlanie komunikatu o stabilności
        print(b)
        stable = check_stability(b)
        print(stable)
        if stable == 0:
            text = "Układ jest NIESTABILNY"
            color = 'red'
        elif stable == 2:
            text = "Układ jest na GRANICY STABILNOŚCI"
            color = 'orange'
        else:
            text = "Układ jest STABILNY"
            color = 'green'
        plt.figtext(0.5, 0.01, text, ha='center', fontsize=14, color=color)
        plt.tight_layout()

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()