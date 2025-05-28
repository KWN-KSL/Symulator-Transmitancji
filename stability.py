import numpy as np
import cmath as cm
   
def cbrt(z):
    """Pierwiastek sześcienny liczby zespolonej"""
    return z**(1/3) if z.imag != 0 else (z.real)**(1/3) if z.real >= 0 else -(-z.real)**(1/3)

def round_complex(z, tol=1e-6):
    """Zaokrągla bardzo małe części rzeczywiste i urojone do zera"""
    real = 0 if abs(z.real) < tol else z.real
    imag = 0 if abs(z.imag) < tol else z.imag
    return complex(real, imag)

def unique_roots(roots, tol=1e-6):
    """Usuwa duplikaty pierwiastków z uwzględnieniem tolerancji numerycznej"""
    unique = []
    for r in roots:
        if not any(abs(r - u) < tol for u in unique):
            unique.append(r)
    return unique

def find_roots(coeffs):
    degree = len(coeffs) - 1

    if degree == 0:
        c = coeffs[0]
        if c == 0:
            return ["Nieskończenie wiele rozwiązań"]
        else:
            return ["Brak pierwiastków (równanie stałe != 0)"]

    elif degree == 1:
        a, b = coeffs
        if a == 0:
            return ["Brak pierwiastków (a=0)"]
        root = -b / a
        return [round_complex(root)]

    elif degree == 2:
        a, b, c = coeffs
        if a == 0:
            return find_roots([b, c])
        delta = cm.sqrt(b**2 - 4*a*c)
        x1 = (-b + delta) / (2*a)
        x2 = (-b - delta) / (2*a)
        roots = [x1, x2]

    elif degree == 3:
        a, b, c, d = coeffs
        if a == 0:
            return find_roots([b, c, d])

        p = (3*a*c - b**2) / (3*a**2)
        q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)
        delta = (q**2) / 4 + (p**3) / 27

        shift = -b / (3*a)

        if abs(delta) < 1e-12:
            if abs(q) < 1e-12:
                t = 0
                roots = [shift + t] * 3
            else:
                u = cbrt(-q / 2)
                t1 = 2*u
                t2 = -u
                roots = [shift + t1, shift + t2, shift + t2]
        elif delta > 0:
            u = cbrt(-q/2 + cm.sqrt(delta))
            v = cbrt(-q/2 - cm.sqrt(delta))
            t1 = u + v
            t2 = -(u + v)/2 + (u - v)*cm.sqrt(3)*1j/2
            t3 = -(u + v)/2 - (u - v)*cm.sqrt(3)*1j/2
            roots = [shift + t1, shift + t2, shift + t3]
        else:
            r = cm.sqrt(-p**3 / 27)
            phi = cm.phase(-q / 2 + r)
            m = 2 * cm.sqrt(-p / 3)
            t1 = m * cm.cos(phi / 3)
            t2 = m * cm.cos((phi + 2*cm.pi) / 3)
            t3 = m * cm.cos((phi + 4*cm.pi) / 3)
            roots = [shift + t1, shift + t2, shift + t3]

    elif degree == 4:
        roots = np.roots(coeffs).tolist()

    else:
        return ["Obsługiwane są tylko stopnie 0–4"]

    # Zaokrąglanie i redukcja duplikatów
    roots = [round_complex(r) for r in roots]
    return unique_roots(roots)
    
def check_stability(coeffs):
    """
    Sprawdza stabilność układu metodą Routha–Hurwitza.
    Zwraca:
    1 - stabilny
    0 - niestabilny
    2 - na granicy stabilności
    Dodatkowo wyświetla tablicę Routha.
    """
    n = len(coeffs)
    if coeffs[0] == 0:
        raise ValueError("Wiodący współczynnik nie może być zerem.")

    # Tworzenie pustej tablicy Routha
    columns = (n + 1) // 2
    routh = np.zeros((n, columns))

    # Pierwsze dwa wiersze
    routh[0, :len(coeffs[::2])] = coeffs[::2]
    routh[1, :len(coeffs[1::2])] = coeffs[1::2]

    marginal = False

    # Wypełnianie tablicy
    for i in range(2, n):
        if np.all(np.abs(routh[i - 1]) < 1e-6):
            # Wiersz zerowy -> granica stabilności
            marginal = True
            # Zostawiamy zerowy wiersz, nie modyfikujemy
            continue

        for j in range(columns - 1):
            a = routh[i - 2][0]
            b = routh[i - 2][j + 1]
            c = routh[i - 1][0]
            d = routh[i - 1][j + 1]

            if abs(c) < 1e-6:
                c = 1e-6  # unikamy dzielenia przez zero

            routh[i][j] = (c * b - a * d) / c

    # WYŚWIETLENIE TABLICY ROUTHA
    print("\nTablica Routha:")
    for i in range(n):
        print(f"{'s^' + str(n - i - 1):>5}: ", end='')
        for val in routh[i]:
            print(f"{val:>10.4f}", end=' ')
        print()

    # Sprawdzenie zmian znaków w pierwszej kolumnie
    first_col = [row[0] for row in routh if abs(row[0]) > 1e-6]
    sign_changes = sum(1 for i in range(1, len(first_col)) if first_col[i] * first_col[i - 1] < 0)

    if sign_changes > 0:
        return 0  # niestabilny
    elif marginal:
        return 2  # na granicy stabilności
    else:
        return 1  # stabilny