# Parameteridentifikation mit MATLAB, Python und Fortran 90/95

[![Tests](https://github.com/Liebau/parameteridentifikation-matlab-python-fortran/actions/workflows/ci.yml/badge.svg)](https://github.com/Liebau/parameteridentifikation-matlab-python-fortran/actions/workflows/ci.yml)

Dieses Repository enthält drei numerisch äquivalente Implementierungen eines
Modellproblems zur Parameteridentifikation. Es begleitet den GWR-Fachbeitrag
[Parameteridentifikation](https://gwr-mbh.de/parameteridentifikation/).

Im Mittelpunkt stehen nicht syntaktische Unterschiede der Sprachen, sondern
eine konsistente Diskretisierung, Tikhonov-Regularisierung und reproduzierbare
Vergleiche der numerischen Ergebnisse.

## Mathematisches Modell

Auf den äquidistanten Stützstellen `x_i = i/n` wird der Integraloperator mit
der summierten Trapezregel diskretisiert:

```text
A(i,j) = w_j exp(x_i x_j)
y(t)   = (exp(t+1)-1)/(t+1)
lambda(t) = exp(t)
```

Die direkte Aufgabe lautet

```text
A alpha = y.
```

Die zum gewichteten diskreten Modell konsistente Regularisierung wird als

```text
(beta I + A A) alpha_beta = A y
```

gelöst. Daneben bleibt die euklidische Variante mit `transpose(A)` als
ausdrücklich gewählter Vergleich verfügbar.

## Implementierungen

| Verzeichnis | Inhalt |
|---|---|
| [`matlab/`](matlab/) | MATLAB-Implementierung, Tests, Benchmarks und Niedrigrang-Erweiterung |
| [`python/`](python/) | Python-/NumPy-Paket, Tests, Benchmarks und Niedrigrang-Erweiterung |
| [`fortran/`](fortran/) | modulare Fortran-90/95-Referenz mit BLAS/LAPACK, CMake und CTest |

Die dichte Referenzrechnung ist in allen drei Sprachen vorhanden. Die
separierbare Niedrigrangformulierung ist eine weiterführende MATLAB-/Python-
Ergänzung; sie verändert die gemeinsame dichte Referenz nicht.

## Python ausführen

Vorausgesetzt wird Python 3.13.

```bash
cd python
python -m venv .venv
. .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
python -m gwr_parameteridentifikation.run_reference
```

## Fortran bauen und prüfen

Benötigt werden ein Fortran-Compiler, CMake sowie BLAS und LAPACK. Getestet
ist die Implementierung mit GNU Fortran und OpenBLAS.

```bash
cmake -S fortran -B fortran/build -DCMAKE_BUILD_TYPE=Release
cmake --build fortran/build --parallel
ctest --test-dir fortran/build --output-on-failure
```

Die Fortran-Tests prüfen unter anderem den Modellaufbau, stabile direkte
Fälle, beide Adjungierten-Konventionen und den Speicher-Lebenszyklus. Ein
zusätzlicher CTest vergleicht die Fortran-Ausgabe mit der unabhängigen
Python-Implementierung.

## MATLAB ausführen

In MATLAB in das Verzeichnis `matlab` wechseln und ausführen:

```matlab
run_tests
run_reference
```

MATLAB ist nicht Bestandteil der automatisierten GitHub-CI. Die gleichen
fachlichen Referenzfälle werden dort durch Python und Fortran geprüft.

## Sprachübergreifende Validierung

Für 29 stabile direkte und regularisierte Ergebniswerte beträgt die größte
beobachtete absolute Fehlerdifferenz zwischen Fortran und Python
`1.310e-9`, die größte relative Differenz `4.818e-10`. Feinere direkte
Systeme sind extrem schlecht konditioniert und werden deshalb nicht als
sprachübergreifende Genauigkeitsreferenz interpretiert.

## Fortran-Threadskalierung

Auf `turing` erreicht der vollständige dichte Fortran-Arbeitsablauf bei
`n=8192` mit OpenBLAS 0.3.34 sein Minimum bei 16 Threads: Die Gesamtzeit
sinkt von `22.512 s` auf `3.992 s`, entsprechend einem Speedup von `5.64`.
Damit liegt Fortran in dieser Messung praktisch auf dem Niveau von
MATLAB/MKL und Python/NumPy.

Eine frühere Messung mit dem AlmaLinux-Systempaket OpenBLAS 0.3.29 skalierte
in `DGESV` nicht. Ursache war ein inzwischen korrigierter Schwellwertfehler
dieser OpenBLAS-Version und keine Grenze von Fortran 95. Messwerte,
Kontrollversuche und Reproduktionshinweise stehen in
[OPENBLAS-THREADING.md](OPENBLAS-THREADING.md).

## Lizenz

Der Quellcode steht unter der [MIT-Lizenz](LICENSE).

