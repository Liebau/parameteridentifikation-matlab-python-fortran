# Fortran-90/95-Referenz

Diese Implementierung reproduziert die dichte diskrete
Parameteridentifikation aus den vorhandenen MATLAB- und Python-Fassungen. Sie
ist eine neu entwickelte numerische Referenz desselben Modells und keine
Migration eines historischen Fortran-Programms.

## Modell und Varianten

Die Stuetzstellen sind `x_i = i/n`, die Randgewichte `1/(2n)` und alle
inneren Trapezgewichte `1/n`. Damit gilt

```text
A(i,j) = w_j exp(x_i x_j),
y(t)   = (exp(t+1)-1)/(t+1),
lambda(t) = exp(t).
```

Die gewichtete, zum Modell konsistente Regularisierung wird im Quellcode
ausdruecklich als

```text
(beta I + A*A) alpha = A*y
```

aufgebaut. Zum dokumentierten Vergleich gibt es separat die euklidische
Variante

```text
(beta I + transpose(A)*A) alpha = transpose(A)*y.
```

## Build und Tests

Vorausgesetzt werden ein Fortran-Compiler sowie BLAS und LAPACK. Getestet ist
der Build mit `gfortran` und CMake:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

CTest prueft Matrix und Trapezgewichte fuer `n=2`, den stabilen direkten Fall
`n=4`, die gewichtete Regularisierung fuer `n=32`, `beta=1e-6`, den deutlichen
Unterschied der beiden Adjungierten, wiederholten Aufbau und Abbau der
Pointer-basierten Problemdaten sowie alle ausgegebenen stabilen
Rekonstruktionsfehler gegen Python.

Die Referenztabelle und CSV-Datei entstehen mit:

```bash
(cd results && ../build/gwr_reference)
```

Der Python-Vergleich kann auch separat ausgefuehrt werden:

```bash
python3 tools/compare_with_python.py \
    --fortran-executable build/gwr_reference \
    --python-source ../python/src \
    --output build/reference_fortran.csv
```

## Dichte Benchmarks

`gwr_benchmark` misst dieselben vier Phasen wie die vorhandenen MATLAB- und
Python-Benchmarks: Matrixaufbau, direkte Loesung, Aufbau des konsistenten
Tikhonov-Systems und regularisierte Loesung. Jede Phase erhaelt einen Warm-up;
aus neun Wiederholungen wird der Median gebildet:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    taskset -c 0 ./build/gwr_benchmark
```

Die sieben Programme `gwr_thread_benchmark_<threads>` messen `n=8192` mit
drei internen Wiederholungen. Threadzahl und CPU-Affinitaet werden wie bei
MATLAB und Python ausserhalb des Fortran-Codes gesteuert. Nur das Zeitmodul
verwendet fuer einen hochaufgeloesten 64-Bit-Wandzeitzaehler die von gfortran
unterstuetzte neuere `SYSTEM_CLOCK`-Variante. Es liegt deshalb in der
separaten Bibliothek `gwr_benchmark_support`; die numerische Bibliothek
`gwr_parameteridentifikation` bleibt mit strengem Fortran 95 uebersetzt.

Fuer reproduzierbare Messungen kann ein bestimmtes kombiniertes
BLAS-/LAPACK-Backend beim Konfigurieren explizit gewaehlt werden. Auf
`turing` wurde OpenBLAS 0.3.34 benutzerlokal installiert und so eingebunden:

```bash
cmake -S . -B build-openblas -DCMAKE_BUILD_TYPE=Release \
    -DGWR_BLAS_LAPACK_LIBRARY=/path/to/openblas-0.3.34/lib/libopenblas.so
cmake --build build-openblas --parallel
```

Das auf dem Testsystem vorhandene OpenBLAS 0.3.29 besitzt im `DGESV`-Treiber
einen inzwischen korrigierten Schwellwertfehler, durch den der hier verwendete
Fall mit einer rechten Seite seriell ausgefuehrt wird. Fuer Threadvergleiche
sollte deshalb OpenBLAS 0.3.32 oder neuer verwendet werden. Die genaue
Diagnose und die Messreihe mit OpenBLAS 0.3.34 sind in
[`../OPENBLAS-THREADING.md`](../OPENBLAS-THREADING.md) dokumentiert.

## Lineare Algebra

- `DGESV` loest direkte und regularisierte lineare Gleichungssysteme.
- `DGESVD` liefert die Singulaerwerte fuer die 2-Norm-Konditionszahl.
- `DGEMM` baut die beiden unterschiedlichen Tikhonov-Matrizen auf.
- `DGEMV` bildet die zugehoerigen rechten Seiten.

Eingabematrizen fuer `DGESV` und `DGESVD` werden kopiert, da LAPACK sie
ueberschreibt. Alle LAPACK-Statuscodes werden geprueft. BLAS-Routinen besitzen
keinen Rueckgabestatus; Dimensionen und Modus werden vor ihrem Aufruf
validiert.

## Toleranzen und Grenzen

Fuer Fehler und relative Fehler verwendet der Sprachvergleich `rtol=5e-7`
und `atol=2e-9`. Das laesst kleine Unterschiede durch Rechenreihenfolge und
BLAS/LAPACK-Backend zu, liegt aber deutlich unter den ausgewiesenen
Rekonstruktionsfehlern. Konditionszahlen werden nur fuer die stabilen direkten
Faelle `n=1,2,4` verglichen. Bei feineren Gittern liegen kleinste Singulaerwerte
an oder unter der Rundungsgrenze und sind backendabhaengig.

Die Niedrigrangapproximation und eine spaetere Java-Migration gehoeren nicht
zu dieser Referenzimplementierung. Dichte Performance- und Threadmessungen
sind als getrennte Benchmark-Werkzeuge enthalten; sie aendern den
Referenzalgorithmus nicht.

## Zuordnung zur Python-Referenz

Die Kommentare in den Fortran-Quellen verweisen bewusst auf die fachlich
entsprechenden Python-Funktionen. Die Zuordnung lautet:

| Fortran | Python |
|---|---|
| `src/precision.f90` | `numpy.float64` |
| `src/model.f90: build_problem` | `model.py: build_problem` |
| `src/linear_algebra.f90: solve_linear_system` | `numpy.linalg.solve` |
| `src/linear_algebra.f90: condition_number_2` | `numpy.linalg.svd` |
| `src/linear_algebra.f90: build_tikhonov_system` | regularisierter Zweig von `model.py: analyze` |
| `src/evaluation.f90: analyze` | `model.py: analyze` |
| `src/run_reference.f90` | `run_reference.py` |
| `tests/test_reference.f90` | die dichten Referenztests in `tests/test_model.py` |

Die Verweise dokumentieren mathematische Aequivalenz, nicht eine automatische
oder zeilenweise Uebersetzung von Python nach Fortran.


