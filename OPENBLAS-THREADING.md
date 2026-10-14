# Fortran-Threadskalierung mit OpenBLAS

Stand: 14. Oktober 2026

## Kurzfassung

Der unveränderte Fortran-95-Code skaliert mit OpenBLAS 0.3.34 erwartungsgemäß.
Bei `n=8192` beziehungsweise einer Matrixdimension von `N=8193` sinkt die
Gesamtzeit von `22.512 s` bei einem Thread auf `3.992 s` bei 16 Threads. Das
entspricht einem Speedup von `5.64`.

Die zunächst beobachtete geringe Skalierung war auf den `DGESV`-Treiber des
auf `turing` installierten Systempakets OpenBLAS 0.3.29 zurückzuführen.

## Ursache in OpenBLAS 0.3.29

`DGESV` löst ein lineares Gleichungssystem im Wesentlichen durch eine
LU-Zerlegung mit `DGETRF` und das anschließende Rückwärtseinsetzen mit
`DGETRS`. OpenBLAS 0.3.29 entschied über den parallelen Pfad anhand von

```text
N * NRHS
```

statt anhand der quadratischen Matrixgröße. Für den Benchmark gilt:

```text
N = 8193
NRHS = 1
N * NRHS = 8193 < 10000
```

Dadurch wählte `DGESV` trotz korrekt angeforderter Threads den seriellen
Pfad. Der direkte Aufruf von `DGETRF` skalierte dagegen bereits in OpenBLAS
0.3.29. Ein zusätzlicher Grenztest machte den Fehler sichtbar:

| Aufruf mit 20 Threads | Wandzeit | CPU-/Wandzeit |
|---|---:|---:|
| `DGESV`, N=9999 | 8.610 s | 1.00 |
| `DGESV`, N=10000 | 1.540 s | 19.89 |

Die Upstream-Korrektur ersetzt die Prüfung durch die quadratische
Matrixgröße und ist ab OpenBLAS 0.3.32 enthalten:
[OpenBLAS-Commit dc8b16c](https://github.com/OpenMathLib/OpenBLAS/commit/dc8b16c57cffcb28b53ea3749d99e52d2b97f956).

## Vergleich der vollständigen Messreihen

Gemessen wurde der dichte Arbeitsablauf aus Matrixaufbau, direkter Lösung,
Aufbau des Tikhonov-Systems und regularisierter Lösung. Jeder Tabellenwert ist
Median aus drei unabhängigen Prozessen; jede Phase verwendet innerhalb eines
Prozesses einen Warm-up und den Median aus drei Wiederholungen. Bis 20 Threads
liegen alle Kerne auf dem ersten CPU-Sockel, 40 Threads verwenden beide
Sockel. Hyperthreading bleibt ausgeschlossen.

| Threads | OpenBLAS 0.3.29 [s] | OpenBLAS 0.3.34 [s] | Speedup 0.3.34 |
|---:|---:|---:|---:|
| 1 | 22.705 | 22.512 | 1.00 |
| 2 | 16.700 | 11.849 | 1.90 |
| 4 | 13.802 | 6.440 | 3.50 |
| 8 | 12.453 | 4.210 | 5.35 |
| 16 | 12.082 | 3.992 | 5.64 |
| 20 | 12.059 | 4.130 | 5.45 |
| 40 | 12.271 | 7.027 | 3.20 |

Bei 16 Threads beschleunigt sich der `DGEMM`-Anteil um den Faktor `7.71`.
Die direkte und die regularisierte `DGESV`-Lösung erreichen Speedups von
`4.53` beziehungsweise `4.76`. Der explizite Matrixaufbau bleibt
einthreadig, benötigt aber nur rund `0.23 s`.

Der 40-Thread-Punkt überschreitet die NUMA-Grenze zwischen den beiden
Prozessorsockeln und ist deshalb langsamer als die Messpunkte mit 16 und
20 Threads.

## Verwendeter OpenBLAS-Build

OpenBLAS 0.3.34 wurde auf `turing` ohne Administratorrechte gebaut:

```bash
git clone --branch v0.3.34 --depth 1 \
    https://github.com/OpenMathLib/OpenBLAS.git OpenBLAS-0.3.34
make -C OpenBLAS-0.3.34 -j20 \
    DYNAMIC_ARCH=1 USE_THREAD=1 NO_AFFINITY=1 NUM_THREADS=80
make -C OpenBLAS-0.3.34 \
    PREFIX="$HOME/opt/openblas-0.3.34" install
```

Der Fortran-Build kann anschließend explizit gegen diese Bibliothek gelinkt
werden:

```bash
cmake -S fortran -B fortran/build-openblas \
    -DCMAKE_BUILD_TYPE=Release \
    -DGWR_BLAS_LAPACK_LIBRARY="$HOME/opt/openblas-0.3.34/lib/libopenblas.so"
cmake --build fortran/build-openblas --parallel
ctest --test-dir fortran/build-openblas --output-on-failure
```

Beispiel für den 16-Thread-Benchmark:

```bash
OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 \
taskset -c 0-15 fortran/build-openblas/gwr_thread_benchmark_16
```

Die Fortran-Komponententests und der unabhängige Vergleich mit Python wurden
mit OpenBLAS 0.3.34 erfolgreich ausgeführt. Der Wechsel betrifft nur das
numerische Backend; der Fortran-Quellcode blieb unverändert.
