# Software

## Python

Für numerische Berechungen benutzen wir im Kurs Python. Eine gute Einführung in Python ist im interaktiven Jupyter-Book [Python Programming And Numerical Methods: A Guide For Engineers And Scientists](https://pythonnumericalmethods.berkeley.edu/notebooks/Index.html) zu finden.

Im wesentlichen werden wir fundamentale Funktionalität von Python benutzen. Für numerische Berechnungen werden wir auf NumPy und SciPy zurückgreifen. Die Visualisierung kann sehr matlab nahe mit Hilfe von Matplotlib umgesetzt werden.

- [Python](https://python.org)
- [NumPy](https://numpy.org)
- [SciPy](https://scipy.org)
- [Matplotlib](https://matplotlib.org)

Falls Sie Python mit den entsprechenden Module noch nie benutzt haben, ist es sehr sinnvoll sich vorab mit dem wesentlichen auseinander zu setzen.


## NGSolve

Im Rahmen der numerischen Methoden für partielle Differentialgleichungen werden wir die Methode der finiten Elemente kennen lernen. Wir werden für die Umsetzung auf eine C++ Bibliothek [NGSolve](https://ngsolve.org) zurückgreifen, welche eine umfangreiche Python-Schnittstelle zur Verfügung stellt.

Netgen/NGSolve ist eine hochleistungsfähige Multiphysik-Finite-Elemente-Software. Sie wird häufig zur Analyse von Modellen aus den Bereichen Festkörpermechanik, Strömungsmechanik und Elektromagnetik eingesetzt. Dank der flexiblen Python-Schnittstelle können neue physikalische Gleichungen und Lösungsalgorithmen leicht implementiert werden.

Das Ziel des Moduls besteht in der Vermittlung der mathematischen Grundlagen und insbesondere auch numerischen Anwendung derer auf konkrete Beispiele aus der Ingenieur Praxis.

```{note}
Installieren Sie NGSolve auf ihrem Rechner. Die dazu notwendigen Anleitungen finden Sie auf (https://ngsolve.org) unter dem Abschnitt INSTALLATION:

- [Downloads NGSolve](https://ngsolve.org/downloads)
- [Using Jupyter notebook](https://docu.ngsolve.org/latest/install/usejupyter.html)
```

## Jupyter Hub

Auf dem Rechner `CLT-DSK-T-7307.zhaw.ch` ist ein JupyterHub eingerichtet. Der Rechner ist nur im Intranet verfügbar (VPN Verbindung notwendig). 

[JupyterHub](http://CLT-DSK-T-7307.zhaw.ch)

Ein Account kann selbständig kreiert werden, wobei dieser authentifiziert wird. Bitte benutzen Sie als User-Name Ihr ZHAW Kürzel und ein sicheres Passwort.

```{figure} native_auth_flow.png
---
align: left
height: 250px
name: native_auth_flow
---
Authetication workflow
```