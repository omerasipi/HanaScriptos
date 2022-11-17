---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Eindimensionaler Fall

Im eindimensionalen Fall ist eine Zerlegung des Rechengebietes gegeben durch die Teilintervalle $T_i = [x_i,x_{i+1}]$, wobei die Punkte

$$a=x_0 < x_1 < \ldots < x_{n-1} < x_n = b$$

gegeben sind. Es gilt daher 

$$[a,b] = \bigcup_{i=1}^n T_i.$$

Wir betrachten nun die Matrix $A$ {eq}`eq:linBilinForm1dproblem`. Die praktische Berechnung führt nicht über alle Kombinationen $j, k = 1, \ldots, n$. Die Basis Funktionen sind abhängig von der Zerlegung und können daher nicht vorab definiert werden. Ziel des Assembling ist, die Integration einmal mit Hilfe eines Referenz-Intervalls / Referenz-Elements durchzuführen und auf diese Berechnung zurück zu greifen.

Dazu betrachtet man das Integral in jedem Teilintervall und summiert diese auf

$$
A_{j,k} = \underbrace{\int_a^b\varphi_j'(x)\cdot \varphi_k'(x) dx}_{\text{global}} = \sum_{i=1}^n \underbrace{\int_{T_i} \varphi_j'(x)\cdot \varphi_k'(x) dx}_{\text{lokal (pro Element)}}.$$

Mit der Transformation

$$\begin{array}{rl}\sigma_i : [0,1] & \to  T_i\\
t & \mapsto \sigma_i(t) = x_i + (x_{i+1}-x_i)\cdot t\end{array}\quad\text{bzw.}\quad
\begin{array}{rl}\sigma_i^{-1} : T_i & \to [0,1]\\
x & \mapsto \sigma_i^{-1}(x) = \frac{x-x_i}{x_{i+1}-x_i} \quad\end{array} 
$$

folgt für das Integral

$$
\int_{T_i} \varphi_j'(x)\cdot \varphi_k'(x)\, dx = \int_{x_i}^{x_{i+1}} \varphi_j'(x)\cdot \varphi_k'(x)\, dx  = \int_{\sigma_i^{-1}(x_i)}^{\sigma_i^{-1}(x_{i+1})} \varphi_j'(\sigma_i(t))\cdot \varphi_k'(\sigma_i(t))\cdot \underbrace{\dot{\sigma}_i(t)\cdot dt}_{= dx}.
$$ (eq:femlokalglobal)

Sei nun $\tilde{\varphi}(t)$ die auf dem Einheitsintervall $[0,1]$ definierte Basisfunktion, dann gilt

$$\tilde{\varphi}(t) = \varphi(\sigma_i(t)).$$ (eq:femlokalbasisfkt)

Für die Ableitung folgt somit

$$\frac{d}{dt}\tilde{\varphi}(t) = \varphi'(\sigma_i(t))\cdot \dot{\sigma}_i(t)\quad\Rightarrow\quad 
\varphi'(\sigma_i(t)) = \frac{\dot{\tilde{\varphi}}(t)}{\dot{\sigma}_i(t)}.$$

Schreibt man das Integral {eq}`eq:femlokalglobal` mit den Basisfunktionen auf dem Einheitsintervall, so folgt

$$
\int_{T_i} \varphi_j'(x)\cdot \varphi_k'(x)\, dx = \int_{\sigma_i^{-1}(x_i)}^{\sigma_i^{-1}(x_{i+1})} \frac{\dot{\tilde{\varphi}}_j(t)}{\dot{\sigma}_i(t)} \cdot \frac{\dot{\tilde{\varphi}}_k(t)}{\dot{\sigma}_i(t)}\cdot \dot{\sigma}_i(t) dt = \int_{0}^{1} \dot{\tilde{\varphi}}_j(t) \cdot \dot{\tilde{\varphi}}_k(t) \frac{1}{\dot{\sigma}_i(t)} dt.$$

Da die Transformation eine affine Funktion ist ($\sigma, \sigma^{-1}$ sind Polynome ersten Grades), gilt $\dot{\sigma}_i=x_{i+1}-x_i$. Für die globale Steifigkeitsmatrix $A$ folgt somit

$$A_{j,k} =  \sum_{i=1}^n \frac{1}{h_i} \int_{0}^{1}\dot{\tilde{\varphi}}_j(t)\cdot \dot{\tilde{\varphi}}_k(t) dt\quad\text{mit}\ h_i = x_{i+1}-x_i,$$

wobei das Integral stets über das Einheitsintervall geht. Wir können daher den finite Element Ansatz auf dem Einheitsintervall studieren und die globale Matrix durch Addition an den entsprechenden Matrixeinträgen zusammenstellen (man spricht hier vom `Assembling`).

Das Produkt $\tilde{\varphi}_j(t)\cdot \tilde{\varphi}_k(t)$ ist für fast alle $j, k$ Kombinationen null. In der Abbildung {numref}`fig-fembasissupport` ist der Fall für Polynome mit Grad 1 dargestellt. In der Abbildung {numref}`fig-FEMlocalBasisProductp13` sind die Kombinationen der *lokalen* Basisfunktionen für unterschiedliche Polynomgrade dargestellt.

```{figure} fembasissupport.png
---
align: center
height: 400px
name: fig-fembasissupport
---

Produkt $\varphi_j \cdot \varphi_k$ auf dem Intervall $[0,1]$ für stückweise Polynome mit Grad 1.
```

```{figure} FEMlocalBasisProductp13.png
---
align: center
height: 400px
name: fig-FEMlocalBasisProductp13
---

Lokale Kombinationen $\varphi_j \cdot \varphi_k$ der nodalen FEM Basisfunktionen auf dem Einheitsintervall $[0,1]$.
```

Praktisch können wir die globale Matrix mit Hilfe der lokalen Elementmatrix (farbige identische Matrizen in der Abbildung {numref}`fig-fembasissupport`) berechnen.

```{prf:definition} Elementmatrizen
:label: my-def-Elementmat

Man nennt die Matrix

$$A_e = \int_0^1 \tilde{\varphi}_j'(x)\cdot \tilde{\varphi}_k'(x) dx \quad \text{für}\ j,k = 0, \ldots, p$$

**Elementsteifigkeitsmatrix**

und 

$$M_e = \int_0^1 \tilde{\varphi}_j(x)\cdot \tilde{\varphi}_k(x) dx  \quad \text{für}\ j,k = 0, \ldots, p$$

**Elementmassenmatrix**, wobei $p$ der maximale Polynomgrad der Basisfunktionen sei. 
```

Die Berechnung der globalen Matrix $A$ erfolgt nun mit Hilfe der im Voraus berechneten lokalen Elementmatrizen $A_e$. Dabei muss spezifiziert werden, wo in der globalen Matrix $A$ die Einträge der lokalen Matrix $A_e$ abhängig vom Element $T_i$ gespeichert werden müssen. Sei $T$ die Abbildung der lokalen zu den globalen Freiheitsgrade

$$T : \{0,1, \ldots,n\} \times \{0,1,\ldots, p\} \to \{0, 1, \ldots, N\}$$ (eq:lokalglobalMapping)

wobei $n$ die Anzahl Elemente, $p$ den Polynom Grad der lokalen Basisfunktionen und $N$ die globale Anzahl Freiheitsgrade bezeichne. Für den eindimensionalen Fall mit Polynome ersten Grades die Abbildung $T$ gegeben durch

$$\begin{array}{c|cccccc}
j,i & 0 & 1 & 2 & 3 & \ldots & n\\\hline
0   & 0 & 1 & 2 & 3 & \ldots & N-1\\
1   & 1 & 2 & 3 & 4 & \ldots & N.\\
\end{array}$$

Die **globalen Matrizen** können daher mit lokalen Elementmatrizen wie folgt berechnet werden

```{prf:algorithm} Assembling Bilinearform
:label: my-alg-AssemnlingBilinearform

1. for $i$ in range($n$):
    1. compute local Matrix $A_e$
    2. for $j$ in range($p$):
        1. for $k$ in range($p$):
            1. $A_{T(i,j),T(i,k)}$ += $\frac{1}{h_i}\cdot {A_e}_{j,k}$
```

```{prf:remark}
:label: my-rm-numpde1

* Im Fall, dass die zu lösende Differentialgleichung konstante Koeffizienten hat, kann die Elementmatrix vorab berechnet werden.
* Sind ortsabhängige Koeffizienten vorhanden, muss die Elementmatrix für jedes Element lokal berechnet werden. Die Definition der Basisfunktionen (shape functions) erfolgt auch in dem Fall nur auf dem Referenzelement.
```

Die Dimension der Elementmatrix ist vom Polynomgrad $p$ abhängig. In den Tabellen {numref}`chap:1dFEMElemente` sind die $S_e, M_e$ für die Polynomgrade $p=1,2,3$ berechnet, analog zur Abbildung {numref}`fig-FEMlocalBasisProductp13`.

Das Assembling der Systemmatrix (Matrix der Bilinearform) des Einstiegsbeispiels {eq}`eq:schwacheGleichungBeispiel` mit Hilfe Elemente erster Ordnung kann wie folgt umgesetzt werden:

```{code-cell} ipython3
:tags: [hide-cell, remove-output]
import numpy as np
from pandas import DataFrame
def highlight_ortho(s):
    is_ortho = np.abs(s) < 1e-13
    return ['background-color: yellow' if v else '' for v in is_ortho]
import matplotlib.pyplot as plt
from myst_nb import glue

def uanalytic(x):
    return -0.5*x*(x-1)
```
```{code-cell} ipython3
# Elementsteifigkeitsmatrix
order = 1
Ae = np.array([[1,-1],[-1,1]])

# Zerlegung des Gebiets
n=5;a=0;b=1
xi = np.linspace(a,b,n+1)
h = xi[1:]-xi[:-1]

# lokal - global mapping
N = n+1
T = np.array([[i,i+1] for i in range(N)])

# Globale Steiffigkeitsmatrix
A = np.zeros((N,N))

# Loop über Elemente (Assembling)
for i in range(n):
    for j in range(order+1):
        for k in range(order+1):
            A[T[i,j],T[i,k]] += Ae[j,k]/h[i]
A
```

Für die Linearform folgt analog

$$\int_{T_i} f(x)\, \varphi_k(x) dx = \int_{x_i}^{x_{i+1}} f(x) \, \varphi_k(x) dx  = \int_{\sigma_i^{-1}(x_i)}^{\sigma_i^{-1}(x_{i+1})} f(\sigma_i(t)) \, \varphi_k(\sigma_i(t)) \dot{\sigma}_i(t) dt$$

Mit Hilfe des lokal-global Mapping {eq}`eq:lokalglobalMapping` und den Referenz Basisfunktionen folgt

```{prf:algorithm} Assembling Linearform
:label: my-alg-AssemblingLinform

1. for $i$ in range($n$):
    1. compute local Vector $f_e$
    2. for $j$ in range($p$):
        1. $f_{T(i,j)}$ += $h_i\cdot {f_e}_{j}$
```

Für das Beispiel folgt:

```{code-cell} ipython3
from scipy.integrate import quad
# Referenz Elementfunktionen
def phi(t,i):
    if i == 0:
        return 1-t
    else:
        return t

# Rechteseite der DGL
def func(x):
    return np.ones_like(x)

# Koordinaten Transformation
def sigma(t,i):
    return xi[i]+t*h[i]

# Globale Vektor der Linearform
f = np.zeros(N)

# Loop über Elemente (Assembling)
for i in range(n):
    # berechne lokaler Vektor
    fe = [quad(lambda t: func(sigma(t,i))*phi(t,j), 0,1)[0] for j in range(order+1)]
    for j in range(order+1):
        f[T[i,j]] += h[i]*fe[j]
f
```

Mit der Berücksichtigung der Dirichlet Randbedingung folgt die numerische Lösung, abgebildet in der Abbildung {numref}`FEM_1d_p1_solutionexmp_fig2`.

```{code-cell} ipython3
:tags: [hide-cell, remove-output]
ui = np.zeros_like(xi)
ui[1:-1] = np.linalg.solve(A[1:-1,1:-1],f[1:-1])

xp = np.linspace(0,1,400)
fig, ax = plt.subplots(figsize=(6, 2))
ax.plot(xi,ui,label='FEM Lösung')
ax.plot(xp,uanalytic(xp),label='exakte Lösung')
ax.legend()
glue("FEM_1d_p1_solutionexmp_fig2", fig, display=False)
```

```{glue:figure} FEM_1d_p1_solutionexmp_fig2
---
figwidth: 400px
name: FEM_1d_p1_solutionexmp_fig2
---

FEM 1d affine Basisfunktionen
```
