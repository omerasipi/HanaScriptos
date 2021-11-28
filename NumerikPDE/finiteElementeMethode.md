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

# Methode der finiten Elemente

## Assembling

### Eindimensionaler Fall

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
align: left
height: 400px
name: fig-fembasissupport
---

Produkt $\varphi_j \cdot \varphi_k$ auf dem Intervall $[0,1]$ für stückweise Polynome mit Grad 1.
```

```{figure} FEMlocalBasisProductp13.png
---
align: left
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

### Zweidimensionaler Fall

Eine reguläre Triangulierung $\mathcal{T} = \{T_1, \ldots, T_M\}$ eines Gebiets $\Omega$ ist die Zerlegung in Dreiecke $T_i$ so, dass $\bar{\Omega} = \cup_i T_i$ und $T_i\cap T_j$ ist
* entweder leer
* oder hat eine gemeinsame Kante.

In einem erweiterten Sinne kann die Triangulierung aus verschiedenen Elemente bestehen: Dreiecke, Vierecke, (Tetraeder, Hexeder, Prismen, Pyramiden im dreidimensionalen). Die finite Elemente werden typischerweise, wie wir es im eindimensionalen gemacht haben, auf einem Referenzelement definiert. Die einzelnen Elemente der Zerlegung können mit Hilfe einer affinen Transformation und dem Referenzelement beschrieben werden. 

Die Koordinaten Transformation im zweidimensionalen erfordert etwas mehr Rechenaufwand.

```{figure} KoordinatenTransformation2d.png
---
align: left
height: 250px
name: fig-KoordinatenTransformation2d
---

Transformation auf Einheitsdreieck in $2D$
```

Ein Dreieck $T_i$ in allgemeiner Lage mit den Eckpunkten $P_1(x_1,y_1)$, $P_2(x_2,y_2)$ und  $P_3(x_3,y_3)$, welche im Gegenuhrzeigersinn fortlaufend numeriert seien, wie dies in Abb. {numref}`fig-KoordinatenTransformation2d` erfolgte, kann mittels der linearen Transformation

$$\begin{split}
x & = x_1 + (x_2-x_1)\, \xi + (x_3-x_1)\, \eta\\
y & = y_1 + (y_2-y_1)\, \xi + (y_3-y_1)\, \eta
\end{split}$$

bijektiv auf das Einheitsdreieck $T$ abgebildet werden. Für das Dreieck $T_i$ gegeben durch die drei Punkte $\vec{p}_{i,1}, \vec{p}_{i,2}, \vec{p}_{i,3}$ definieren wir die Matrix

$$A_i = (\vec{p}_{i,2}-\vec{p}_{i,1},\ \vec{p}_{i,3}-\vec{p}_{i,1}) \in \mathbb{R}^{2\times 2}.$$

Damit können wir die Transformation wie folgt schreiben

$$\begin{array}{rl}\sigma_i : T & \to  T_i\\
\vec{t} & \mapsto \begin{pmatrix}x\\ y\end{pmatrix} = \sigma_i(\xi,\eta) = \vec{p}_{i,1} + A_i\cdot \begin{pmatrix}\xi\\\eta\end{pmatrix}\end{array}$$

bzw. die inverse Transformation

$$\begin{array}{rl}\sigma_i^{-1} : T_i & \to T\\
\vec{x} & \mapsto \begin{pmatrix}\xi\\ \eta\end{pmatrix} = \sigma_i^{-1}(x,y) = A_i^{-1}\cdot\left(\begin{pmatrix}x\\ y\end{pmatrix} - \vec{p}_{i,1}\right)\end{array} 
$$


Das Flächenelement $dx dy$ kann mit der Jacobi-Determinante $J = \det A_i$ der Transformation durch $dx\, dy = J\, d\xi\, d\eta$ ersetzt werden. Die Basisfunktion $\varphi$ auf dem Element $T_i$ kann (analog zum eindimensionalen Fall) mit Hilfe der Basisfunktion $\tilde{\varphi}$ auf dem Referenzelement beschrieben werden. Es gilt

$$\varphi(\vec{x}) = \tilde{\varphi}(\sigma_i^{-1}(\vec{x})).$$ (eq:femlokalbasisfkt2d)

analog zu {eq}`eq:femlokalbasisfkt`. Für die Jacobimatrix folgt mit Hilfe der Kettenregel

$$D\varphi(\vec{x}) = D\tilde{\varphi}(\sigma_i^{-1}(\vec{x}))\cdot \underbrace{D\sigma_i^{-1}(\vec{x})}_{=A_i^{-1}}$$

und damit für den Gradient

$$\nabla \varphi(\vec{x}) = D\varphi(\vec{x})^T = (D\tilde{\varphi}(\sigma_i^{-1}(\vec{x}))\cdot A_i^{-1})^T = {A_i^{-1}}^T\cdot \nabla \tilde{\varphi}(\sigma_i^{-1}(\vec{x})).$$

Für die Elementsteifigkeitsmatrix folgt somit

$$\int_{T_i} \nabla \varphi_j(\vec{x})\cdot \nabla \varphi_k(\vec{x}) dx dy = \int_{T} \nabla \tilde{\varphi}_j(\vec{t})\cdot C_i\cdot \nabla \tilde{\varphi}_k(\vec{t})\, J_i\, d\xi d\eta \quad \text{für}\ j,k = 0,1,\ldots, N_e$$ (eq:Elementsteifigkeitmatrix2d)

mit 

$$C_i = {A_i^{-1}}\cdot {A_i^{-1}}^T$$

und für die Elementmassenmatrix

$$\int_{T_i} \varphi_j(\vec{x})\cdot \varphi_k(\vec{x}) dx dy = \int_{T} \tilde{\varphi}_j(\vec{t}) \tilde{\varphi}_k(\vec{t})\, J_i\, d\xi d\eta \quad \text{für}\ j,k = 0,1,\ldots, N_e.$$ (eq:Elementmassenmatrix2d)

**Beispiel**: Betrachten wir Basisfunktionen erster Ordnung 

$$\begin{split}
\tilde{\varphi}_0(\xi,\eta) & = 1-\xi-\eta\\
\tilde{\varphi}_1(\xi,\eta) & = \xi\\
\tilde{\varphi}_2(\xi,\eta) & = \eta\end{split}$$

auf dem Einheitsdreieck, gegeben durch die Punkte $(0,0), (1,0), (0,1)$. In dem Fall gilt $N_e = 2$.

```{code-cell} ipython3
# Basisfunktionen auf dem Referenzelement
def myshape(t, j):
    xi, eta = t
    if j == 0:
        return 1-xi-eta
    elif j == 1:
        return xi
    else:
        return eta

# Gradienten der Basisfunktionen auf dem Referenzelement
def Dmyshape(t, j):
    if j == 0:
        return np.array([-1,-1])
    elif j == 1:
        return np.array([1,0])
    else:
        return np.array([0,1])
```

Die inverse Jacobi-Matrix der Transformation ist gegeben durch

```{code-cell} ipython3
def invDsigma(p1,p2,p3):
    return np.linalg.inv(np.array([(p2-p1),(p3-p1)]).T)
```

Damit können wir für ein beliebiges Dreieck die Elementsteifigkeitsmatrix {eq}`eq:Elementsteifigkeitmatrix2d` wie folgt berechnen

```{code-cell} ipython3
# Berechnung der C_i Matrix und Jacobi Determinante
def CnJ(p1,p2,p3):
    A = np.array([(p2-p1),(p3-p1)]).T
    invA = np.linalg.inv(A)
    return invA@invA.T, np.linalg.det(A)

# Gebietsintegration über Einheitsdreieck
from scipy.integrate import dblquad
def quadT(f):
    return dblquad(f, 0, 1, 0, lambda x: 1-x)[0]

# Integration über Dreieck gegeben durch (1,1), (1.5,-1), (2,1.2)
Ci, Ji = CnJ(np.array([1,1]),np.array([1.5,-1]),np.array([2,1.2]))

Ai = np.array([[ quadT(lambda x,y: Dmyshape([x,y], j)@Ci@Dmyshape([x,y], k)*Ji)
           for k in range(3)] for j in range(3)])
Ai = DataFrame(Ai)
Ai.style.\
    apply(highlight_ortho).\
    set_table_attributes('style="font-size: 12px"')
```

Die Matrix $C_i$ kann analytisch berechnet werden so, womit die Invertierung hinfällig wird. Für $C_i\cdot J_i^2$ erhalten wir

```{code-cell} ipython3
:tags: [hide-input]
from sympy import symbols, Matrix
p1x,p1y, p2x,p2y, p3x,p3y, = symbols('p_1_x,p_1_y,p_2_x,p_2_y,p_3_x,p_3_y')
p1 = np.array([p1x,p1y])
p2 = np.array([p2x,p2y])
p3 = np.array([p3x,p3y])
A=Matrix(np.array([p2-p1,p3-p1]).T).as_immutable()
J = A.det()
invA=A.inverse().simplify()
CJ2=invA@invA.T*J**2
CJ2.simplify()
```

und damit

$$C_i = \frac{1}{J_i^2} \begin{pmatrix}a_i & -c_i\\ -c_i & b_i\end{pmatrix}$$ (eq:Transmatrix2d1)

mit

$$\begin{split}
a_i & = \|\vec{p}_{i,3}-\vec{p}_{i,1}\|^2 = (p_{i,3,x}-p_{i,1,x})^2+(p_{i,3,y}-p_{i,1,y})^2\\
b_i & = \|\vec{p}_{i,2}-\vec{p}_{i,1}\|^2 = (p_{i,2,x}-p_{i,1,x})^2+(p_{i,2,y}-p_{i,1,y})^2\\
c_i & = (p_{i,2,x}-p_{i,1,x})(p_{i,3,x}-p_{i,1,x}) + (p_{i,2,y}-p_{i,1,y})(p_{i,3,y}-p_{i,1,y})\\
J_i & = (p_{i,2,x}-p_{i,1,x})(p_{i,3,y}-p_{i,1,y}) - (p_{i,3,x}-p_{i,1,x})(p_{i,2,y}-p_{i,1,y}).
\end{split}$$ (eq:Transmatrix2d2)

Für die Elementmassenmatrix {eq}`eq:Elementmassenmatrix2d` folgt

```{code-cell} ipython3
Bi = np.array([[ quadT(lambda x,y: myshape([x,y], j)*myshape([x,y], k)*Ji)
           for k in range(3)] for j in range(3)])
Bi = DataFrame(Bi)
Bi.style.\
    apply(highlight_ortho).\
    set_table_attributes('style="font-size: 12px"')
```

Das Assembling des Gesamtsystems erfolgt mit dem Algorithmus {prf:ref}`my-alg-AssemnlingBilinearform` und der Abbildung der Freiheitsgrade $T$ analog zu {eq}`eq:lokalglobalMapping` für den zweidimensionalen Fall

$$T: \{0,1,\ldots, n\} \times \{0,1,2\} \to \{0,1,\ldots, N\},$$

wobei $N$ im Fall von Elemente erster Ordnung die Anzahl Knoten und $n$ die Anzahl Elemente (Dreiecke).

```{prf:remark}
:label: my-rm-numpde2

Es erweist sich mathematisch wie auch Software technisch als bedeutend effizienter, die Berechnung der System Matrizen über die Triangulierung zu berechnen und die einzelnen Beiträge in der globalen Matrix aufzukummulieren. Diesen Prozess nennt man **Assembling**.

Allgemein können wir dies in der Form

$$A = \sum_{T \in \mathcal{T}} C_T A_T C_T^T$$

und 

$$f = \sum_{T \in \mathcal{T}} C_T f_T$$

schreiben, wobei $C_T$ die Verknüpfung zwischen lokalen und globalen Funktionen darstellt.
```
