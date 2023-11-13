---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Finite Elemente Methode Einstiegsbeispiele

## Eindimensionaler Fall

Zum Einstieg in die Methode der finite Elemente kommen wir auf das skalare Randwertproblem {eq}`eq:eindimrwp`, dem Poisson Problem zurück

$$\begin{split}
-u''(x) & = f(x)\quad\forall\ x\in (0,1)\\
u(0) & = u(1) = 0,
\end{split}$$ (eq:starkeGleichungBeispiel)

mit $f(x) = 1$.

Die analytische Lösung erhalten wir in dem Fall leicht. Durch zweimaliges Integrieren der rechten Seite erhalten wir ein Polynom  2. Grades

$$u(x) = \frac{1}{2} x^2 + C_1 x + C_2.$$

Durch Einsetzen der Randbedingungen folgt die analytische Lösung

$$u(x) = -\frac{1}{2} x (x-1).$$

```{code-cell} ipython3
def uanalytic(x):
    return -0.5*x*(x-1)
```

Für die numerische Lösung multiplitzieren wir die Differentialgleichung mit einer beliebigen Testfunktion $v(x)\in C_0^\infty(0,1)$ und integrieren über das Intervall $(0,1)$

$$-\int_0^1 u''(x) v(x)\,dx = \int_0^1 f(x) v(x)\,dx \quad\text{für alle}\ v\in C_0^\infty(0,1).$$

Mit Hilfe der partiellen Integration erhalten wir die schwache Gleichung

$$\int_0^1 u'(x) v'(x)\,dx - \underbrace{\big[u'(x) v(x)\big]_0^1}_{=0\, \text{da}\, v(0)=v(1)=0} = \int_0^1 f(x) v(x)\,dx \quad\text{für alle}\ v\in C_0^\infty(0,1).$$

Gesucht ist eine Funktion $u(x)$ so, dass $u(0)=u(1)=0$ und

$$\int_0^1 u'(x) v'(x)\,dx  = \int_0^1 f(x) v(x)\,dx \quad\text{für alle}\ v\in C_0^\infty(0,1).$$ (eq:schwacheGleichungBeispiel)

Um der Gleichung zu genügen, muss die Lösung $u$ nicht zwingend eine zweimal stetig differenzierbare Funktion sein. Anstelle dessen werden $u$ und $v$ jeweils einmal differenziert. Die Funktionen liegen im Sobolev-Raum $H_0^1(0,1)$. Die rechte Seite muss ebenfalls nicht mehr zwingend stetig sein. Es reicht, wenn die Funktion $f$ quadratisch integrierbar ist, daher $f\in L_2(0,1)$. Die finite Elemente Methode benutzt die eingeführte verallgemeinerte Ableitung {eq}`eq:VerallgemeinerteAbleitung`. Wir erhalten Lösungen, welche nicht immer zwingend auch Lösung der starken Gleichung {eq}`eq:starkeGleichungBeispiel` sein müssen.

Wir diskretisieren nun das Intervall $(0,1)$ in Teilintervalle $(x_{i-1},x_i)$ für $i=1,\ldots, n$. Auf dieser Zerlegung definieren wir die Basisfunktionen

$$\varphi_i(x) = \begin{cases}
\frac{x-x_{i-1}}{x_i-x_{i-1}}\quad \text{für}\ x \in (x_{i-1},x_i)\\
\frac{x_{i+1}-x}{x_{i+1}-x_{i}}\quad \text{für}\ x \in (x_{i},x_{i+1})\\
0\quad \text{sonst,}\end{cases}$$

mit welchen wir den Sobolevraum $H_0^1(0,1)$ approximieren.

### Direkte (numpy) Lösung

```{code-cell} ipython3
:tags: [hide-cell, remove-output]

import numpy as np
import matplotlib.pyplot as plt
from myst_nb import glue

n=5
a=0
b=1
xi = np.linspace(a,b,n+1)

def phi(i,x,xi=xi):
    """
    1d affine FEM Basisfunktionen
    
    Parameters
    ----------
    i : int, Nummer der Basisfunktion
    x : nparray, Koordinaten für Auswertung
    xi : nparray, Knoten

    Return
    ------
    y : nparray, Funktionswerte für x
    """

    y = np.zeros_like(x)
    if i > 0:
        ind = (xi[i-1]<=x)*(x<=xi[i])
        y[ind] = (x[ind]-xi[i-1])/(xi[i]-xi[i-1])
    if i < n:
        ind = (xi[i]<=x)*(x<xi[i+1]) 
        y[ind] = (xi[i+1]-x[ind])/(xi[i+1]-xi[i])
    return y

xp = np.linspace(a,b,400)
fig, ax = plt.subplots(figsize=(6, 2))
for i in range(n+1):
    ax.plot(xp,phi(i,xp),label=r'$\varphi_'+str(i)+'$')
ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1))
glue("FEM_1d_p1_fig", fig, display=False)
glue("FEM_1d_p1_n",n,display=False)
```

In der Abbildung {numref}`fig-FEM1dBasisFkt` ist das Intervall in {glue:}`FEM_1d_p1_n` Teilintervalle zerlegt. Entsprechend haben wir {glue:}`FEM_1d_p1_n` +1 Basisfunktionen. Abgesehen von der ersten und letzten Basisfunktion erstreckt sich der Support jeweilen über zwei Teilintervalle. Die Basisfunktionen sind genau in einem Knoten $x_i$ eins, in allen anderen beträgt der Funktionswert 0. 

```{glue:figure} FEM_1d_p1_fig
---
figwidth: 400px
name: fig-FEM1dBasisFkt
---

FEM 1d affine Basisfunktionen
```

Für die schwache Gleichung {eq}`eq:schwacheGleichungBeispiel` erhalten wir das endlich dimensionale mit finite Elemente diskrete Problem

Gesucht ist eine Funktion

$$u_h(x) = \sum_{i=0}^n u_i \, \varphi_i(x)$$

so, dass $u(0) = u(1) = 0$ und 

$$\int_0^1 u_h'(x) \varphi_j'(x)\,dx  = \int_0^1 f(x) \varphi_j(x)\,dx \quad\text{für alle}\ j = 1, \ldots, n-1.$$ (eq:diskreteschwacheGleichungBeispiel)

Setzt man die Darstellung für $u_h$ ein, so folgt

$$\sum_{i=0}^n u_i\ \int_0^1 \varphi_i'(x) \varphi_j'(x)\,dx  = \int_0^1 f(x) \varphi_j(x)\,dx \quad\text{für alle}\ j = 1, \ldots, n-1.$$ 

Wir definieren die Matrix $A$ und den Vektor $b$

$$\begin{split}
A_{j,i} & = \int_0^1 \varphi_i'(x) \varphi_j'(x)\,dx\\
b_j & = \int_0^1 f(x) \varphi_j(x)\,dx
\end{split}$$ (eq:linBilinForm1dproblem)

und erhalten so das (reduzierte) lineare Gleichungssystem für die Koeffizienten $u_1, \ldots, u_{n-1}$

$$\sum_{i=1}^{n-1} A_{j,i} u_i = b_j \quad \forall\ j=1, \ldots, n-1.$$(eq:FEM1dLinOrderSys)

Um die Matrix $A$ aufzustellen werden die Ableitungen der Basisfunktionen $\varphi_i$ benötigt. Gehen wir wie im Beispiel oben {numref}`fig-FEM1dBasisFkt` von einer konstanten Unterteilung aus, so gilt

$$\varphi_i'(x) = \begin{cases}
1/h\quad \text{für}\ x\in (x_{i-1},x_i)\\
-1/h\quad \text{für}\ x\in (x_{i},x_{i+1})\\
0 \quad \text{sonst.}\end{cases}$$

```{code-cell} ipython3
:tags: [hide-cell, remove-output]

def dphi(i,x,xi=xi):
    """
    Ableitung 1d affine FEM Basisfunktionen
    
    Parameters
    ----------
    i : int, Nummer der Basisfunktion
    x : nparray, Koordinaten für Auswertung
    xi : nparray, Knoten

    Return
    ------
    y : nparray, Funktionswerte für x
    """
    y = np.zeros_like(x)
    if i > 0:
        h = xi[i]-xi[i-1]
        ind = (xi[i-1]<=x)*(x<=xi[i])
        y[ind] = np.ones_like(x[ind])/h
    if i < n:
        h = xi[i+1]-xi[i]
        ind = (xi[i]<=x)*(x<xi[i+1]) 
        y[ind] = -np.ones_like(x[ind])/h
    return y

fig, ax = plt.subplots(figsize=(6, 2))
for i in range(n+1):
    ax.plot(xp,dphi(i,xp),label=r'$\varphi_'+str(i)+'\'(x)$')
ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.05))
glue("FEM_1d_p1_deriv_fig", fig, display=False)
```

```{glue:figure} FEM_1d_p1_deriv_fig
---
figwidth: 400px
name: fig-FEM1dDerivBasisFkt
---

Ableitung FEM 1d affine Basisfunktionen
```

Für das konkrete Beispiel erhalten wir für die Matrixkoeffizienten die Matrix $A$:

```{code-cell} ipython3
:tags: [hide-input]

from scipy.integrate import fixed_quad

A = []
for i in range(0,n+1):
    ai = []
    for j in range(0,n+1):
        # Integration über die Elemente mit Hilfe der Gauss-Quadratur
        aij = np.sum([fixed_quad(lambda x:dphi(i,np.array(x))*dphi(j,np.array(x)), xi[k], xi[k+1],n=2)[0]
                      for k in range(n)])
        ai.append(aij)
    A.append(ai)
A = np.array(A,dtype=float)
A
```

```{admonition} Aufgabe
Berechne die Koeffizienten von Hand.
```

Analog zur Matrix $A$ folgt für die rechte Seite mit konkreter Funktion $f(x) = 1$:

```{code-cell} ipython3
:tags: [hide-input]

def f(x):
    x = np.array(x)
    return np.ones_like(x)

b = []
for j in range(0,n+1):
    b.append(np.sum([fixed_quad(lambda x:f(np.array(x))*phi(j,np.array(x)),
                                 xi[k], xi[k+1],n=2)[0]
                      for k in range(n)]))
b = np.array(b,dtype=float)
b
```

Die FEM Lösung $u(x)$ ist somit gegeben durch die Lösung des Gleichungssystems {eq}`eq:FEM1dLinOrderSys` mit der berechneten Matrix und Vektor. Da die Randwerte gegeben sind, werden nur die inneren Freiheitsgrade benutzt. Es folgt

```{code-cell} ipython3
u = np.zeros_like(xi)

# um die Wahl der linearen Gleichungslöser kümmern wir uns später:
u[1:-1] = np.linalg.solve(A[1:-1,1:-1],b[1:-1])
u
```

```{code-cell} ipython3
:tags: [hide-cell, remove-output]

xp = np.linspace(0,1,400)
fig, ax = plt.subplots(figsize=(6, 2))
ax.plot(xi,u,label='FEM Lösung')
ax.plot(xp,uanalytic(xp),label='exakte Lösung')
ax.legend()
glue("FEM_1d_p1_solutionexmp_fig", fig, display=False)
```

```{glue:figure} FEM_1d_p1_solutionexmp_fig
---
figwidth: 400px
name: fig-FEM1dSolutionExmp
---

Lösung des Randwertproblem mit Hilfe 1. Ordnung FEM
```

```{admonition} Aufgaben
* Wie lautet das Gleichungssystem, wenn für die Integration der rechten Seite die Trapezregel benutzt wird?
* Berechne die Lösung mit Hilfe der finiten Differenzen Methode und vergleiche die beiden Systeme.
```

### NGSolve Lösung

Das selbe nun mit NGSolve

* Wir erstellen als erstes ein eindimensionales Mesh.

```{code-cell} ipython3
from netgen.meshing import Mesh as NGMesh # Vorsicht es gibt Mesh auch in ngsolve!
from netgen.meshing import MeshPoint, Pnt, Element1D, Element0D
from ngsolve import *
```

```{code-cell} ipython3
m = NGMesh(dim=1)

# Anzahl Teilintervalle
N = 5

# Punkte für die Zerlegung auf dem Intervall [0,1]
pnums = []
for i in range(0, N+1):
    pnums.append (m.Add (MeshPoint (Pnt(i/N, 0, 0))))

# Jedes 1D-Element (Teilintervall) kann einem Material zugeordnet
# werden. In unserem Fall gibt es nur ein Material.
idx = m.AddRegion("material", dim=1)
for i in range(0,N):
    m.Add (Element1D ([pnums[i],pnums[i+1]], index=idx))

# Linkes und Rechtes Ende sind Randwertpunkte (0D-Elemente)
idx_left = m.AddRegion("left", dim=0)
idx_right = m.AddRegion("right", dim=0)

m.Add (Element0D (pnums[0], index=idx_left))
m.Add (Element0D (pnums[N], index=idx_right))

# Damit haben wir das Mesh definiert
mesh = Mesh(m)
```

Nun erstellen wir einen $H^1$ Funktionenraum mit Hilfe dieses 1D Mesh und den Dirichlet Randpunkte `left` und `right`.

```{code-cell} ipython3
V = H1(mesh,order = 1, dirichlet='left|right')
u,v = V.TnT()
```

$u,v$ sind Trial und Test Funktionen für die Definition der Linear- und Bilinearfunktion

```{code-cell} ipython3
a = BilinearForm(V)
a += grad(u)*grad(v)*dx

f = CoefficientFunction(1)
b = LinearForm(V)
b += f*v*dx
```

Damit sind die beiden Operatoren definiert, jedoch noch nicht berechnet. Das Berechnen nennt man auch **assembling** Zusammenstellen. Wir werden später sehen, was damit gemeint ist.

```{code-cell} ipython3
a.Assemble()
b.Assemble();
```

Die Matrix und der Vektor der Bilinear- und Linearform beinhaltet sämtliche Freiheitsgrade, insbesondere also auch die Randpunkte. Diese sind jedoch durch die Dirichletrandwerte gegeben und müssen nicht berechnet werden. Dem werden wir beim Lösen des Systems rechnungtragen.

```{code-cell} ipython3
print(a.mat)
```

```{code-cell} ipython3
print(b.vec)
```

Die Lösung selber wird in einer `GridFunction` gespeichert. Die Trial und Test Functions haben zwar die gleiche Struktur, jedoch keinen Memory. Hier sind nur die Freiheitsgrade etc. des FE-Raumes gespeichert. Für die Lösung benötigen wir eine GridFunction. In der Auswertung dieser wird Linearkombination der Basisfunktionen automatisch berechnet.

```{code-cell} ipython3
gfu = GridFunction(V)
```

Berechnung der FEM Lösung mit NGSolve:

```{code-cell} ipython3
gfu.vec.data = a.mat.Inverse(freedofs=V.FreeDofs())*b.vec
```

Bei der Ausführung des Befehls wird nicht die Matrix Invertiert und von links an den Vektor multipliziert. Das wäre numerisch viel zu aufwändig. Auch wenn die Notation anderes behauptet, es wird nur das Gleichungssystem gelöst.

Es folgt das selbe Resultat wie oben:

```{code-cell} ipython3
:tags: [hide-input]

xp = np.linspace(0,1,400)
plt.plot(xp,[gfu(mesh(xi,0)) for xi in xp],label='numerische Lösung')
plt.plot(xp, uanalytic(xp),label='exakte Lösung')
plt.legend()
plt.grid()
plt.show()
```

## Zweidimensionaler Fall

Wir betrachten nun das Beispiel Randwertproblem

$$\begin{split}
-\Delta u & = 1 \quad\text{für}\ x\in \Omega\\
u & = 0 \quad\text{für}\ x\in \partial\Omega
\end{split}$$

auf dem einheits Rechteck im $\mathbb{R}^2$. Analog zum Einstiegsbeispiel {numref}`ref:IntroPoissonSchwacheGleichung` folgt die schwache Gleichung

$$\int_\Omega \nabla u \cdot \nabla \varphi_j dx = \int_\Omega f(x) \varphi_j dx \quad \text{für alle}\ j\in J$$

wobei mit $J$ die Basisfunktionen bezeichnet sind, deren Maximum im Innern des Einheitsquadrats angenommen werden.

```{code-cell} ipython3
from ngsolve import *
from ngsolve.webgui import Draw
```

Von `netgen` importieren wir das Einheits Quadrat. Mit Hilfe des netgen Moduls können wir 2d und 3d Geometrien definieren und vernetzen.

```{code-cell} ipython3
from netgen.geom2d import unit_square
```

Wir generieren nun unstrukturiertes Mesh mit der maximalen Kantenlänge von 0.25:

```{code-cell} ipython3
mesh = Mesh(unit_square.GenerateMesh(maxh=0.25))
```

Damit erhalten wir das Mesh, wobei die Feinheit über den `maxh` Parameter gesteuert werden kann.

```{code-cell} ipython3
Draw(mesh);
```

Wie sehen im zweidimensionalen die Basisfunktionen aus? Dazu erstellen wir einen FEM Funktionenraum und visualisieren die Basisfunktionen. Der Rand des Einheitsquadrats besteht aus 4 Linien, welche verschiedene Labels haben:

```{code-cell} ipython3
mesh.GetBoundaries()
```

Wir definieren den `H1` FEM Funktionenraum mit der Dirichletrandbedingung und initialisieren eine `GridFunction`, eine FEM Funktion aus dem Funktionenraum:

```{code-cell} ipython3
V = H1(mesh,order=1,dirichlet='bottom|right|top|left')
gfu = GridFunction(V)
```

In der Gridfunction wird der Lösungsvektor, die Koeffizienten der Linearkombination der Basisfunktionen gespeichert. In diesen Vektor schreiben wir nun an unterschiedlichen Stellen eine 1, womit wir die verschiedenen Basisfunktionen visualisieren können.

```{code-cell} ipython3
gfu.vec.FV()[:] = 0 # alle Einträge mit 0 initialisieren
gfu.vec.FV()[20] = 1
```

```{code-cell} ipython3
Draw(gfu,mesh,'u');
```

Die inneren Freiheitsgrade sind in unserem Fall die freien Freiheitsgrade des FEM Raumes, daher in der Zahl 9 Stück:

```{code-cell} ipython3
freedofs = V.FreeDofs()
print (freedofs)
freedofsnp = np.array([i for i in freedofs])
```

Wir setzen nun alle inneren Freiheitsgrade auf 1 und erhalten damit die Linearkombination der Basisfunktionen:

```{code-cell} ipython3
gfu.vec.FV().NumPy()[freedofsnp] = 1.
Draw(gfu);
```

Die Ableitungen der Basisfunktionen sind stückweise konstante Funktionen. Überschlagen wir hier kurz die Steigung: $1/4$-tel ist die Seitenlänge der Zerlegung. Das bedeutet, dass die Steigung in $x/y$ Richtung 4 beträgt.

Visualisieren wir von $\nabla u$ die $x$-Komponente, so sollte, das einen Graph mit stückweise $\pm 4$ in den Elementen ergeben, welche in $x$-Richtung steigen bzw. fallen.

```{code-cell} ipython3
Draw(grad(gfu)[0],mesh);
```

und in $y$ Richtung

```{code-cell} ipython3
Draw(grad(gfu)[1],mesh);
```

Mit Hilfe von `ngsolve` können wir nun die Bilinearfunktion $A: V \times V \to \mathbb{R}$ und die Linearfunktion $f: V \to \mathbb{R}$ sehr einfach berechnen.

Es folgt für 

$$\begin{split}
A : V \times V & \to \mathbb{R}\\
    (u,v) & \mapsto A(u,v) = \int_\Omega \nabla u\cdot \nabla v dx,
\end{split}$$

wobei hier $u,v$ stellvertretend für $\varphi_i, \varphi_j\in V$ steht

```{code-cell} ipython3
u,v = V.TnT()
```

```{code-cell} ipython3
A = BilinearForm(V)
A += grad(u) * grad(v)*dx
```

und analog für die Linearform

$$\begin{split}
b : V & \to \mathbb{R}\\
    v & \mapsto b(v) = \int_\Omega f(x) v dx,
\end{split}$$

```{code-cell} ipython3
f = CoefficientFunction(1)
b = LinearForm(V)
b += f*v*dx
```

Die gegebene Funktion $f(x)$ der rechten Seite können wir mit sogenannten `CoefficientFunction` definieren.

Damit haben wir die Bilinearform, welche im endlichdimensionalen mit Hilfe einer Matrix und die Linearform, welche mit einem Vektor beschrieben werden kann definiert, aber noch nicht berechnet. Die Berechnung derer nennt man auch **Assembling**.

```{code-cell} ipython3
A.Assemble()
b.Assemble();
```

Die Matrix $A$ ist "sparse" gespeichert. Das bedeutet, dass nur die Matrix Einträge gespeichert werden, für welche wir potentiell einen Eintrag erhalten. Für den ganzen Rest der Matrix wird der Speicher gar nicht allokiert. Grundsätzlich haben wir hier alle Freiheitsgrade:

```{code-cell} ipython3
print(A.mat)
```

Wir können diese Matrix (zumindest solange sie klein ist!) auch als vollbesetzte "dense" Matrix betrachten:

```{code-cell} ipython3
rows,cols,vals = A.mat.COO()

denseA = np.zeros((np.max(rows)+1,np.max(rows)+1))
k=0
for i,j in zip(rows,cols):
    denseA[i,j] = vals[k]
    k+=1

plt.spy(denseA)
plt.show()
```

Das Pattern der Matrixeinträge hängt von der Nummerierung der Knoten (aus dem Meshing) und damit der Freiheitsgrade für den FEM Ansatz erster Ordnung ab.

```{code-cell} ipython3
for e in mesh.edges:
    line = np.array([mesh.vertices[v.nr].point for v in e.vertices])
    plt.plot(line[:,0],line[:,1],c='gray',alpha=0.75)
for v in mesh.vertices:
    plt.text(*v.point,v,color='red')
plt.gca().set_axis_off()
plt.gca().set_aspect(1)
plt.show()
```

```{code-cell} ipython3
ind = np.arange(freedofsnp.shape[0])[freedofsnp]
plt.spy(denseA[np.ix_(ind,ind)])
plt.show()
```

Lösen wir das System für die inneren Freiheitsgrade:

```{code-cell} ipython3
gfu.vec.data = A.mat.Inverse(freedofs=V.FreeDofs())*b.vec
```

```{code-cell} ipython3
Draw(gfu,mesh,'u');
```
