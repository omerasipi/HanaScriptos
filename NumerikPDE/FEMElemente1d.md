---
jupytext:
  formats: ipynb,md:myst
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

(chap:1dFEMElemente)=
# Element Matrizen 1d Fall

## Lagrange Polynome als Basisfunktionen

```{code-cell} ipython3
:tags: [hide-cell]

import numpy as np
import matplotlib.pyplot as plt

from sympy import integrate
from sympy.abc import x

from pandas import DataFrame
from IPython.display import display

def highlight_ortho(s):
    is_ortho = np.abs(s) < 1e-13
    return ['background-color: yellow' if v else '' for v in is_ortho]
```

Wir benutzen Lagrange Polynome

$$\varphi_j(x) = \prod_{i\not=j}^n \frac{x-x_i}{x_j-x_i}$$

verschiedener Ordnung als FEM Basis Funktionen:

```{code-cell} ipython3
def lagrangePoly(x,j,order):
    xi = np.linspace(0,1,order+1)
    J = np.delete(np.arange(order+1),j)
    return np.prod([(x-xi[i])/(xi[j]-xi[i]) for i in J],axis=0)
```

Die Polynome sind gegeben durch:

```{code-cell} ipython3
:tags: [hide-input]

maxorder = 4
for order in range(1,maxorder+1):
    print('Order = ',order)
    for i in range(order+1):
        display(lagrangePoly(x,i,order).expand())
```

```{code-cell} ipython3
:tags: [hide-input]

xp = np.linspace(0,1,400)
col = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink']
for order in range(1,maxorder+1):
    xi = np.linspace(0,1,order+1)
    for j in range(order+1):
        plt.plot(xp,lagrangePoly(xp,j,order),c=col[j],label=r'$\varphi_'+str(j)+'(x)$')
        plt.plot(xi,lagrangePoly(xi,j,order),'o',c=col[j])
    plt.legend(loc=5)
    plt.grid()
    plt.title('Lagrange Polynome order = '+str(order))
    plt.show()
```

Für die Steifigkeit-Elementmatrizen

$$A_{i,j} = \int_0^1 \varphi_i'(x) \varphi_j'(x) dx$$

erhalten wir

```{code-cell} ipython3
:tags: [hide-input]

for order in range(1,maxorder+1):
    print('order = ',order)
    m = [[integrate(lagrangePoly(x,i,order).diff()*lagrangePoly(x,j,order).diff(),(x,0,1))
                     for i in range(order+1)] for j in range(order+1)]
    df = DataFrame(m)
    display(df.style.\
        apply(highlight_ortho).\
        set_table_attributes('style="font-size: 12px"').\
        format('{:.3f}'))
```

Für die Massen-Elementmatrizen

$$M_{i,j} = \int_0^1 \varphi_i(x) \varphi_j(x) dx$$

erhalten wir

```{code-cell} ipython3
:tags: [hide-input]

for order in range(1,maxorder+1):
    print('order = ',order)
    m = [[integrate(lagrangePoly(x,i,order)*lagrangePoly(x,j,order),(x,0,1))
                     for i in range(order+1)] for j in range(order+1)]
    df = DataFrame(m)
    display(df.style.\
        apply(highlight_ortho).\
        set_table_attributes('style="font-size: 12px"').\
        format('{:.3f}'))
```

## Hierarchische Basis Polynome

NGSolve benutzt für die finite Elemente Räume höherer Ordnung immer die Basisfunktionen aus den Räumen niederer Ordnung und erweitert diese entsprechend.

```{code-cell} ipython3
:tags: [hide-input]

from netgen.meshing import Mesh as NGMesh # Vorsicht es gibt Mesh auch in ngsolve!
from netgen.meshing import MeshPoint, Pnt, Element1D, Element0D
from ngsolve import *

m = NGMesh(dim=1)

# Punkte für die Zerlegung auf dem Intervall [0,1]
pnums = []
pnums.append (m.Add (MeshPoint (Pnt(0, 0, 0))))
pnums.append (m.Add (MeshPoint (Pnt(1, 0, 0))))

# Jedes 1D-Element (Teilintervall) kann einem Material zugeordnet
# werden. In unserem Fall gibt es nur ein Material.
idx = m.AddRegion("material", dim=1)
m.Add (Element1D ([pnums[0],pnums[1]], index=idx))

# Linkes und Rechtes Ende sind Randwertpunkte (0D-Elemente)
idx_left = m.AddRegion("left", dim=0)
idx_right = m.AddRegion("right", dim=0)

m.Add (Element0D (pnums[0], index=idx_left))
m.Add (Element0D (pnums[1], index=idx_right))

# Damit haben wir das Mesh definiert
mesh = Mesh(m)

xp = np.linspace(0,1,300)
for order in range(1,maxorder+1):
    V = H1(mesh,order = order, dirichlet='left|right')
    gfu = GridFunction(V)

    for k in range(V.ndof):
        gfu.vec[:] = 0
        gfu.vec[k] = 1
        plt.plot(xp,gfu(mesh(xp)),label=r'$\varphi_'+str(k)+'$')
    plt.legend()
    plt.title('order = '+str(order))
    plt.grid()
    plt.show()
```

Für die Steifigkeit-Elementmatrizen

$$A_{i,j} = \int_0^1 \varphi_i'(x) \varphi_j'(x) dx$$

erhalten wir

```{code-cell} ipython3
:tags: [hide-input]

for order in range(1,maxorder+1):
    print('order = ',order)
    V = H1(mesh,order = order, dirichlet='left|right')
    gfu = GridFunction(V)

    phii = GridFunction(V)
    phij = GridFunction(V)

    a = []
    for i in range(V.ndof):
        phii.vec[:]=0
        phii.vec[i]=1
        ai = []
        for j in range(V.ndof):
            phij.vec[:]=0
            phij.vec[j]=1
            ai.append(Integrate(grad(phii)*grad(phij),mesh,order=2*order))
        a.append(ai)

    df = DataFrame(a)
    display(df.style.\
        apply(highlight_ortho).\
        set_table_attributes('style="font-size: 12px"').\
        format('{:.4f}'))
```

Für die Massen-Elementmatrizen

$$M_{i,j} = \int_0^1 \varphi_i(x) \varphi_j(x) dx$$

erhalten wir

```{code-cell} ipython3
:tags: [hide-input]

for order in range(1,maxorder+1):
    print('order = ',order)
    V = H1(mesh,order = order, dirichlet='left|right')
    gfu = GridFunction(V)

    phii = GridFunction(V)
    phij = GridFunction(V)

    b = []
    for i in range(V.ndof):
        phii.vec[:]=0
        phii.vec[i]=1
        bi = []
        for j in range(V.ndof):
            phij.vec[:]=0
            phij.vec[j]=1
            bi.append(Integrate(phii*phij,mesh,order=2*order))
        b.append(bi)

    df = DataFrame(b)
    display(df.style.\
        apply(highlight_ortho).\
        set_table_attributes('style="font-size: 12px"').\
        format('{:.3f}'))
```
