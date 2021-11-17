---
jupytext:
  formats: ipynb,md:myst
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

(chap:1dFEMElemente)=
# Element Matrizen 1d Fall

```{code-cell} ipython3
:tags: [hide-cell]
import numpy as np
import matplotlib.pyplot as plt

from sympy import integrate
from sympy.abc import x

from pandas import DataFrame
from IPython.display import display
```

Wir benutzen Lagrange Polynome

$$\varphi_j(x) = \prod_{i\not=j}^n \frac{x-x_i}{x_j-x_i}$$

verschiedener Ordnung als FEM Basis Funktionen:

```{code-cell} ipython3
def lagrangePoly(x,j,order):
    xi = np.linspace(0,1,order+1)
    J = np.delete(np.arange(order+1),j)
    return np.product([(x-xi[i])/(xi[j]-xi[i]) for i in J],axis=0)
```

Die Polynome sind gegeben durch:

```{code-cell} ipython3
:tags: [hide-input]
for order in [1,2,3]:
    print('Order = ',order)
    for i in range(order+1):
        display(lagrangePoly(x,i,order).expand())
```

```{code-cell} ipython3
:tags: [hide-input]
xp = np.linspace(0,1,400)
col = ['tab:blue','tab:orange','tab:green','tab:red']
for order in [1,2,3]:
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
for order in [1,2,3]:
    print('order = ',order)
    m = [[integrate(lagrangePoly(x,i,order).diff()*lagrangePoly(x,j,order).diff(),(x,0,1))
                     for i in range(order+1)] for j in range(order+1)]
    df = DataFrame(m)
    display(df.style.format('{:.4f}'))
    print()
```

Für die Massen-Elementmatrizen

$$A_{i,j} = \int_0^1 \varphi_i(x) \varphi_j(x) dx$$

erhalten wir

```{code-cell} ipython3
:tags: [hide-input]
for order in [1,2,3]:
    print('order = ',order)
    m = [[integrate(lagrangePoly(x,i,order)*lagrangePoly(x,j,order),(x,0,1))
                     for i in range(order+1)] for j in range(order+1)]
    df = DataFrame(m)
    display(df.style.format('{:.4f}'))
    print()
```
