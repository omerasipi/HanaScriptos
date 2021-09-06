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

# Maxwell Gleichungen

```{code-cell} ipython3
from netgen.csg import *
from ngsolve import *
from ngsolve.webgui import Draw
```

Die Quelle des Beispiels ist in der NGSolve Dokumentation {cite}`schoeberlNGSolveDoc`.

+++

## Geometrie

+++

Mit Hilfe der `netgen` Bibliothek können auch einfach 3D Geometrien beschrieben werden. Für komplexere kann auf die OCC Bibliothek zurückgegriffen werden.

```{code-cell} ipython3
def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(2,1,2)).bc("outer")

    core = OrthoBrick(Pnt(0,-0.05,0),Pnt(0.8,0.05,1))- \
           OrthoBrick(Pnt(0.1,-1,0.1),Pnt(0.7,1,0.9))- \
           OrthoBrick(Pnt(0.5,-1,0.4),Pnt(1,1,0.6)).maxh(0.2).mat("core")

    coil = (Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.3) - \
            Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.15)) * \
            OrthoBrick (Pnt(-1,-1,0.3),Pnt(1,1,0.7)).maxh(0.2).mat("coil")

    geometry.Add ((box-core-coil).mat("air"), transparent=True)
    geometry.Add (core, col=(0.5,0.5,0))
    geometry.Add (coil, col=(0,1,0))
    return geometry

geo = MakeGeometry()
# Draw (geo)
```

```{code-cell} ipython3
mesh = Mesh(geo.GenerateMesh(maxh=0.5))
mesh.Curve(5)
Draw (mesh, clipping = { "pnt" : (0,0,0), "vec" : (0,1,0) })
```

```{figure} image.png
---
align: left
height: 250px
name: MagfeldGeometrie
---
Geometrie
```

+++

## Magnetostatisches Problem

+++

Mit Hilfe Maxwell Gleichungen folgt das partielle Differentialgleichungssystem für statische Magnetfelder, gegeben durch

$$\begin{split}
\mathop{curl} H & = j\\
\mathop{div} B & = 0.
\end{split}$$

Mit Hilfe des Vektorpotential Ansatz $B = \mathop{curl} A$, motiviert durch die Forderung $\mathop{div} B = 0$ folgt die PDE

$$\mathop{curl} (\mu^{-1}(x) \mathop{curl} A) = j,$$

wobei $A$ das gesuchte Vektorpotential und $j$ eine externe Stromdichte sei.

Wir können wiederum eine schwache Gleichung berechnen, wobei der geeignete Funktionenraum in dem Fall durch den $H(\mathop{curl})$ gegeben ist. Die schwache Gleichung lautet

$$\int_\Omega \mu^{-1}(x) \mathop{curl} A \cdot \mathop{curl} \Psi dx = \int_\Omega j(x)\cdot \Psi dx\quad \forall\  \Psi\in H(\mathop{curl}, \Omega).$$ (eq:magnetostatic)


```{code-cell} ipython3
V = HCurl(mesh,order=3)
u,v = V.TnT()
gfu = GridFunction(V)
```

Das Rechengebiet $\Omega$ setzt sich hier aus den Teilgebiete
* Luft (air)
* Kern (core)
* Spule (coil)
zusammen.

```{code-cell} ipython3
mesh.GetMaterials()
```

Auf den jeweiligen Teilgebiete haben wir unterschiedliche Materialien, welche sich in der relativen Permeabilität unterscheiden. Wir definieren die relative Permeabilität wie folgt:

```{code-cell} ipython3
mur = { "core" : 1000, "coil" : 1, "air" : 1 }
mu0 = 1.257e-6
nu_coef = [ 1/(mu0*mur[mat]) for mat in mesh.GetMaterials() ]
nu = CoefficientFunction(nu_coef)
```

Um die eingeprägte Stromdichte in der Spule beschreiben zu können, ist der Richtungsvektor $w$ der Stromdichte erforderlich. Wir haben im Beispiel eine zylindrische Spule. Entsprechend definieren wir den Richtungsvektor $w$

```{code-cell} ipython3
w = CoefficientFunction((y,0.05-x,0))/sqrt(x*x+y*y)
```

Damit können wir nun die Bilinearform und Linearform des Systems {eq}`eq:magnetostatic` definieren und berechnen:

```{code-cell} ipython3
a = BilinearForm(V)
a += nu*curl(u)*curl(v)*dx + 1e-6*nu*u*v*dx

f = LinearForm(V)
f += w * v * dx("coil")
```

Da das Problem abhängig von Mesh Grösse und Polynomordnung der FEM Basisfunktionen schon recht gross werden kann, benutzen wir einen iterativen Solver:

```{code-cell} ipython3
c = Preconditioner(a, type="multigrid")
```

Nun berechnen wir die Systemmatrix, den rechten Vektor und die Lösung des Systems parallel mit shared Memory:

```{code-cell} ipython3
with TaskManager():
    a.Assemble()
    f.Assemble()
    solver = CGSolver(mat=a.mat, pre=c.mat)
    gfu.vec.data = solver * f.vec
```

```{code-cell} ipython3
Draw (curl(gfu), mesh, "B-field", draw_surf=False, \
      clipping = { "pnt" : (0,0,0), "vec" : (0,1,0), "function" : False },
      vectors = { "grid_size" : 100 })
```

```{figure} bfield.png
---
align: left
height: 250px
name: Magfeld
---
B-Feld
```
