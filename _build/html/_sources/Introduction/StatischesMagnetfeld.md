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

## Geometrie

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

```{code-cell} ipython3

```
