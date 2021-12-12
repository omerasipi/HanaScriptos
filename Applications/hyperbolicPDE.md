---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Hyperbolische Gleichungen zweiter Ordnung: Wellen Gleichung

In der Struktur Dynamik und auch anderen Anwendungen trifft man oft auf Systeme von der semidiskreten Form

$$M \ddot{d} + C \dot{d} + K d = F,$$ (eq:SemidiscrethyperbolicPDE)

wobei $M$ die Massenmatrix, $C$ eine viskose Dämpfung, $K$ die Steifigkeitsmatrix und $F$ ein Vektor gegeben durch die angewandten Kräfte ist. Wir nehmen an, dass $M, C, K$ symmetrisch sind. Die Massenmatrix $M$ sei positiv definit und $C, K$ positiv-semidefinit. Für die vollständige Beschreibung des Anfangswertproblems {eq}`eq:SemidiscrethyperbolicPDE` benötigen wir Anfangswerte für

$$d(0) = d_0\quad \text{und}\quad \dot{d}(0) = v_0.$$

Die wohl am meisten benutzten Verfahren gehöhren in die Familie der **Newmark Verfahren**:

$$\begin{split}
M\,a_{n+1} + C\,v_{n+1} + K\,d_{n+1} & = F_{n+1}\\
d_{n+1} & = d_n + \Delta t\,v_n + \frac{\Delta t^2}{2} ((1-2\beta) a_n + 2\beta a_{n+1})\\
v_{n+1} & = v_n + \Delta t\,((1-\gamma)\,a_n + \gamma\,a_{n+1}),
\end{split}$$ (eq:NewmarkSchema)

wobei $d_n, v_n, a_n$ Approximationen von $d(t_n), \dot{d}(t_n), \ddot{d}(t_n)$ sind. Die Parameter $\beta$ und $\gamma$ bestimmen die Stabilitäts- und Genauigkeits-Charakteristik. 


Es gibt verschiedene Formen für die Implementierung. Oft wird die sogenannte $a$-Form benutzt, welche auf der Beschleunigung $a$ aufbaut. Analoges mit der Deformation $d$ führt zu äquivalenten Implementationen. Wir betrachten hier die **a**-Form und definieren folgende **Predictors**:

$$\begin{split}
\tilde{d}_{n+1} & = d_n + \Delta t\,v_n + \frac{\Delta t^2}{2} (1-2\beta) a_n\\
\tilde{v}_{n+1} & = v_n + (1-\gamma) \Delta t\, a_n
\end{split}$$ (eq:NewmarkPredicators)

und für $d_{n+1}$ und $v_{n+1}$ aus {eq}`eq:NewmarkSchema` folgend die **Correctors**:

$$\begin{split}
d_{n+1} & = \tilde{d}_{n+1} + \beta\,\Delta t^2\, a_{n+1}\\
v_{n+1} & = \tilde{v}_{n+1} + \gamma\, \Delta t\, a_{n+1}.\end{split}$$ (eq:NewmarkCorrectors)

Die initiale Beschleunigung $a_0$ kann aus

$$M a_0 = F - C v_0 - K d_0$$

berechnet werden. Die noch unbekannte Beschleunigung zum neuen Zeitschritt $a_{n+1}$ ist nun gegeben durch die Lösung des Gleichungssystems

$$(M + \gamma\, \Delta t\, C + \beta\, \Delta t^2\, K)\,a_{n+1} = F_{n+1} - C\,\tilde{v}_{n+1} - K\,\tilde{d}_{n+1}.$$

Mit den Gleichungen {eq}`eq:NewmarkCorrectors` folgt $d_{n+1}$ und $v_{n+1}$, die Schätzwerte $\tilde{d}_{n+1}, \tilde{v}_{n+1}$ werden korrigiert. Entsprechend nennt man die Gleichungen {eq}`eq:NewmarkCorrectors` **Correctors**.

```{prf:remark}
:label: my-rm-trapezmethode

Für $\beta = \frac{1}{4}$ und $\gamma = \frac{1}{2}$ wird quasi die mittlere Beschleunigung benutzt. Wir erhalten in dem Fall die **Trapezmethode**, welche bedingungslos stabil ist.

```

## Anwendung

Wir betrachten als Anwendung die Wellengleichung

$$\begin{split}\ddot{u}(t,x) - \Delta u(t,x) & = f(t)\quad\text{für}\ x\in\Omega = [0,1]^2\\
u(t, x) & = 0\quad \text{für}\ x\in\partial \Omega\\
u(0) & = 0\quad \text{für}\ x\in \Omega\\
\dot{u}(0) & = 0\quad \text{für}\ x\in \Omega\\
\end{split}$$

```{code-cell} ipython3
:tags: [hide-cell, remove-output]

from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.webgui import Draw

import numpy as np
import matplotlib.pyplot as plt
```

```{code-cell} ipython3
from netgen.geom2d import CSG2d, Circle
circle = Circle((0,0),1,bc='wall')
geo = CSG2d()
geo.Add(circle)
```

```{code-cell} ipython3
#mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
mesh = Mesh(geo.GenerateMesh(maxh=0.1))# etwas grob 0.05
V = H1(mesh,order = 2, dirichlet='wall')#'bottom|right|top|left')
u,v = V.TnT()
gfu = GridFunction(V)
gfv = GridFunction(V)# 1. Zeitableitung
gfa = GridFunction(V)# 2. Zeitableitung
```

```{code-cell} ipython3
K = BilinearForm(V)
K += 0.1*grad(u)*grad(v)*dx
K.Assemble()

M = BilinearForm(V)
M += u*v*dx
M.Assemble();
```

Als Quelle benutzen wir einen Puls von der Zeitlänge $\Delta t$ mit einer beschränkten örtlichen Ausdehnung:

$$f(t) = e^{-((x-x_0)^2+(y-y_0)^2)/\tau^2}$$

```{code-cell} ipython3
x0=0
y0=0
f = 5000*exp(-((x-x0)**2+(y-y0)**2)/0.001)

F = LinearForm(V)
F += f*v*dx
F.Assemble();
```

Sei $M^* = M + \beta\, \Delta t^2\, K$. Wir wählen die Trapezmethode:

```{code-cell} ipython3
beta = 1/4
gamma = 1/2

dt = 1/125#0
mstar = K.mat.CreateMatrix()
mstar.AsVector().data = M.mat.AsVector() + dt * K.mat.AsVector()
```

Wir berechnen nun die zeitabhängige Lösung.

```{code-cell} ipython3
scene = Draw(gfu,mesh,'u', min=-2,max=2,autoscale=False)
```

```{code-cell} ipython3
t = 0
gfu.vec[:] = 0
gfv.vec[:] = 0
gfa.vec[:] = 0
T = 10
```

```{code-cell} ipython3
# predictors
utilde = gfu.vec.CreateVector()
vtilde = gfu.vec.CreateVector()
b = gfu.vec.CreateVector()
mip = mesh(0,0)
ui = [gfu(mip)]
ti = [t]
while t < T:
    #F.Assemble() # könnte hier vermieden werden!
    # predictors berechnne
    utilde.data = gfu.vec + dt * gfv.vec + dt**2/2*(1-2*beta)*gfa.vec
    vtilde.data = gfv.vec + (1-gamma)*dt*gfa.vec
    # a neu berechnen
    if t < dt: # ein Impuls
        b.data = F.vec - K.mat*utilde
    else:
        b.data = - K.mat*utilde
    gfa.vec.data = mstar.Inverse(freedofs=V.FreeDofs())*b
    # corrector berechnen
    gfu.vec.data = utilde + beta*dt**2*gfa.vec
    gfv.vec.data = vtilde + gamma*dt*gfa.vec
    # Redraw
    scene.Redraw()
    ui.append(gfu(mip))
    ti.append(t)
    t+=dt
```
