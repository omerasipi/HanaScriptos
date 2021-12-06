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

# Parabolische partielle Differentialgleichungen

Wir betrachten die parabolische Differentialgleichung

$$\begin{split}\partial_t u(t,x) - \operatorname{div}(k(x)\, \nabla u(t,x)) & = q(t,x)\quad \text{für}\ x\in \Omega\\
u(t,x) & = u_D(x)\quad\text{für}\ x\in \Gamma_D\\
k(x)\,\partial_n u(t,x) & = g(t,x)\quad\text{für}\ x\in \Gamma_N\\
u(0,x) & = u_0(t,x)\quad\text{für}\ x\in \Omega
\end{split}$$ (eq:parabolischeDGLApplication)

mit den gegebenen Rand $\partial\Omega = \Gamma_D \cup \Gamma_N$. In dem Fall muss $\Gamma_D\not=\emptyset$ nicht gefordert werden.

Als Anwendung der parabolischen Gleichung betrachten wir den Temperaturverlauf in einem Raum (2d Querschnitt) mit einem Radiator, wobei wir für Raumseite beim Radiator eine Temperatur von 5°C als Dirichletrandwert vorgeben. Die restlichen Seiten des Raumes seien optimal isoliert, was wir mit der Neumann Randbedingung beschreiben können.

```{code-cell} ipython3
:tags: [hide-cell, remove-output]

from netgen.geom2d import CSG2d, Circle, Rectangle
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import matplotlib.pyplot as plt
```

```{code-cell} ipython3
geo = CSG2d()
raum = Rectangle( pmin=(0,0), pmax=(4,2.5), mat="raum", bc="bc_raum", right='right')
radiator = Rectangle( pmin=(3.8,0.2), pmax=(3.9,1), mat="radiator", bc="bc_radiator" )
air = raum-radiator

# add top level objects to geometry
geo.Add(air)
geo.Add(radiator)

mesh = Mesh( geo.GenerateMesh(maxh=0.2))
mesh.Curve(3) # Polynome 3. Grad für gekrümmte Geometrie
Draw(mesh);
```

Die schwache Gleichung für {eq}`eq:parabolischeDGLApplication` lautet

$$\int_\Omega \partial_t u\,v\,dx + \int_\Omega k\,\nabla u\cdot\nabla v\,dx = \int_\Omega f\,v\,dx\quad \forall \ v \in H^1(\Omega).$$

In einem ersten Schritt betrachten wir die stationäre Lösung, daher die Lösung, für welche wir keine Zeitabhängigkeit (mehr) haben:

$$\int_\Omega k\,\nabla u\cdot\nabla v\,dx = \int_\Omega f\,v\,dx\quad \forall \ v \in H^1(\Omega)$$

Die Lösung unabhängig von der Zeit liegt im Sobolev-Raum $H^1(\Omega)$.

Im zweiten Schritt suchen wir die zeitabhängige Lösung, daher die Funktion $u: [0,T] \to H^1(\Omega)$ mit

$$\int_\Omega\partial_t u\,v\,dx + \int_\Omega k\,\nabla u\cdot\nabla v\,dx = \int_\Omega f\,v\,dx\quad \forall \ v \in H_0^1(\Omega),$$ (eq:weakparabolischeDGLApplication)

wobei der Anfangswert $u(t=0) = u_0(x)$ gegeben ist.


## Stationäre Lösung


Randbedingungen

```{code-cell} ipython3
mesh.GetBoundaries()
```

```{code-cell} ipython3
V = H1(mesh, order=3, dirichlet="right")
u,v = V.TnT()
```

Für den Radiator nehmen wir als Wärmequelle an:

```{code-cell} ipython3
f = LinearForm(V)
f += 100*v*dx(definedon=mesh.Materials('radiator'))
f.Assemble();
```

Für die Wärmeleitung haben wir abhängig vom Gebiet stückweise unterschiedliche Werte:

```{code-cell} ipython3
mesh.GetMaterials()
```

```{code-cell} ipython3
# Luft, Wasser in der Reihenfolge der subdomains:
k = CoefficientFunction([0.0262,0.5562])
```

Damit können wir die Bilinearform definieren

```{code-cell} ipython3
a = BilinearForm(V)
a += k*grad(u)*grad(v)*dx
a.Assemble();
```

Randbedingung an der Rechten Seite des Raums ist die Front auf 5°C, die restlichen Seiten sind optimal isoliert. Ist der Radiator ausgeschaltet, sollte die Temperatur im Raum 5° konstant betragen:

```{code-cell} ipython3
gfu_D = GridFunction(V)
gfu_D.Set(CoefficientFunction(5))
```

```{code-cell} ipython3
gfu_stat = GridFunction(V)
```

Die Lösung ist gegeben durch $u = u_0 + u_D$, wobei wir das $u_0$ für gegebenes $u_D$ berechnen.

```{code-cell} ipython3
gfu_stat.vec.data = gfu_D.vec
gfu_stat.vec.data += a.mat.Inverse(freedofs=V.FreeDofs())*(f.vec-a.mat*gfu_D.vec)
Draw(gfu_stat,mesh,'T stationär');
```

Für die stationäre Temperatur in der Mitte des Raumes erhalten wir

```{code-cell} ipython3
mip = mesh(4/2,2.5/2)
Tstat = gfu_stat(mip)
Tstat
```

Die mittlere Temperatur in der Luft beträgt

```{code-cell} ipython3
Aair = Integrate(1, mesh, definedon = mesh.Materials('raum'))
TstatMean = 1/Aair*Integrate(gfu_stat, mesh, definedon = mesh.Materials('raum'))
TstatMean
```

## Zeitabhängige Lösung


Mit Hilfe der stationären Gleichung erhalten wir die Temperaturverteilung nach langer (genau genommen unendlich langer) Zeit. Nun insteressiert uns der zeitliche Verlauf. Wir lösen daher das parabolische Problem.

Wir suchen also $u(t)$, die Funktion, welche uns für jedes Zeit die örtliche Verteilung der Temperatur liefert, sprich $u(t) \in H^1(\Omega)$.

Dazu ersetzen wir die partielle Ableitung nach der Zeit durch den Differenzenquotient:

$$\partial_t u(t) \approx \frac{u(t+\Delta t) - u(t)}{\Delta t}.$$

Es folgt für {eq}`eq:weakparabolischeDGLApplication`

$$\int_\Omega \frac{u(t+\Delta)-u(t)}{\Delta t}\,v\,dx + \int_\Omega k\,\nabla u\cdot\nabla v\,dx = \int_\Omega f\,v\,dx\quad \forall \ v \in H_0^1(\Omega).$$

Es stellt sich die Frage, für welches $t$ wir $\nabla u$ auswerten.

* Für $\nabla u(t)$ folgt das **explizite** Euler Verfahren. Durch Separation der bekannten Grössen auf die rechte Seite und der unbekannten auf die linke Seite folgt

  $$\int_\Omega u(t+\Delta)\,v\,dx = \int_\Omega u(t) v dx + \Delta t \left(\int_\Omega f(t+\Delta t)\,v\,dx - \int_\Omega k\,\nabla u(t)\cdot\nabla v\,dx\right)\quad \forall \ v \in H_0^1(\Omega)$$

  In dem Fall können wir durch einfache Vektor Addition den zeitlichen Verlauf berechnen. Das Verfahren ist jedoch für grosse $\Delta t$ nicht stabil. Wir benutzen daher diesen Ansatz nicht.
* Für $\nabla u(t+\Delta t)$ folgt das **implizite** Euler Verfahren

  $$\int_\Omega (u(t+\Delta)-u(t))\,v\,dx + \Delta t \int_\Omega k\,\nabla u(t+\Delta t)\cdot\nabla v\,dx = \Delta t \int_\Omega f\,v\,dx\quad \forall \ v \in H_0^1(\Omega).$$ (eq:weakparabolischeDGLApplicationImplizit)

  Separieren wir wiederum bekanntes und unbekanntes auf die rechte bzw. linke Seite der Gleichung so folgt das Gleichungssystem
  
  $$\int_\Omega u(t+\Delta)\,v\,dx + \Delta t \int_\Omega k\,\nabla u(t+\Delta t)\cdot\nabla v\,dx = \int_\Omega u(t)\,v\,dx + \Delta t \int_\Omega f\,v\,dx\quad \forall \ v \in H_0^1(\Omega),$$

  welches wir in jedem Zeitschritt lösen müssen. Diskret können wir das System auch als
  
  $$\underbrace{(M + \Delta t\ A)}_{=:M^*}\cdot u_{n+1} = M\cdot u_n + f$$
  
  mit der Massenmatrix $M$ für die Bilinearform $\int u\,v dx$ und der Steiffigkeitsmatrix $A$ für die Bilinearform $\int \nabla u\cdot \nabla v\, dx$.

```{code-cell} ipython3
m = BilinearForm(V)
m += u*v*dx
m.Assemble()
```

Für das $\Delta t$ setzen wir

```{code-cell} ipython3
time = 0
dt = 1
```

Damit folgt für $M^*$

```{code-cell} ipython3
mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
invmstar = mstar.Inverse(freedofs=V.FreeDofs())
```

Der Anfangswert ist gegeben durch

```{code-cell} ipython3
gfu = GridFunction(V)
gfu.vec.data = gfu_D.vec
```

```{code-cell} ipython3
scene = Draw(gfu,mesh,"u");
```

```{code-cell} ipython3
res = gfu.vec.CreateVector()
tstep = 600 # time that we want to step over within one block-run
```

$M^*\cdot(u_{0,n+1}+u_D) = M\cdot u_n + f$

```{code-cell} ipython3
t_intermediate=0 # time counter within one block-run
while t_intermediate < tstep - 0.5 * dt:
    res.data = m.mat*gfu.vec + dt * f.vec - mstar*gfu_D.vec
    gfu.vec.data = gfu_D.vec
    gfu.vec.data += invmstar * res 
    t_intermediate += dt
    scene.Redraw()
time+=t_intermediate
```

Oft benutzt man auch eine inkrementelle Berechnung. Mit $\delta u = u(t+\Delta t)-u(t)$ folgt für {eq}`eq:weakparabolischeDGLApplicationImplizit`

$$\int_\Omega \delta u\,v\,dx + \Delta t \int_\Omega k\,\nabla \delta u\cdot\nabla v\,dx + \Delta t \int_\Omega k\,\nabla u(t)\cdot\nabla v\,dx = \Delta t \int_\Omega f\,v\,dx\quad \forall \ v \in H_0^1(\Omega).$$

und damit

$$\int_\Omega \delta u\,v\,dx + \Delta t \int_\Omega k\,\nabla \delta u\cdot\nabla v\,dx = \Delta t \int_\Omega f\,v\,dx - \Delta t \int_\Omega k\,\nabla u(t)\cdot\nabla v\,dx\quad \forall \ v \in H_0^1(\Omega).$$

Für das diskrete System folgt
  
$$\underbrace{(M + \Delta t\ A)}_{=:M^*}\cdot \delta u_{n+1} = f - A\cdot u_n.$$

Die Form hat den entscheidenden Vorteil, dass die Dirichletrandwerte nicht in jedem Zeitschritt berechnet werden müssen. Wir sparen uns daher eine Matrix Multiplikation und Vektor Initialisierung.

```{code-cell} ipython3
time = 0
gfu2 = GridFunction(V)
gfu2.vec.data = gfu_D.vec
```

```{code-cell} ipython3
scene2 = Draw(gfu2,mesh,'u');
```

```{code-cell} ipython3
t_intermediate=0 # time counter within one block-run
# Temperatur in Raummitte
T = [gfu2(mip)]
Tmean = [1/Aair*Integrate(gfu2, mesh, definedon = mesh.Materials('raum'))]
while t_intermediate < tstep - 0.5 * dt:
    res.data = dt * f.vec - dt * a.mat * gfu2.vec
    gfu2.vec.data += invmstar * res
    t_intermediate += dt
    T.append(gfu2(mip))
    Tmean.append(1/Aair*Integrate(gfu2, mesh, definedon = mesh.Materials('raum')))
    scene2.Redraw()
time+=t_intermediate
```

```{code-cell} ipython3
Draw(gfu-gfu_stat,mesh,'Differenz');
```

## Erweiterung mit Konvektion

Wir erweitern das parabolische Problem noch durch einen konvektiven Term. Daher eine Strömung, welche die Temperatur im Raum verteilt.

$$\begin{split}\partial_t u(t,x) - \operatorname{div}(k(x)\, u(t,x)) + b\cdot \nabla u & = q(t,x)\quad \text{für}\ x\in \Omega\\
u(t,x) & = u_D\quad\text{für}\ x\in \partial\Omega
\end{split}$$

Das Strömungsfeld $b$ nennt man wind, wir definieren dieses gegeben als

$$b(x,y) = v_0 \begin{pmatrix}
-\frac{(y-1.25)}{1.25}\, \left(1-\left(\frac{x}{4}\right)\,\frac{x}{4}\right)\\
\frac{(x-2)}{2} \frac{y}{2.5}\,\left(1-\left(\frac{y}{2.5}\right)\right)\end{pmatrix}$$

```{code-cell} ipython3
v0 = 2
b = v0*CoefficientFunction((-(y-1.25)/1.25*(1-(x/4))*x/4, (x-2)/2*y/2.5*(1-y/2.5)))
Draw(b,mesh,"wind");
from ngsolve.internal import visoptions
visoptions.scalfunction = "wind:0"
```

```{code-cell} ipython3
a = BilinearForm(V, symmetric=False)
a += k*grad(u)*grad(v)*dx + b*grad(u)*v*dx
a.Assemble();
```

Für die stationäre Lösung folgt

```{code-cell} ipython3
gfu_stat2 = GridFunction(V)
gfu_stat2.vec.data = gfu_D.vec
gfu_stat2.vec.data += a.mat.Inverse(freedofs=V.FreeDofs())*(f.vec - a.mat*gfu_D.vec)
```

Wie beobachtet werden kann ist die Lösung unphysikalisch. Wird $v_0$ zu gross gewählt, wird das numerische Problem instabil. In dem Fall müssen wir entweder
* die Diskretisierung stabilisieren
* oder das Mesh feiner wählen.

```{code-cell} ipython3
Draw(gfu_stat2,mesh,'u2 stationär');
```

```{code-cell} ipython3
Tstat2 = gfu_stat2(mip)
print('stationäre Temperatur in Raum Mitte: ', Tstat2)

Tstat2Mean = 1/Aair*Integrate(gfu_stat2, mesh, definedon = mesh.Materials('raum'))
print('stationäre mittlere Temperatur: ', Tstat2)
```

```{code-cell} ipython3
mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
invmstar = mstar.Inverse(freedofs=V.FreeDofs())
```

```{code-cell} ipython3
gfu3 = GridFunction(V)
```

```{code-cell} ipython3
time = 0
gfu3.vec.data = gfu_D.vec
```

```{code-cell} ipython3
scene3 = Draw(gfu3,mesh,'u');
```

```{code-cell} ipython3
t_intermediate=0 # time counter within one block-run
# Temperatur in Raummitte
Tkonv = [gfu3(mip)]
TkonvMean = [1/Aair*Integrate(gfu3, mesh, definedon = mesh.Materials('raum'))]
while t_intermediate < tstep - 0.5 * dt:
    res.data = dt * f.vec - dt * a.mat * gfu3.vec
    gfu3.vec.data += invmstar * res
    t_intermediate += dt
    Tkonv.append(gfu3(mip))
    TkonvMean.append(1/Aair*Integrate(gfu3, mesh, definedon = mesh.Materials('raum')))
    scene3.Redraw()
time+=t_intermediate
```

Im Bereich vor dem Radiator ist es bedeutend kühler, hingegen durch den Abtransport der Wärme oberhalb bedeutend wärmer.

```{code-cell} ipython3
Draw(gfu3-gfu2,mesh,'Differenz');
```

Für die zeitliche Entwicklung der Tempertur der beiden Ansätze erhalten den folgenden Verlauf

```{code-cell} ipython3
:tags: [hide-input]
plt.axhline(Tstat,c='gray',label='stationäre Temperatur ohne Konvektion')
plt.axhline(Tstat2,c='gray',label='stationäre Temperatur mit Konvektion')
plt.plot(dt*np.arange(len(T)),T,label='Temperatur ohne Konvektion')
plt.plot(dt*np.arange(len(Tkonv)),Tkonv,label='Temperatur mit Konvektion')
plt.legend()
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('T [°]')
plt.title('zeitlicher Temperatur Verlauf in Raum Mitte')
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
plt.axhline(TstatMean,c='gray',label='stationäre Temperatur ohne Konvektion')
plt.axhline(Tstat2Mean,c='gray',label='stationäre Temperatur mit Konvektion')
plt.plot(dt*np.arange(len(Tmean)),Tmean,label='Temperatur ohne Konvektion')
plt.plot(dt*np.arange(len(TkonvMean)),TkonvMean,label='Temperatur mit Konvektion')
plt.legend()
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('T [°]')
plt.title('zeitlicher Verlauf der mittleren Temperatur')
plt.show()
```
