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

# Einstiegsbeispiel: Espresso in Isolations- und Keramiktasse

## Modellierung

Wir betrachten die parabolische Differentialgleichung

$$\begin{split}C_p(x) \rho(x)\, \partial_t T(t,x) - \operatorname{div}(k(x)\, \nabla T(t,x)) & = 0\quad \text{für}\ x\in \Omega\\
T(t,x) & = T_D(x)\quad\text{für}\ x\in \partial\Omega\\
T(0,x) & = T_0(t,x)\quad\text{für}\ x\in \Omega
\end{split}$$ (eq:parabolischeDGLApplicationNoSource)

mit den gegebenen Rand $\partial\Omega = \Gamma_D \cup \Gamma_N$. In dem Fall muss $\Gamma_D\not=\emptyset$ nicht gefordert werden. Die Materialparameter sind vom Ort abhängig und wie folgt definiert: $C_p$ ist die spezifische Wärmekapazität, $\rho$ die Dichte des Mediums und $k$ der Wärmeleitungskoeffizient.

Wir betrachten nun das in der Einführung zur Vorlesung betrachtete Beispiel des Espressos in zwei verschiedenen Tassen, wobei wir nun an der zeitlichen Entwicklung des Temperaturfeldes $T$ interessiert sind.

```{figure} Espressotasse.jpg
---
align: center
height: 200px
name: Kaffeetassen2
---
Espressotasse
```

Wir gehen davon aus, dass abgesehen vom Kaffee die Temperatur durch die Umgebungstemperatur $T_0=20^\circ$ C gegeben ist. Am Rand nehmen wir an, dass die Temperatur konstant bleibt. Der Kaffee wird sich daher abkühlen und für $t\to\infty$ gegen die Umgebungstemperatur konvergieren.

```{code-cell} ipython3
:tags: [hide-input]
from netgen.occ import WorkPlane, OCCGeometry, X, Y, Z, Glue, Compound
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import matplotlib.pyplot as plt

wp = WorkPlane()
sPnts = [(0.024, 0.008),
         (0.028, 0.016),
         (0.032, 0.04),
         (0.022, 0.02),
         (0.004, 0.006),
         (-0.004, 0.006),
         (-0.022, 0.02),
         (-0.032, 0.04),
         (-0.028, 0.016),
         (-0.024, 0.008),
         (-0.016, 0)]
wp.MoveTo(4*4e-3,0).Spline(sPnts)
face = wp.Close().Face()

wp2 = WorkPlane()
sPnts2 = [(0.0224, 0.008),
          (0.0264, 0.016),
          (0.0312, 0.036),
          (0.022, 0.0184),
          (0.004, 0.0048),
          (-0.004, 0.0048),
          (-0.022, 0.0184),
          (-0.0312, 0.036),
          (-0.0264, 0.016),
          (-0.0224, 0.008),
          (-0.0152, 0.0012)]
wp2.MoveTo(4*3.8e-3,4*.3e-3).Spline(sPnts2)
isolator = wp2.Close().Face()

glas = face-isolator

wp3 = WorkPlane()
wp3.MoveTo(-4*3.9e-3,2*1.5e-3).LineTo(4*3.9e-3,2*1.5e-3).LineTo(4*7.2e-3,3*9.5e-3).LineTo(-4*7.2e-3,3*9.5e-3)
face3 = wp3.Close().Face()
liquid = face3-face

wp4 = WorkPlane()
wp4.MoveTo(-.1,0).LineTo(.1,0).LineTo(.1,.15).LineTo(-.1,.15)
box = wp4.Close().Face()
air = box-face-face3

air.faces.name = 'air'
glas.faces.name = 'glas'
isolator.faces.name = 'isolator'
liquid.faces.name = 'liquid'

air.maxh = 30e-3
isolator.maxh = 5e-3
glas.maxh = 1e-3
liquid.maxh = 2.5e-3

glas.edges.Min(Y).name='bottom'
air.edges.Min(X).name='outer'
air.edges.Max(X).name='outer'
air.edges.Max(Z).name='outer'
air.edges.Max(Y).name='outer'
air.edges.Nearest((-.05,0)).name='bottom'
air.edges.Nearest((.05,0)).name='bottom'

model = Glue([air, glas, isolator, liquid])

geo = OCCGeometry(model, dim=2) # wichtig: dim=2 muss hier zwingend definiert werden!
mesh = Mesh(geo.GenerateMesh(maxh=0.015))

Draw(mesh);
```

Das initiale Temperaturfeld $T_0$ ist gegeben durch

$$T_0(x) = \begin{cases}80 & \quad \text{für}\ x\in\Omega_{\text{liquid}}\\
20  & \quad \text{sonst.}\end{cases}$$

Wir benötigen daher einen $H^1$-finite Elementraum, wobei wir den Rand als Dirichlet Rand mit $T(\cdot,x) = T_D = 20^\circ$ C definieren.

```{code-cell} ipython3
V = H1(mesh, order=2, dirichlet="bottom|outer")
T = V.TrialFunction()
v = V.TestFunction()
print('V.ndof = ',V.ndof)
```

Die Materialparameter sind vom Gebiet abhängig und gegeben durch:

```{code-cell} ipython3
# Isolationstasse
k_iso_mat = {'glas': 20, 'isolator': 0.0262, 'liquid': 0.597, 'air': 0.0262}
# Keramiktasse
k_ker_mat = {'glas': 20, 'isolator': 20, 'liquid': 0.597, 'air': 0.0262}
k_iso = CoefficientFunction([k_iso_mat[mat] for mat in mesh.GetMaterials()])
```

Mit einer `mesh.MaterialCF` kann dies noch einfacher definiert werden:

```{code-cell} ipython3
k_iso = mesh.MaterialCF(k_iso_mat)
k_ker = mesh.MaterialCF(k_ker_mat)
```

Für die Wärmekapazität und Dichte wählen wir

```{code-cell} ipython3
# Isolationstasse
rho_iso_mat = {'glas': 2500, 'isolator': 1.25, 'liquid': 998, 'air': 1.25}
Cp_iso_mat = {'glas': 800, 'isolator': 1010, 'liquid': 4184, 'air': 1010}
rho_iso = mesh.MaterialCF(rho_iso_mat)
Cp_iso = mesh.MaterialCF(Cp_iso_mat)
# Keramiktasse
rho_ker_mat = {'glas': 2500, 'isolator': 2500, 'liquid': 998, 'air': 1.25}
Cp_ker_mat = {'glas': 800, 'isolator': 800, 'liquid': 4184, 'air': 1010}
rho_ker = mesh.MaterialCF(rho_ker_mat)
Cp_ker = mesh.MaterialCF(Cp_ker_mat)
```

Die initiale Temperatur $T_0$ definieren wir analog:

```{code-cell} ipython3
T0_mat = {'glas': 20, 'isolator': 20, 'liquid': 80, 'air': 20}
T0 = mesh.MaterialCF(T0_mat)

Draw(T0,mesh,'T_0');
```

Da wir die Temperatur $T_D$ in jedem Zeitschritt benötigen, definieren wir diese als `GridFunction` und definieren am Rand $20^\circ$ C.

```{code-cell} ipython3
gfTD = GridFunction(V)
gfTD.Set(20,BND)
```

Für die Temperatur ca. in der Mitte des Espressos erhalten wir

```{code-cell} ipython3
mip = mesh(0,.02)
mip2 = mesh(0,.04)
Tstat = T0(mip)
Tstat
```

Die mittlere Temperatur im Kaffee ist gegeben durch

```{code-cell} ipython3
Aliquid = Integrate(1, mesh, definedon = mesh.Materials('liquid'))
TstatMean = 1/Aliquid*Integrate(T0, mesh, definedon = mesh.Materials('liquid'))
TstatMean
```

## Zeitabhängige Lösung

Wir betrachten nun das parabolische Problem und suchen $T(t)$, die Funktion, welche uns für jedes Zeit die örtliche Verteilung der Temperatur liefert, sprich $T(t) \in H^1(\Omega)$.

$$\int_\Omega C_p\, \rho\,(\partial_t T)\,v\,dx + \int_\Omega k\,\nabla T\cdot\nabla v\,dx = 0\quad \forall \ v \in H^1(\Omega).$$ (eq:weakparabolischeDGLApplicationNoSource)

Dazu ersetzen wir die partielle Ableitung nach der Zeit durch den **Differenzenquotient**:

$$\partial_t T(t) \approx \frac{T(t+\Delta t) - T(t)}{\Delta t}.$$

Es folgt für {eq}`eq:weakparabolischeDGLApplicationNoSource`

$$\int_\Omega \frac{T(t+\Delta)-T(t)}{\Delta t}\,v\,dx + \int_\Omega k\,\nabla T\cdot\nabla v\,dx = 0\quad \forall \ v \in H_0^1(\Omega).$$

Es stellt sich die Frage, für welches $t$ wir $\nabla T$ auswerten.

* Für $\nabla T(t)$ folgt das **explizite** Euler Verfahren. Durch Separation der bekannten Grössen auf die rechte Seite und der unbekannten auf die linke Seite folgt

  $$\int_\Omega T(t+\Delta t)\,v\,dx = \int_\Omega T(t) v dx - \Delta t \int_\Omega k\,\nabla u(t)\cdot\nabla v\,dx\quad \forall \ v \in H_0^1(\Omega)$$

  In dem Fall können wir durch einfache Vektor Addition den zeitlichen Verlauf berechnen. Das Verfahren ist jedoch für grosse $\Delta t$ nicht stabil. Wir benutzen daher diesen Ansatz nicht.
* Für $\nabla T(t+\Delta t)$ folgt das **implizite** Euler Verfahren

  $$\int_\Omega (T(t+\Delta t)-T(t))\,v\,dx + \Delta t \int_\Omega k\,\nabla T(t+\Delta t)\cdot\nabla v\,dx = 0\quad \forall \ v \in H_0^1(\Omega).$$ (eq:weakparabolischeDGLApplicationImplizitNoSource)

  Separieren wir wiederum bekanntes und unbekanntes auf die rechte bzw. linke Seite der Gleichung so folgt das Gleichungssystem
  
  $$\int_\Omega T(t+\Delta t)\,v\,dx + \Delta t \int_\Omega k\,\nabla T(t+\Delta t)\cdot\nabla v\,dx = \int_\Omega T(t)\,v\,dx\quad \forall \ v \in H_0^1(\Omega),$$

  welches wir in jedem Zeitschritt lösen müssen. Diskret können wir das System auch als
  
  $$\underbrace{(M + \Delta t\ A)}_{=:M^*}\cdot T_{n+1} = M\cdot T_n$$ (eq:eulerimplizitParabolischNoSource)
  
  mit der Massenmatrix $M$ für die Bilinearform $\int T\,v dx$ und der Steiffigkeitsmatrix $A$ für die Bilinearform $\int \nabla T\cdot \nabla v\, dx$.

```{code-cell} ipython3
a_iso = BilinearForm(V)
a_iso += k_iso*grad(T)*grad(v)*dx
a_iso.Assemble()

m_iso = BilinearForm(V)
m_iso += rho_iso*Cp_iso*T*v*dx
m_iso.Assemble();
```

Für das $\Delta t$ setzen wir

```{code-cell} ipython3
time = 0
dt = 5
```

Damit folgt für $M^*$

```{code-cell} ipython3
mstar_iso = m_iso.mat.CreateMatrix()
mstar_iso.AsVector().data = m_iso.mat.AsVector() + dt * a_iso.mat.AsVector()
invmstar_iso = mstar_iso.Inverse(freedofs=V.FreeDofs())
```

Der Anfangswert ist gegeben durch die `CoefficientFunction T0`

```{code-cell} ipython3
gfT_iso = GridFunction(V)
gfT_iso.Set(T0)
gfT_iso.AddMultiDimComponent(gfT_iso.vec)
```

Um die inhomogene Dirichlet-Randbedingung zu erfüllen, benutzen wir wieder den Ansatz

$$T_{n+1}(x) = T_{0,n+1}(x) + T_D(x)$$

mit $T_0 \in H_0^1(\Omega)$. Für ({eq}`eq:eulerimplizitParabolischNoSource`) folgt

$$M^*\cdot(T_{0,n+1}+T_D) = M\cdot T_n$$

und damit für jeden Zeitschritt

1. $M^*\cdot T_{0,n+1} = M\cdot T_n - M^*\cdot T_D$
2. $T_{n+1}(x) = T_{0,n+1}(x) + T_D$

```{code-cell} ipython3
res = gfT_iso.vec.CreateVector()
tstep = 60*60

t_intermediate=0 # time counter within one block-run
while t_intermediate < tstep - 0.5 * dt:
    res.data = m_iso.mat*gfT_iso.vec - mstar_iso*gfTD.vec
    gfT_iso.vec.data = gfTD.vec
    gfT_iso.vec.data += invmstar_iso * res
    if t_intermediate % 50 == 0:
        gfT_iso.AddMultiDimComponent(gfT_iso.vec)
    t_intermediate += dt
time+=t_intermediate
```

```{code-cell} ipython3
Draw(gfT_iso,mesh,"T",min=20,max=80);
```

Oft benutzt man auch eine inkrementelle Berechnung. Mit $\delta T = T(t+\Delta t)-T(t)$ folgt für {eq}`eq:weakparabolischeDGLApplicationImplizitNoSource`

$$\int_\Omega \delta T\,v\,dx + \Delta t \int_\Omega k\,\nabla \delta T\cdot\nabla v\,dx + \Delta t \int_\Omega k\,\nabla T(t)\cdot\nabla v\,dx = 0\quad \forall \ v \in H_0^1(\Omega).$$

und damit

$$\int_\Omega \delta T\,v\,dx + \Delta t \int_\Omega k\,\nabla \delta T\cdot\nabla v\,dx = - \Delta t \int_\Omega k\,\nabla T(t)\cdot\nabla v\,dx\quad \forall \ v \in H_0^1(\Omega).$$

Für das diskrete System folgt
  
$$\underbrace{(M + \Delta t\ A)}_{=:M^*}\cdot \delta T_{n+1} = - \Delta t\, A\cdot T_n.$$

Die Form hat den entscheidenden Vorteil, dass die Dirichletrandwerte nicht in jedem Zeitschritt berechnet werden müssen. Wir sparen uns daher eine Matrix Multiplikation und Vektor Initialisierung.

```{code-cell} ipython3
time = 0
gfT2_iso = GridFunction(V)
gfT2_iso.Set(T0)

T_iso = [gfT2_iso(mip)]
T2_iso = [gfT2_iso(mip2)]
ti_iso = [time]
Tmean_iso = [1/Aliquid*Integrate(gfT2_iso, mesh, definedon = mesh.Materials('liquid'))]
```

```{code-cell} ipython3
scene2 = Draw(gfT2_iso,mesh,'T',min=20,max=80);
```

```{code-cell} ipython3
t_intermediate=0 # time counter within one block-run
while t_intermediate < tstep - 0.5 * dt:
    res.data = - dt * a_iso.mat * gfT2_iso.vec
    gfT2_iso.vec.data += invmstar_iso * res
    t_intermediate += dt
    time += dt
    # Speichern der Resultate
    ti_iso.append(time)
    T_iso.append(gfT2_iso(mip))
    T2_iso.append(gfT2_iso(mip2))
    Tmean_iso.append(1/Aliquid*Integrate(gfT2_iso, mesh, definedon = mesh.Materials('liquid')))
    scene2.Redraw()
```

```{code-cell} ipython3
:tags: [hide-input]
plt.plot(np.array(ti_iso)/60,T_iso, label='Temperatur im Punkt $p = (0,20)$  [mm]')
plt.plot(np.array(ti_iso)/60,T2_iso, label='Temperatur im Punkt $p_2=(0,40)$ [mm]')
plt.plot(np.array(ti_iso)/60,Tmean_iso, label='mittlere Temperatur')
plt.legend()
plt.xlabel('$t$ [Min]')
plt.ylabel('$T$')
plt.title('Temperaturverlauf') 
plt.grid()
plt.show()
```

```{admonition} Aufgabe
Erweitern Sie das Notebook mit der Berechnung des Temperaturverlaufs für die Keramiktasse und vergleichen Sie die Lösung.
```

+++

```{prf:remark}
Die Temperatur nimmt hier relativ langsam ab. In der Realität können wir auch mit einer Isolationstasse eine schnellere Abkühlung beobachten. Das bedeutet, dass die im Modell berücksichtigte Wärmeleitung nicht alleine für die Abkühlung verantwortlich ist. Die Konvektion müsste hier ebenfalls berücksichtigt werden. Man könnte das Modell mit Hilfe der Boussinesq Approximation erweitern.
```

```{code-cell} ipython3

```
