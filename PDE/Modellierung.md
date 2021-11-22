# Modellierung

## Elliptische partielle Differentialgleichungen

### Randbedingungen

Sei $\Omega \subset \mathbb{R}^d$ ein beschränktes Gebiet. Eine elliptische PDE benötigt eine Randbedingung für alle Punkte des Randes von $\Omega$. Mögliche Randbedingungen sind

* **Dirichlet-Randbedingung**: Bedingung an den Funktionswert

  $$u(x) = u_d(x)\quad \forall x\in\Gamma_D$$

* **Neumann-Randbedingung**: Bedingung an die Normalableitung

  $$\frac{\partial u}{\partial n}(x) = \nabla u\cdot n = g(x)\quad \forall x\in\Gamma_N$$

* **Robin-Randbedingung**: Gemischte Randbedingung

  $$\frac{\partial u}{\partial n}(x) + \alpha\, u(x) = g(x)\quad \forall x\in\Gamma_R$$

### Bilanzierung

Am Beispiel der stationären Wärmeleitung in einer dünnen Platte mit konstanter Dicke modellieren wir einen **stationären Diffusionsprozess**. Alle Grössen seien über die Dicke konstant. Deshalb reicht es die Mittelfläche $\Omega \subset \mathbb{R}^2$ zu betrachten. Die involvierten physikalischen Grössen sind

$$\begin{align}
T(x) & \quad\text{Temperatur [K]} & T_0(x) & \quad\text{Umgebungstemperatur [K]}\\
q(x) & \quad\text{Wärmefluss [J/(m s)]} & 
f(x) & \quad\text{Wärmequelle [J/(m$^2$ s)]}
\end{align}$$

Sei $L$ ein Linienstück mit beliebig orientiertem Normalvektor $n$, dann ist die Wärmemenge, die pro Zeiteinheit durch das Linienstück in die durch $n$ definierte Richtung fliesst, definiert durch

$$\int_L q(x)\cdot n\, dx.$$

Für den Wärmefluss $q$ (Vektor) nehmen wir das Fourier'sche Gesetz an

$$q(x) = -k\, \nabla T.$$

Der Energieabfluss in die Umgebung unterhalb und oberhalb der Platte sei durch ein lineares Wärmeübergangsgesetz (mit Wärmeübergngskoeffizienten $\alpha > 0$) beschrieben

$$r(x) = \alpha\, (T-T_0).$$

Sein nun $B \subset \Omega$ ein beliebiges Kontrollvolumen, dann lautet die Energiebilanz:

$$\boxed{\begin{array}{l}\text{Erzeugte Energie}\\\text{in $B$ pro Zeiteinheit}\end{array}} = \boxed{\begin{array}{l}\text{Energieabfluss aus $B$}\\\text{in Umgebung pro Zeiteinheit}\end{array}} + \boxed{\begin{array}{l}\text{Energieabfluss durch}\\\text{$\partial B$ pro Zeiteinheit}\end{array}}$$

also

$$\int_B f(x) dx = \int_B r(x) dx + \int_{\partial B} q(x)\cdot n\, ds.$$

Der Satz von Gauss erlaubt uns das Randintegral in ein Volumenintegral umzuschreiben

$$\int_{\partial B} q(x)\cdot n\, ds = \int_B \mathop{div} q(x) dx.$$

Damit folgt

$$\int_B f(x)\, dx = \int_B r(x) + \int_{B} \mathop{div} q(x)\,dx.$$

Da das Kontrollvolumen $B$ beliebig ist, erhält man unter der Voraussetzung, dass die Funktionen stetig sind

$$f(x) = r(x) +  \mathop{div} q(x).$$

Setzt man die Terme für $q(x)$ und $r(x)$ ein, so erhalten wir die partielle Differentialgleichung zweiter Ordnung

$$-\mathop{div}(k\,\nabla T(x)) + \alpha\, T(x) = f(x) + \alpha\, T_0(x)\quad \forall\,x\in\Omega.$$

Die möglichen Randbedingungen sind

* Dirichlet-Randbedingung

  $$T(x) = T_d(x)\quad \forall x\in\Gamma_D,$$

* Neumann-Randbedingung

  $$-k\,\frac{\partial T}{\partial n}(x) = g(x)\quad \forall x\in\Gamma_N,$$

* Robin-Randbedingung

  $$-k\,\frac{\partial T}{\partial n}(x) + \alpha\, T(x) = g(x)\quad \forall x\in\Gamma_R.$$

```{prf:remark}
:label: my-rm-mod1

Die hier vorgestellte Bilanzierung lässt sich analog in beliebig dimensionale Räume übertragen und gilt daher identisch auch im $\mathbb{R}^3$.
```

## Partielle Differentialgleichungen erster Ordnung

Im Zusammenhang mit Transportgleichungen treffen wir typischerweise partielle Differentialgleichungen erster Ordnung an. Gegeben sei die Strömung $b(x)$ und Quelle $f(x)$. Der Rand $\Gamma = \partial\Omega$ wird je nach Vorzeichen von $b\cdot n$ in Einflussrand $\Gamma_{\text{in}}$ oder Ausflussrand $\Gamma_{\text{out}}$ unterteilt:

$$\begin{split}
\Gamma_{\text{in}} & = \{x\in\Gamma\,|\, b\cdot n < 0\}\\
\Gamma_{\text{out}} & = \{x\in\Gamma\,|\, b\cdot n \ge 0\}.
\end{split}$$

Gesucht ist eine Konzentration $u(x)$, welche die Gleichung

$$\begin{split}
\mathop{div}(b\,u(x)) & = f(x)\quad \text{in}\,\Omega\\
u & = u_g \quad \text{in}\,\Gamma_{\text{in}}
\end{split}$$

erfüllt. Die Randbedingung kann nur am Einflussrand $\Gamma_{\text{in}}$ gestellt werden. Oft ist das gegebene Strömungsfeld divergenzfrei, daher $\mathop{div} b = 0$. In der Strömungsmechanik bedeutet das, dass die Strömung inkompressibel ist oder Volumen erhaltend.

Wir betrachten die Modellierung wiederum im zweidimensionalen. Die Transportgleichung beschreibt die Ausbreitung einer Konzentration in einem gegebenen Strömungsfeld (z.B. Farbe in einem Fluss, Milch im Kaffee, $\ldots$). Die involvierten Grössen sind

$$\begin{align}
u(x) & \quad\text{Massenkonzentration [kg/m$^2$]} & f(x) & \quad\text{Quelle [kg/(s m$^2$)]}\\
j(x) & \quad\text{Massenstrom [kg/(s m)]} & b(x) & \quad\text{Strömungsfeld [m/s]}.
\end{align}$$

Wir wenden wieder die Massenerhaltung in einem beliebigen Kontrollvolumen $B$ an:

$$\boxed{\begin{array}{l}\text{Zugeführte Konzentration}\\\text{in $B$ pro Zeiteinheit}\end{array}} = \boxed{\begin{array}{l}\text{Massenabfluss durch $\partial B$}\\\text{pro Zeiteinheit}\end{array}}$$

Formal bedeutet das

$$\int_B f(x) dx = \int_{\partial B} j(x)\cdot n ds.$$

Nehmen wir an, dass die beteiligten Funktionen glatt sind, so folgt mit dem Satz von Gauss

$$\int_B f(x) dx = \int_{B} \mathop{div} j(x) dx.$$

Da dies für alle Testvolumina $B$ gilt, erhalten wir punktweise die Gleichung

$$f(x) = \mathop{div} j(x)\quad \text{für}\ x\in\Omega.$$

Mit dem **Transportmodell** $j(x) = b(x)\, u(x)$ folgt die Gleichung

$$\mathop{div}(b(x) u(x)) = f(x)\quad \text{für}\ x\in\Omega.$$

Die Konzentration muss als Randbedingung am Einflussrand vorgegeben werden.

$$u(x) = u_g(x)\quad \text{für}\ x\in\Gamma_{\text{in}}.$$

Am Ausflussrand ergibt sich die Konzentration aus der Differentialgleichung.

```{prf:remark}
:label: my-rm-mod2

Für die numerische Berechnung von Lösungen dieses Types von partiellen Differentialgleichungen bieten sich Finite Volumen Methoden (FVM) an. Die Methoden haben den Vorteil, dass die Massenerhalten per Konstruktion exakt erfüllt wird, was bei gewöhnlichen Finite Element Methoden nicht per se gegeben ist. Neue diskontinuierliche FEM Methoden lösen dieses Problem. In diesem Kurs können wir auf die Thematik nicht eingehen.
```

```{seealso}
[Discontinuous Galerkin Discretizations for linear transport](https://docu.ngsolve.org/latest/i-tutorials/unit-3.3-scalardg/scalardg.html)
```


## Parabolische partielle Differentialgleichungen

Sei $\Omega \subset \mathbb{R}^d$ und $T > 0$. Der Prototyp einer parabolischen Differentialgleichung ist: Gesucht ist $u(x,t) : \Omega \times [0,T] \to \mathbb{R}$ so, dass

$$\begin{split}
\frac{\partial u}{\partial t}(x,t) - \Delta u(x,t) & = f(x, t)\quad \text{für}\ x\in\Omega,\,t\in(0,T)\\
u(x,t) & = u_d(x) \quad \text{für}\ x\in\partial\Omega,\,t\in(0,T)\\
u(x,0) & = u_0(x) \quad \text{für}\ x\in\partial\Omega.
\end{split}$$

Eine parabolische Differentialgleichung benötigt für jeden Zeitpunkt $t>0$ Randbedingungen. Für den Startzeitpunkt $t=0$ wird eine Anfangsbedingung $u(x,0)$ benötigt.

Für die Modellierung betrachten wir die **instationäre Wärmeleitung**. Die Energiebilanz für das Zeitintervall $(t,t+\tau)$ und das Kontrollvolumen $B\subset \Omega$ ist gegeben durch

$$\boxed{\begin{array}{l}\text{Energie in $B$}\\\text{zum Zeitpunkt $t+\tau$}\end{array}} - \boxed{\begin{array}{l}\text{Energie in $B$}\\\text{zum Zeitpunkt $t$}\end{array}} = \boxed{\begin{array}{l}\text{Erzeugte Energie}\\\text{in $B$ während $(t,t+\tau)$}\end{array}} - \boxed{\begin{array}{l}\text{Energieabfluss durch}\\\text{$\partial B$ während $(t,t+\tau)$}\end{array}}$$

Die thermische Energie $W$ sei proportional zur Temperatur $T$, daher

$$W = c\,T.$$

Die Energiebilanz lautet damit formal:

$$\int_B c\, T(x,t+\tau) dx - \int_B c\, T(x,t) dx = \int_t^{t+\tau} \int_B f(x,\tilde{t}) dx\, d\tilde{t} - \int_t^{t+\tau} \int_{\partial B} (-k\,\nabla T(x,\tilde{t})\cdot n\, ds\, d\tilde{t}$$

Da die Bilanz für alle Zeitintervalle $(t,t+\tau)$ und für alle Kontrollvolumina $B\subset \Omega$ gilt, gilt unter Annahme glatter Funktionen die Gleichheit punktweise. Wir erhalten die parabolische PDE

$$c\,\frac{\partial u}{\partial t}(x,t) - \mathop{div}(k\, \nabla T(x,t)) = f(x, t).$$

## Hyperbolische partielle Differentialgleichungen

Sei $\Omega \subset \mathbb{R}^d$ und $T > 0$. Der Prototyp einer hyperbolischen Differentialgleichung ist die Wellengleichung: Gesucht ist $u(x,t) : \Omega \times [0,T] \to \mathbb{R}$ so, dass

$$\begin{split}
\frac{\partial^2 u}{\partial t^2}(x,t) - \Delta u(x,t) & = f(x, t)\quad \text{für}\ x\in\Omega,\,t\in(0,T)\\
u(x,t) & = u_d(x) \quad \text{für}\ x\in\partial\Omega,\,t\in(0,T)\\
u(x,0) & = u_0(x) \quad \text{für}\ x\in\partial\Omega,\\
\frac{\partial u}{\partial t}(x,0) & = v_0(x) \quad \text{für}\ x\in\partial\Omega.
\end{split}$$

Eine hyperbolische PDE benötigt für jeden Zeitpunkt $t>0$ Randbedingungen. Zur Anfangszeit $t=0$ wird eine Anfangsbedingung für $u$ und $\dot{u}$ benötigt.

Betrachten wir die Modellierung wiederum an einem physikalischen Problem, in dem Fall der akustischen Wellenausbreitung (Schallwellen). Dazu sei

$$\begin{split}
v & \quad\text{Teilchengeschwindigkeit [m/s]}\\
p & \quad\text{Luftdruck [Pa]}\\
\rho & \quad\text{Dichte [kg/m$^3$]}
\end{split}$$

Wir stellen die Dichte als $\rho(x,t) = \rho_0 + \rho_s(x,t)$ mit einer stationären Dichte $\rho_0$ und einer kleinen Schwankung $\rho_s$ dar. Analog für den Druck $p(x,t) = p_0 + p_s(x,t)$.

* Für die Beschleunigung der Luftteilchen setzen wir den zum Druckgradienten proportionalen Ansatz:

$$\frac{\partial \rho v}{\partial t} = \nabla p.$$

* Wir nehmen an, dass Druckänderungen proportional zu Dichteänderungen sind

  $$p_s(x,t) = c^2\,\rho_s(x,t).$$

* Eine Divergenz der Geschwindigkeit führt zu einer Dichteänderung, dh.

  $$\frac{\partial \rho}{\partial t} = -\mathop{div}(\rho\,v).$$

  Eine positive Divergenz bedeutet Expansion, was namensgebend für die Divergenz ist (lateinisch divergere „auseinanderstreben“) und was in der Realität mit einer Abnahme der Dichte einher geht.

Vernachlässigt man kleine Schwankungen gegenüber dem grossen stationären Wert in der Dichte, so erhalten wir

$$\frac{\partial}{\partial t}(\rho_0 + \rho_s(x,t)) v \approx \frac{\partial}{\partial t}\rho_0 v = \nabla (p_0 + p_s(x,t)) = \nabla p_s(x,t)$$

und 

$$\frac{\partial p_s}{\partial t} = c^2 \frac{\partial\rho_s(x,t)}{\partial t} = c^2 \mathop{div}(\rho\,v) \approx c^2 \rho_0 \mathop{div}(v)$$

Daher folgt

$$\begin{split}
\frac{\partial v}{\partial t} & = \frac{1}{\rho_0} \nabla p_s(x,t)\\
\frac{\partial p_s}{\partial t} & = c^2 \rho_0 \mathop{div}(v)\end{split}$$

Bildet man die Divergenz der ersten Gleichung, und die Zeitableitung der zweiten Gleichung, so erhält man

$$\begin{split}
\mathop{div}\frac{\partial v}{\partial t} & = \frac{1}{\rho_0} \Delta p_s(x,t)\\
\frac{\partial^2 p_s}{\partial t^2} & = c^2 \rho_0 \frac{\partial}{\partial t}\mathop{div}(v)\end{split}$$

Beide Gleichungen liefern somit

$$\frac{\partial^2 p_s}{\partial t^2} = c^2 \Delta p_s(x,t).$$

die Wellengleichung für $p_2$, wobei für die Quelle $f=0$ gilt. Eine harte, Schall reflektierende Wand kann durch die Randbedingung $v\cdot n = 0$ modelliert werden.

``````{seealso}
Ein paar Beispiele mit unterschiedlichen Randbedingungen für die Wellengleichung

* ```{figure} ../movies/Wellengleichung1D.mp4
  ---
  align: left
  height: 250px
  name: fig-Wellengleichung1d
  ---
  1d Beispiel Vergleich Neumann / Dirichlet Randbedingungen
  ```
* ```{figure} ../movies/movieDirichletRB2D.mpg
  ---
  align: left
  height: 250px
  name: fig-Wellengleichung2dDirichlet
  ---
  2d Beispiel mit Dirichlet Randbedingungen
  ```
* ```{figure} ../movies/movieNeumannRB2D.mpg
  ---
  align: left
  height: 250px
  name: fig-Wellengleichung2dNeumann
  ---
  2d Beispiel mit Neumann Randbedingungen
  ```
``````
