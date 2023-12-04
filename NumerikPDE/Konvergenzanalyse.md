(chap:konvergenzanalyse)=
# Existenz und Konvergenz

## Existenz einer Lösung

Betrachten wir die Konvergenz der Methode der finiten Elemente in abstrakter Form für elliptische Probleme. Dazu definieren wir eine wichtige Eigenschaft der Bilinearform

```{prf:definition} koerziv oder elliptisch
:label: my-def-koerz

Sei $V$ ein Hilbertraum und $A$ eine stetige Bilinearform auf $V$. Dann heisst $A$ **koerziv**, falls

$$\exists C > 0: \quad |A(v,v)| \ge C \|v\|_V^2\quad \forall\ v\in V.$$ (eq:koerzivUngleichung)
```

Der folgende Satz von Lax-Milgram liefert uns die Existenz von Lösungen für unsere schwache Gleichung {eq}`eq:weakProblemPoisson2`. Der Satz von Lax-Milgram basiert auf dem Darstellungssatz von Riesz, welchen wir als Abschluss im Kapitel {numref}`chap:lineareFunktionale` angetroffen haben.


```{prf:theorem} Satz von Lax-Milgram
:label: my-thm-LaxMilgram

Sei $V$ ein Hilbertraum, $f: V\to \mathbb{K}$ eine **stetige** Linearform und $A: V\times V\to \mathbb{K}$ eine **stetige** Bilinearform. Zu dem sei $A$ auf einem Teilraum $V_0\subset V$ **koerziv**. Dann hat das Variationsproblem: gesucht $u\in V_g$ so, dass

$$A(u,v) = f(v)\quad \forall v\in V_0$$ (eq:schwacheLoesungVariationsproblem)

und $V_g$ eine lineare Mannigfaltigkeit, gegeben durch

$$V_g = u_g + V_0 := \{u_g + v_0\,|\, v_0\in V_0\}$$

eine **eindeutige** Lösung.
```

Wir betrachten das etwas allgemeinere Problem als in der Gleichung {eq}`eq:weakProblemPoisson2` zugrunde liegende: 

$$\begin{split}
-\Delta u & = f(x)\quad\text{für}\ x\in \Omega\\
u(x) & = u_g(x)\quad\text{für}\ x\in \Gamma_D\\
\frac{\partial u}{\partial n}(x) & = g(x)\quad\text{für}\ x\in \Gamma_N.\end{split}$$

Der Rand $\partial\Omega$ des Gebietes $\Omega$ ist gegeben durch $\partial\Omega = \Gamma_D \cup \Gamma_N$, wobei $\Gamma_D \not= \emptyset$. Wir lassen daher anstelle der homogenen Dirichlet Randwerte beliebige Dirichlet Randwerte auf einem nicht verschwindenden Teilrand $\Gamma_D$ zu. 
Am Dirichletrand fordert man die Randbedingung $u(x) = u_g(x)$ explizit, entsprechend erlaubt man nur Testfunktionen mit $v(x)=0$ auf $\Gamma_D$. Am Neumannrand wird die Randbedingung $\frac{\partial u}{\partial n}(x) = g(x)$ eingesetzt. Testen und partiell Integrieren liefert

$$\int_\Omega \nabla u\cdot\nabla v\, dx - \int_{\Gamma_D} \frac{\partial u}{\partial n} \underbrace{v}_{=0}\ ds - \int_{\Gamma_N} \underbrace{\frac{\partial u}{\partial n}}_{=g} v\ ds = \int_\Omega f(x) v(x) dx$$ (eq:VariationsFormulierungPoissongleichung)

Wir erhalten somit das schwache Problem: finde $u\in H^1(\Omega)$ mit $u(x) = u_g(x)$ für $x\in\partial\Gamma_D$ und

$$\int_\Omega \nabla u\cdot\nabla v\, dx = \int_\Omega f(x) v(x) dx + \int_{\Gamma_N} g(x) v(x)\ ds \quad \forall\ v\ \text{mit $v(x)=0$ auf $\Gamma_D$}.$$

Um den Satz von Lax-Milgram anwenden zu können, wählen wir

* $V = H^1(\Omega)$ mit der Norm $\|v\|_V = \left(\|v\|_{L_2}^2 + \|\nabla v\|_{L_2}^2\right)^{1/2}$
* $V_0 = \{v\in V\ |\ v=0\, \text{auf}\ \Gamma_D\}$
* $V_g = \{v\in V\ |\ v=u_g\, \text{auf}\ \Gamma_D\}$

  Für ein beliebiges $\tilde{u}_g\in V$ mit $\tilde{u}_g = u_g$ auf dem Rand $\Gamma_D$ ist $V_g = \tilde{u}_g + V_0$.
* Linearform $b: V\to\mathbb{R}$ gegeben durch

  $$b(v)= \int_\Omega f(x) v(x) dx + \int_{\Gamma_N} g(x) v(x)
ds$$
* Bilinearform $A: V\times V\to\mathbb{R}$ gegeben durch

  $$A(u,v) = \int_\Omega \nabla u\cdot\nabla v\, dx.$$

Mit Hilfe des Satz von Riesz erhalten wir die *Stetigkeit der Linearform*. Es gilt

$$\begin{split}|b(v)| & = \left|(f,v)_{L_2(\Omega)} + (g,v)_{L_2(\Gamma_N)}\right|\\
& \le \|f\|_{L_2(\Omega)} \|v\|_{L_2(\Omega)} + \|g\|_{L_2(\Gamma_N)} \|v\|_{L_2(\Gamma_N)}\\
& \le \|f\|_{L_2(\Omega)} \|v\|_V + \|g\|_{L_2(\Gamma_N)} c_{tr} \|v\|_V\\
& = \left(\|f\|_{L_2(\Omega)} + \|g\|_{L_2(\Gamma_N)} c_{tr}\right)\, \|v\|_V.
\end{split}$$

In der Abschätzung haben wir die **Cauchy-Schwarz Ungleichung** in $L_2(\Omega)$ und die **Trace Ungleichung**

$$\|v\|_{L_2(\Gamma_N)} \le c_{tr} \|v\|_V\quad \forall\,v\in V$$

des Sobolev Raum $H^1$ angewandt. Die *Stetigkeit der Bilinearform* ist trivial es gilt

$$A(u,v) = \int_\Omega \nabla u\cdot \nabla v\, dx = (\nabla u, \nabla v)_{L_2} \le \|\nabla u\|_{L_2}\ \|\nabla v\|_{L_2} \le \|u\|_V\,\|v\|_V.$$

Um die Koerzivität zeigen zu können, benötigen wir die **Friedrichsungleichung**

$$\|v\|_{L_2(\Omega)} \le c_F \|\nabla v\|_{L_2(\Omega)}\quad \forall\,v\in V_0.$$

Hier kommt die Voraussetzung $\Gamma_D \not= \emptyset$ ins Spiel. Mit der Friedrichsungleichung gilt

$$\|v\|_V^2 = \|v\|_{L_2}^2 + \|\nabla v\|_{L_2}^2 \le (c_F^2+1) \|\nabla v\|_{L_2}^2 = (c_F+1) A(v,v)\quad \forall\, v\in V_0,$$

also folgt für die Konstante der koerziv Bedingung {eq}`eq:koerzivUngleichung` $C=(c_F+1)^{-1}$. Damit haben wir alle Voraussetzungen des Satzes von Lax-Milgram erfüllt und es folgt eine *eindeutige* Lösung der Variationsformulierung der Poissongleichung {eq}`eq:VariationsFormulierungPoissongleichung` im Raum $H^1$.

## Fehlerabschätzung

Sei $V_h\subset V$ ein endlich-dimensionaler Teilraum gegeben zum Beispiel durch stückweise lineare stetige Funktionen. Wir setzen $V_{0,h} = V_0 \cap V_h$ und $V_{g,h} = V_g \cap V_h.$
Damit lautet das **diskrete Problem**: Gesucht ist $u_h\in V_{g,h}$ mit

$$A(u_h,v_h) = f(v_h)\quad \forall v_h\in V_{0,h}.$$

$V_h$ ist wieder ein Hilbertraum. Daher liefert der Satz von Lax-Milgram auch für das diskrerte Problem eine eindeutige Lösung. Da $V_{0,h} \subset V_0$, gilt das unendlich-dimensionale Problem auch für diskrete Testfunktionen. Wir haben daher

$$\begin{split}
A(u,v_h) & = f(v_h)\quad \forall v_h\in V_{0,h}\\
A(u_h,v_h) & = f(v_h)\quad \forall v_h\in V_{0,h}.\end{split}$$

und entsprechend

$$A(u-u_h, v_h) = 0\quad \forall v_h\in V_{0,h}.$$

Der Fehler $u-u_h$ ist daher $A(\cdot,\cdot)$-orthogonal zum diskreten Testraum $V_{0,h}$. Dies nennt man **Gallerkin-Orthogonalität**.

Mit Hilfe der Korezivität, Galerkin Orthogonalität und der Stetigkeit von $A$ folgt das **Cea's Lemma**

```{prf:theorem} Cea's Lemma
:label: my-thm-Cea

Sei $A$ stetig und koerziv. Dann gilt

$$\|u-u_h\|_V \le \frac{\mu_1}{\mu_2}\,\inf_{w_h\in V_{g,h}} \|u-w_h\|_V,$$

wobei die Konstanten $\mu_1$ aus der Koerzivitäts Bedingung, $\mu_2$ aus der Stetigkeit der Bilinearform $A$ ist und für die schwachen Lösungen des Variationsproblems {eq}`eq:schwacheLoesungVariationsproblem` mit $u\in V_0$ sowie $u_h\in V_{0,h}$ gilt.
```

Damit haben wir die Konvergenzaussage, dass die FEM-Lösung bis auf den Faktor $\frac{\mu_1}{\mu_2}$ die **best-mögliche Approximation** zur exakten Lösung im diskreten Raum ist.

Weitere Abschätzungen erhalten wir aus der Polynom-Interpolation mit Hilfe des Interpolationsoperators $I_h: V \to V_h$ 

$$\inf_{w_h\in V_{g,h}} \|u-w_h\|_V \le \|u-I_hu\|_V$$

und der Interpolationsfehlerabschätzung in Sobolevräumen

$$\|u- I_h u\|_{H^1} \le c\,h\,\|u\|_{H^2}.$$