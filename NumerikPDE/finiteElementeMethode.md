# Methode der finiten Elemente

(chap:konvergenzanalyse)=
## Konvergenzanalyse

### Existenz einer Lösung

Betrachten wir die Konvergenz der Methode der finiten Elemente in abstrakter Form für elliptische Probleme. Dazu definieren wir eine wichtige Eigenschaft der Bilinearform

```{admonition} Definition: koerziv oder elliptisch
Sei $V$ ein Hilbertraum und $A$ eine stetige Bilinearform auf $V$. Dann heisst $A$ **koerziv**, falls

$$\exists C > 0: \quad |A(v,v)| \ge C \|v\|_V^2\quad \forall\ v\in V.$$ (eq:koerzivUngleichung)
```

Der folgende Satz von Lax-Milgram liefert uns die Existenz von Lösungen für unsere schwache Gleichung {eq}`eq:weakProblemPoisson2`. Der Satz von Lax-Milgram basiert auf dem Darstellungssatz von Riesz, welchen wir als Abschluss im Kapitel {numref}`chap:lineareFunktionale` angetroffen haben.


```{admonition} Satz: Satz von Lax-Milgram
Sei $V$ ein Hilbertraum, $f: V\to \mathbb{K}$ eine stetige Linearform und $A: V\times V\to \mathbb{K}$ eine stetige Bilinearform. Zu dem sei $A$ auf einem Teilraum $V_0\subset V$ **koerziv**. Dann hat das Variationsproblem: gesucht $u\in V_g$ so, dass

$$A(u,v) = f(v)\quad \forall v\in V_0$$

und $V_g$ eine lineare Mannigfaltigkeit, gegeben durch

$$V_g = u_g + V_0 := \{u_g + v_0\,|\, v_0\in V_0\}$$

eine *eindeutige* Lösung.
```

Wir betrachten das etwas allgemeinere Problem als in der Gleichung {eq}`eq:weakProblemPoisson2` zugrunde liegende: 

$$\begin{split}
-\Delta u & = f(x)\quad\text{für}\ x\in \Omega\\
u(x) & = u_g(x)\quad\text{für}\ x\in \Gamma_D\\
\frac{\partial u}{\partial n}(x) & = g(x)\quad\text{für}\ x\in \Gamma_D.\end{split}$$

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

### Fehlerabschätzung

Sei $V_h\subset V$ ein endlich-dimensionaler Teilraum gegeben zum Beispiel durch stückweise lineare stetige Funktionen. Wir setzen $V_{0,h} = V_0 \cap V_h$ und $V_{g,h} = V_g \cap V_h.$
Damit lautet das **diskrete Problem**: Gesucht ist $u_h\in V_{g,h}$ mit

$$A(u_h,v_h) = f(v_h)\quad \forall v_h\in V_{0,h}.$$

$V_h$ ist wieder ein Hilbertraum. Daher liefert der Satz von Lax-Milgram auch für das diskrerte Problem eine eindeutige Lösung. Da $V_{0,h} \subset V_0$ gilt, das unendlich-dimensionale Problem auch für diskrete Testfunktionen. Wir haben daher

$$\begin{split}
A(u,v_h) & = f(v_h)\quad \forall v_h\in V_{0,h}\\
A(u_h,v_h) & = f(v_h)\quad \forall v_h\in V_{0,h}.\end{split}$$

und entsprechend

$$A(u-u_h, v_h) = 0\quad \forall v_h\in V_{0,h}.$$

Der Fehler $u-h_h$ ist daher $A(\cdot,\cdot)$-orthogonal zum diskreten Testraum $V_{0,h}$. Dies nennt man **Gallerkin-Orthogonalität**.

Mit Hilfe der Korezivität, Galerkin Orthogonalität und der Stetigkeit von $A$ folgt das **Cea's Lemma**

```{admonition} Satz: Cea's Lemma
Sei $A$ stetig und koerziv. Dann gilt

$$\|u-u_h\|_V \le \frac{\mu_1}{\mu_2}\,\inf_{w_h\in V_{g,h}} \|u-w_h\|_V,$$

wobei $\mu_1$ die Konstante aus der Koerzivitäts Bedingung und $\mu_2$ die Konstante aus der Stetigkeit der Bilinearform $A$ ist.
```

Damit haben wir die Konvergenzaussage, dass die FEM-Lösung bis auf den Faktor $\frac{\mu_1}{\mu_2}$ die **best-mögliche Approximation** zur exakten Lösung im diskreten Raum ist.

Weitere Abschätzungen erhalten wir aus der Polynom-Interpolation mit Hilfe des Interpolationsoperators $I_h: V \to V_h$ 

$$\inf_{w_h\in V_{g,h}} \|u-w_h\|_V \le \|u-I_hu\|_V$$

und der Interpolationsfehlerabschätzung in Sobolevräumen

$$\|u- I_h u\|_{H^1} \le c\,h\,\|u\|_{H^2}.$$

## Assembling

### Eindimensionaler Fall

Im eindimensionalen Fall ist eine Zerlegung des Rechengebietes gegeben durch die Teilintervalle $T_i = [x_i,x_{i+1}]$, wobei die Punkte

$$0=x_0 < x_1 < \ldots < x_{n-1} < x_n = 1$$

gegeben sind. Es gilt daher 

$$[0,1] = \bigcup_{i=1}^n T_i.$$

Wir betrachten nun die Matrix $A$ {eq}`eq:linBilinForm1dproblem`. Die praktische Berechnung führt nicht über alle Kombinationen $j, k = 1, \ldots, n$. Die Basis Funktionen sind abhängig von der Zerlegung und können daher nicht vorab definiert werden. Ziel des Assembling ist, die Integration einmal mit Hilfe eines Referenz-Intervalls / Referenz-Elements durchzuführen und auf diese Berechnung zurück zu greifen.

Dazu betrachtet man das Integral in jedem Teilintervall und summiert diese auf. Entsprechend bezeichnet man $A$ als die globale Matrix und $\tilde{A}_i$ als lokale Matrix: 

$$\begin{split}
A & = (a_{j k}) = \underbrace{\int_0^1\varphi_j'(x)\cdot \varphi_k'(x) dx}_{\text{global}} = \sum_{i=1}^n \underbrace{\int_{T_i} \varphi_j'(x)\cdot \varphi_k'(x) dx}_{\text{lokal (pro Element)}}.
\end{split}$$

Mit der Transformation

$$\begin{array}{rl}\tau_i : T_i & \to [0,1]\\
x & \mapsto \tau_i(x) = \frac{x-x_i}{x_{i+1}-x_i} \quad\end{array}\quad\text{bzw.}\quad 
\begin{array}{rl}\sigma_i : [0,1] & \to  T_i\\
t & \mapsto \sigma_i(t) = x_i + (x_{i+1}-x_i)\cdot t\end{array}$$

folgt für das Integral

$$\begin{split}
\int_{T_i} \varphi_j'(x)\cdot \varphi_k'(x)\, dx & = \int_{x_i}^{x_{i+1}} \varphi_j'(x)\cdot \varphi_k'(x)\, dx\\
&  = \int_{\tau_i(x_i)}^{\tau_i(x_{i+1})} \varphi_j'(\sigma_i(t))\cdot \varphi_k'(\sigma_i(t))\cdot \underbrace{\dot{\sigma}_i(t)\cdot dt}_{= dx}.
\end{split}$$ (eq:femlokalglobal)

Sei nun $\tilde{\varphi}_k(t)$ die auf dem Einheitsintervall $[0,1]$ definierte Basisfunktion, dann gilt

$$\tilde{\varphi}(t) = \varphi(\sigma_i(t)).$$

Für die Ableitung folgt somit

$$\frac{d}{dt}\tilde{\varphi}(t) = \varphi'(\sigma_i(t))\cdot \dot{\sigma}_i(t)\quad\Rightarrow\quad 
\varphi'(\sigma_i(t)) = \frac{\dot{\tilde{\varphi}}(t)}{\dot{\sigma}_i(t)}.$$

Schreibt man das Integral {eq}`eq:femlokalglobal` mit den Basisfunktionen auf dem Einheitsintervall, so folgt

$$\begin{split}
\int_{T_i} \varphi_j'(x)\cdot \varphi_k'(x)\, dx & = \int_{\tau_i(x_i)}^{\tau_i(x_{i+1})} \dot{\tilde{\varphi}}_j(t) \cdot \dot{\tilde{\varphi}}_k(t) \frac{\dot{\sigma}_i(t)}{\dot{\sigma}_i(t)^2} dt = \int_{0}^{1} \dot{\tilde{\varphi}}_j(t) \cdot \dot{\tilde{\varphi}}_k(t) \frac{1}{\dot{\sigma}_i(t)} dt
\end{split}$$

Da die Transformation eine affine Funktion ist ($\tau, \sigma$ sind Polynome ersten Grades), gilt $\dot{\sigma}\equiv \text{const}$. Für die globale Matrix $A$ folgt somit

$$A = \sum_{i=1}^n \frac{1}{x_{i+1}-x_i} \int_{0}^{1}\tilde{\varphi}_j'(t)\cdot \tilde{\varphi}_k'(t) dt = \sum_{i=1}^n \frac{1}{h_i} \underbrace{\int_{0}^{1}\tilde{\varphi}_j'(t)\cdot \tilde{\varphi}_k'(t) dt}_{\text{Steifigkeitselementmatrix}}\quad\text{mit}\ h_i = x_{i+1}-x_i,$$

wobei das Integral stets über das Einheitsintervall geht. Wir können daher den finite Element Ansatz auf dem Einheitsintervall studieren und die globale Matrix durch Addition an den entsprechenden Matrixpositionen zusammenstellen (man spricht hier vom `Assembling`).
Das Produkt $\tilde{\varphi}_j(t)\cdot \tilde{\varphi}_k(t)$ ist für fast alle $j, k$ Kombinationen null. In der Abbildung {numref}`fig-fembasissupport` ist der Fall für Polynome mit Grad 1 dargestellt. In der Abbildung {numref}`fig-FEMlocalBasisProductp13` sind die Kombinationen der *lokalen* Basisfunktionen für unterschiedliche Polynomgrade dargestellt.

```{figure} fembasissupport.png
---
align: left
height: 400px
name: fig-fembasissupport
---

Produkt $\varphi_j \cdot \varphi_k$ auf dem Intervall $[0,1]$ für stückweise Polynome mit Grad 1.
```

```{figure} FEMlocalBasisProductp13.png
---
align: left
height: 400px
name: fig-FEMlocalBasisProductp13
---

Lokale Kombinationen $\varphi_j \cdot \varphi_k$ der nodalen FEM Basisfunktionen auf dem Einheitsintervall $[0,1]$.
```

Die **globalen Matrizen** können daher mit lokalen Elementmatrizen wie folgt berechnet werden

$$\begin{split}
M & = \sum_{i=1}^n h_i \int_0^1 \tilde{\varphi}_j(t)\cdot \tilde{\varphi}_k(t)dt\\
A & = \sum_{i=1}^n \frac{1}{h_i} \int_0^1 \dot{\tilde{\varphi}}_j(t)\cdot \dot{\tilde{\varphi}}_k(t)dt
\end{split}$$

Man nennt die Matrix

$$S_e = \int_0^1 \tilde{\varphi}_j'(x)\cdot \tilde{\varphi}_j'(x) dx$$

**Steiffigkeits-Elementmatrix**

und 

$$M_e = \int_0^1 \tilde{\varphi}_j(x)\cdot \tilde{\varphi}_j(x) dx$$

**Massen-Elementmatrix**.

Die Dimension ist vom Polynomgrad $p$ abhängig. In den Tabellen {numref}`chap:1dFEMElemente` sind die $S_e, M_e$ für die Polynomgrade $p=1,2,3$ berechnet, analog zur Abbildung {numref}`fig-FEMlocalBasisProductp13`.

```{admonition} Aufgabe
Berechne die Elementmatrizen basierend auf den Lagrange Polynome bezogen auf das Einheitsintervall $[0,1]$ für

* Steifigkeitsmatrix

  $$A_{i,j} = \int_0^1 \varphi_i'(x) \varphi_j'(x)\,dx$$

* Massenmatrix

  $$M_{i,1} = \int_0^1 \varphi_i(x) \varphi_j(x)\,dx$$

für die Ordnungen (Polynom Grad der Basisfunktionen) $p = 1, \ldots, 3$.
```

```{admonition} Aufgabe
Berechne die (globale) System Matrix und Vektor das Problem {eq}`eq:schwacheGleichungBeispiel` mit Hilfe der lokalen Matrizen. 
```

### Zweidimensionaler Fall

Eine reguläre Triangulierung $\mathcal{T} = \{T_1, \ldots, T_M\}$ eines Gebiets $\Omega$ ist die Zerlegung in Dreiecke $T_i$ so, dass $\bar{\Omega} = \cup_i T_i$ und $T_i\cap T_j$ ist
* entweder leer
* oder hat eine gemeinsame Kante.

In einem erweiterten Sinne kann die Triangulierung aus verschiedenen Elemente bestehen: Dreiecke, Vierecke, (Tetraeder, Hexeder, Prismen, Pyramiden im dreidimensionalen). Die finite Elemente werden typischerweise, wie wir es im eindimensionalen gemacht haben, auf einem Referenzelement definiert. Die einzelnen Elemente der Zerlegung können mit Hilfe einer affinen Transformation und dem Referenzelement beschrieben werden. 

Die Koordinaten Transformation im zweidimensionalen erfordert etwas mehr Rechenaufwand.

```{figure} KoordinatenTransformation2d.png
---
align: left
height: 250px
name: fig-KoordinatenTransformation2d
---

Transformation auf Einheitsdreieck in $2D$
```

Ein Dreieck $T_i$ in allgemeiner Lage mit den Eckpunkten $P_1(x_1,y_1)$, $P_2(x_2,y_2)$ und  $P_3(x_3,y_3)$, welche im Gegenuhrzeigersinn fortlaufend numeriert seien, wie dies in Abb. {numref}`fig-KoordinatenTransformation2d` erfolgte, kann mittels der linearen Transformation

$$\begin{split}
x & = x_1 + (x_2-x_1)\, \xi + (x_3-x_1)\, \eta\\
y & = y_1 + (y_2-y_1)\, \xi + (y_3-y_1)\, \eta
\end{split}$$

bijektiv auf das Einheitsdreieck $T$ abgebildet werden. Das Flächenelement $dx dy$ kann mit der Jacobi-Determinante

$$J = (x_2-x_1)(y_3-y_1) - (x_3-x_1)(y_2-y_1)$$

der Transformation durch

$$dx\, dy = J\, d\xi\, d\eta$$

ersetzt werden. Die partiellen Ableitungen transformieren sich nach der Kettenregel gemäss

$$\begin{split}
u_x & = u_\xi \xi_x + u_\eta \eta_x\\
u_y & = u_\xi \xi_y + u_\eta \eta_y.
\end{split}$$

Die Gebietesintegrale transformieren sich somit

$$\int_{T_i} \nabla \varphi_j \cdot \nabla \varphi_k dx dy = a \underbrace{\int_T (\partial_\xi\varphi_j) (\partial_\xi\varphi_k) dx}_{=:S_1} + b \underbrace{2 \int_T (\partial_\xi\varphi_j) (\partial_\eta\varphi_k) dx}_{S_2} + c \underbrace{\int_T (\partial_\eta\varphi_j) (\partial_\eta\varphi_k) dx}_{S_3}$$

mit

$$\begin{split}
a & = \frac{1}{J} ((x_3-x_1)^2+(y_3-y_1)^2)\\
b & = -\frac{1}{J} ((x_3-x_1)(x_2-x_1)+(y_3-y_1)(y_2-y_1))\\
c & = \frac{1}{J} ((x_2-x_1)^2 + (y_2-y_1)^2)\\
J & = (x_2-x_1)(y_3-y_1) - (x_3-x_1)(y_2-y_1)
\end{split}$$

Für die Elementmatrizen mit lokalen Basisfunktionen gegeben durch Polynomen 1. Ordnung ergibt sich

$$\begin{split}
S_e & = a S_1 + b S_2 + c S_3\\
M_e & = J S_4
\end{split}$$

mit

$$\begin{split}
S_1 & = \frac{1}{2}\begin{pmatrix}
 1 & -1 & 0\\
 -1 & 1 & 0\\
 0 & 0 & 0	
 \end{pmatrix}
\quad
S_2 = \frac{1}{2}\begin{pmatrix}
 2 & -1 & -1\\
 -1 & 0 & 1\\
 -1 & 1 & 0	
 \end{pmatrix}\\
 S_3 & = \frac{1}{2}\begin{pmatrix}
 1 & 0 & -1\\
 0 & 0 & 0\\
 -1 & 0 & 1	
 \end{pmatrix}
 \quad
 S_4 = \frac{1}{24}\begin{pmatrix}
 2 & 1 & 1\\
 1 & 2 & 1\\
 1 & 1 & 2	
 \end{pmatrix}.
\end{split}$$

> Es erweist sich mathematisch wie auch Software technisch als bedeutend effizienter, die Berechnung der System Matrizen über die Triangulierung zu berechnen und die einzelnen Beiträge in der globalen Matrix aufzukummulieren. Diesen Prozess nennt man **Assembling**.
>
> Allgemein können wir dies in der Form
>
> $$A = \sum_{T \in \mathcal{T}} C_T A_T C_T^T$$
>
> und 
>
> $$f = \sum_{T \in \mathcal{T}} C_T f_T$$
>
> schreiben, wobei $C_T$ die Verknüpfung zwischen lokalen und globalen Funktionen darstellt.
