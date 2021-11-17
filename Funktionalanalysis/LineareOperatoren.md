# Lineare Operatoren

Viele Aufgaben in der Mathematik und den Anwendungen führen auf Gleichungen der Form

$$T(x) = T x = y,$$

wobei $T: V \to W$ eine "lineare Abbildung", $V, W$ normierte Räume sind und $y\in W$ ein gegeben ist.

```{prf:definition} lineare Abbildung
:label: my-def-linabb

Die _Abbildung_ $T$ des normierten Raumes $V$ in den normierten Raum $W$ heisst _linear_, wenn für alle $x,y \in V$ und alle $\alpha \in \mathbb{K}$

$$T(x+y) = T x + T y$$
```

## Beschränkte lineare Operatoren

Analog zum Stetigkeitsbegriff aus der Analysis definiert man die Stetigkeit und Beschränktheit bei linearen Operatoren wie folgt:

```{prf:definition} stetig
:label: my-def-stetig

Sei $T: V\to W$ ein linearer Operator und $V, W$ normierte Räume. Der lineare Operator $T$ heisst *stetig* in $x_0\in V$, wenn es zu jedem $\varepsilon > 0$ ein $\delta = \delta(\epsilon, x_0)>0$ gibt, so dass

$$\| T x- T x_0\|_W < \varepsilon\quad\text{für alle}\ x\in V \quad\text{mit}\ \|x-x_0\|_V < \delta$$

gilt. Man nennt $T$ *stetig in* $V$, wenn $T$ in jedem Punkt von $V$ stetig ist.
```

```{prf:definition} beschränkt
:label: my-def-beschraenkt

Es seien $V, W$ normierte Räume. Der lineare Operator $T: V \to W$ heisst *beschränkt*, wenn es eine Konstante $C>0$ mit

$$\|T x\|_W \le C\, \|x\|_V\quad\forall\,x\in V$$ (eq:beschrlinop)
gibt.
```

Die Menge aller linearen Abbildungen mit der sogenannten Operatornorm versehen, liefert uns wieder einen normierten Vektorraum. Die Operatornorm (vgl. Matrixnorm aus der Numerik oder linearen Algebra) ist gegeben durch

```{prf:definition} Operatornorm
:label: my-def-opnorm

Die kleinste Zahl $C>0$ für die {eq}`eq:beschrlinop` gilt, heisst *Operatornorm von* $T$ und ist durch

$$\|T\|_{V\to W} := \sup_{\substack{x\in V\\x\not= 0}} \frac{\|T x\|_W}{\|x\|_V}$$ (eq:operatornorm)

definiert.
```

```{prf:remark}
Mit der Norm von $T$ lässt sich die Ungleichung {eq}`eq:beschrlinop` auch in der Form

$$\|T x\|_W \le \|T\|_{V\to W}\ \|x\|_V$$ (eq:operatornormungleichung)

schreiben. Die Operatornorm lässt sich anstelle der Schreibweise {eq}`eq:operatornorm` auch durch

$$\|T\|_{V\to W} = \sup_{\substack{x\in V\\x\not= 0}} \frac{\|T x\|_W}{\|x\|_V} = \sup_{\|x\|_V \le 1} \|T x\|_W$$

schreiben.
```

:::{seealso}
[Beispiel zur Operatornorm.](BeispielLineareOperatoren.ipynb)
:::

Zwischen stetigen und beschränkten Operatoren besteht der Zusammenhang

```{prf:theorem}
:label: my-thm-beschrstetig

Sei $T: V\to W$ ein linearer Operator und $V,W$ normierte Räume. Dann gilt

$$T\ \text{ist beschränkt} \quad\Leftrightarrow \quad T\ \text{ist stetig.}$$
```

````{prf:proof}
a) Sei $T$ beschränkt durch $C>0$ und $x_0\in V$ beliebig. Wir zeigen nun, dass es zu jedem $\varepsilon>0$ ein $\delta=\delta(x_0)$ gibt, so dass

$$\| T x - T x_0\|_W < \epsilon\quad \forall \|x -x_0\|_V < \delta.$$
Wähle $\delta = \frac{\epsilon}{C}$, dann folgt

$$\| T x - T x_0\|_W = \|T (x-x_0)\|_W \le C \underbrace{\|x-x_0\|_V}_{< \delta} < C \frac{\varepsilon}{C} = \varepsilon$$
dh. $T$ ist in $x_0$ stetig und da $x_0\in V$ beliebig ist, in ganz $V$.

b) Sei nun $T$ auf $V$ stetig. Wir zeigen im Widerspruch, dass $T$ beschränkt ist. Daher nehmen wir an: $T$ sei nicht beschränkt. Also gibt es eine Folge $\{x_k\} \subset V$ mit $x_k \not=0$ und $\frac{\|T x_k\|_W}{\|x_k\|_V} > k$ für alle $k\in \mathbb{N}$. Setze nun $y_k = \frac{x_k}{k \|x_k\|_V}$, so folgt $y_k\in V$ und

$$\|T y_k\|_W = \left\| T\left(\frac{x_k}{k \|x_k\|_V}\right)\right\|_W = \frac{\|T x_k\|_W}{k \|x_k\|_V} > 1$$ (eq:bewstetigbeschraenkt)

für alle $k\in \mathbb{N}$. Andererseits gilt: $\|y_k\|_V = \frac{1}{k} \to 0$ für $k\to \infty$ bzw. $y_k \to 0$. Aus der Stetigkeit von $T$ in $0$ folgt $T y_k \to 0$ für $k\to\infty$, was im Widerspruch zu {eq}`eq:bewstetigbeschraenkt` steht.
````

```{prf:remark}
Bei linearen Operatoren sind Stetigkeit und Beschränktheit äquivalente Eigenschaften.
```

```{prf:example}
:label: my-example

Sei $V=W=C[a,b]$, $f\in C[a,b]$ und $\|f\| = \max_{a\le x \le b} |f(x)|$. Betrachte den *Integraloperator* $T$ mit

$$(T f)(x) := \int_a^b K(x,y) f(y) dy,\quad x\in [a,b]$$

mit $K: [a,b] \times [a,b] \to \mathbb{R}$ stetigem *Kern*.

* $T$ ist ein linearer Operator, der $C[a,b]$ auf sich abbildet.
* Da $f$ stetig auf $[a,b]$ ist, ist $f\cdot K$ stetig auf $[a,b]\times [a,b]$ und damit ist auch (vgl. Satz 7.17 in {cite:p}`burg_wille_haf_meister_AnalysisI`)

  $$F(x) := \int_a^b K(x,y) f(y) dy$$

  stetig auf $[a,b]$.
* Da $K$ stetig auf $[a,b] \times [a,b]$ ist, existiert

  $$M = \max_{x,y \in [a,b]} |K(x,y)|$$

  und somit

  $$\|T f|| = \max_{x\in[a,b]} |(T f)(x)| \le \max_{x\in[a,b]} \int_a^b |K(x,y)| |f(y)| dy \le M\,\max_{x\in[a,b]} |f(x)|\ \int_a^b 1 dy = M\,(b-a)\,\|f\|.$$

  Damit ist $T$ bezüglich der Maximumsnorm beschränkt. Es gilt

  $$\|T\| = \sup_{\|f\|=1} \|T f\| \le M\, (b-a).$$
```

```{prf:theorem}
:label: my-thm-opnorm

Seien $V, W$ normierte Räume und mit $L(V,W)$ die Menge aller beschränkten linearen Operatoren von $V$ in $W$ bezeichnet. Dann ist $L(V,W)$ bezüglich der Operatornorm

$$\|T\| = \sup_{\substack{x\in V\\x\not= 0}} \frac{\|T x\|_W}{\|x\|_V}$$

wieder ein normierter Raum.

Ist $V$ ein normierter Raum und $W$ ein *Banachraum*, dann ist $L(V,W)$ ein Banachraum.
```

(chap:lineareFunktionale)=
## Lineare Funktionale

Wir betrachten nun spezielle lineare Operatoren, welche in den zugrunde liegenden Zahlenkörper abbilden.

```{prf:definition} lineares Funktional
:label: my-def-linfunktional

Sei $V$ ein normierter Raum. Dann nennt man den linearen Operator $F: V \to \mathbb{K}$ ($\mathbb{R}$ oder $\mathbb{C}$) *lineares Funktional*.
```

Wir werden den Begriff **Linearform** für lineare Funktionale und die **Bilinearform** im Folgenden häufig antreffen. 

```{prf:definition} Linearform, Bilinearform
:label: my-def-linbilinform

* Linearform $f(\cdot)$ ist eine Abbildung

  $$\begin{split} f: V & \to \mathbb{K}\\
  v & \mapsto f(v)\quad \text{linear}\end{split}$$  (eq:Linearform)

* Bilinearform $A(\cdot, \cdot)$ ist eine Abbildung

  $$\begin{split} A: V \times V & \to \mathbb{K}\\
  (u,v) & \mapsto A(u,v)\quad \text{linear in $u$ und $v$.}\end{split}$$ (eq:Bilinearform)
```

```{prf:remark}
* Die Linearform ist ein äquivalenter Begriff für ein lineares Funktional.

* Der Stetigkeitsbegriff linearer Operatoren auf das lineare Funktional bzw. Linearform und auf die Bilinearform angewandt bedeutet:

  * die Linearform $f$ heisst *stetig*, falls

    $$\exists\ C>0:\quad |f(v)| \le C\ \|v\|_V\quad \forall\ v\in V$$

  * die Bilinearform $A(u,v)$ heisst *stetig*, falls

    $$\exists\ C>0:\quad |A(u,v)| \le C\ \|u\|_V\,\|v\|_V\quad \forall\ u,v\in V.$$
```

Normiert man $\mathbb{K}$ mit

$$\|\alpha\| := |\alpha|,\quad\alpha\in\mathbb{K}$$

so ist $\mathbb{K}$ ein Banachraum und damit die Menge aller beschränkten linearen Funktionale auf $V$ ein Banachraum. Dieser Raum ist insbesondere im Zusammenhang mit partiellen Differentialgleichungen sehr wichtig.

```{admonition} Dualraum von $V$
Der Banachraumm $L(V,\mathbb{K})$ aller beschränkten linearen Funktionale auf $V$ heisst der zu $V$ *konjugierte* oder *duale Raum* und wird mit $V^*$ oder $V'$ bezeichnet.
```

**Beispiel**: Sei $V$ ein Hilbertraum und $y_0$ ein beliebiges (festes) Element aus $V$. Für $x\in V$ wird durch

$$\begin{split}
F : V & \to \mathbb{R}\\
 x & \mapsto y = F x := (x,y_0)\end{split}$$

ein lineares Funktional $F$ definiert. (Die Linearität folgt direkt aus den Eigenschaften des Skalarprodukts.) Mit Hilfe der Schwarzschen Ungleichung folgt

$$\|F x\| = |F x| = |(x,y_0)| \le \|x\|\,\|y_0\|\quad\forall\ x\in V.$$

Damit folgt

$$\frac{\|F x\|}{\|x\|} \le \|y_0\|\quad \forall x\in V,\ \text{mit}\,x\not= 0$$

sprich $F$ ist ein beschränktes lineares Funktional, $F\in V^*$ mit $\|F\| \le \|y_0\|$. Da für $x=y_0$

$$\|F y_0\| = |(y_0,y_0)| = \|y_0\|^2 = \|y_0\|\, \|y_0\|$$

gilt, folgt

$$\|F\| = \sup_{\substack{x\in V\\x\not= 0}} \frac{\|F x\|}{\|x\|} = \|y_0\|.$$

Wir kommen nun zum Rieszschen Darstellungssatz: Beschränkte lineare Funktionale eines Hilbertraumes $V$ lassen sich besonders einfach darstellen. Die Darstellung aus obigem Beispiel erfasst **alle** beschränkten linearen Funktionale. Es gilt

```{prf:theorem} Darstellungssatz von Riesz
:label: my-thm-Riesz

Sei $V$ ein Hilbertraum und $F\in V^*$ beliebig. Dann gibt es ein *eindeutig* bestimmtes $y\in V$, so dass $F$ die Darstellung

$$F x = (x,y)\quad\forall x\in V$$

besitzt.
```

```{prf:remark}
Dieser Satz ist das zentrale Ergebnis der Hilbertraum-Theorie. Neben seiner Bedeutung als Darstellungssatz kann er auch als Existenz- und Eindeutigkeitsprinzip aufgefasst werden. Diese Bedeutung des Rieszschen Satzes ist Grundlage für die moderne Theorie der elliptischen partiellen Differentialgleichungen (vgl. auch {numref}`chap:konvergenzanalyse`).
```
