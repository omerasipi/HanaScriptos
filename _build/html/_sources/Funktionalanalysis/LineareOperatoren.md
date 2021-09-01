# Lineare Operatoren

Viele Aufgaben in der Mathematik und den Anwendungen führen auf Gleichungen der Form

$$T(x) = T x = y,$$

wobei $T: X \to Y$ eine "lineare Abbildung", $X, Y$ normierte Räume sind und $y\in Y$ ein gegeben ist.

```{admonition} Definition: lineare Abbildung
Die _Abbildung_ $T$ des normierten Raumes $X$ in den normierten Raum $Y$ heisst _linear_, wenn für alle $x,y \in X$ und alle $\alpha \in \mathbb{K}$

$$T(x+y) = T x + T y$$
```

## Beschränkte lineare Operatoren

Analog zum Stetigkeitsbegriff aus der Analysis definiert man die Stetigkeit und Beschränktheit bei linearen Operatoren wie folgt:

```{admonition} Definition: stetig
Sei $T: X\to Y$ ein linearer Operator und $X, Y$ normierte Räume. Der lineare Operator $T$ heisst *stetig* in $x_0\in X$, wenn es zu jedem $\varepsilon > 0$ ein $\delta = \delta(\epsilon, x_0)>0$ gibt, so dass

$$\| T x- T x_0\|_Y < \varepsilon\quad\text{für alle}\ x\in X \quad\text{mit}\ \|x-x_0\|_X < \delta$$

gilt. Man nennt $T$ *stetig in* $X$, wenn $T$ in jedem Punkt von $X$ stetig ist.
```

```{admonition} Definition: beschränkt
Es seien $X, Y$ normierte Räume. Der lineare Operator $T: X \to Y$ heisst *beschränkt*, wenn es eine Konstante $C>0$ mit

$$\|T x\|_Y \le C\, \|x\|_X\quad\forall\,x\in X$$ (eq:beschrlinop)
gibt.
```

Die Menge aller linearen Abbildungen mit der sogenannten Operatornorm versehen, liefert uns wieder einen normierten Vektorraum. Die Operatornorm (vgl. Matrixnorm aus der Numerik oder linearen Algebra) ist gegeben durch

```{admonition} Definition: Operatornorm
Die kleinste Zahl $C>0$ für die {eq}`eq:beschrlinop` gilt, heisst *Operatornorm von* $T$ und ist durch

$$\|T\|_{X\to Y} := \sup_{\substack{x\in X\\x\not= 0}} \frac{\|T x\|_Y}{\|x\|_X}$$ (eq:operatornorm)

definiert.
```
**Bemerkung**: Mit der Norm von $T$ lässt sich die Ungleichung {eq}`eq:beschrlinop` auch in der Form

$$\|T x\|_Y \le \|T\|_{X\to Y}\ \|x\|_X$$

schreiben. Die Operatornorm lässt sich anstelle der Schreibweise {eq}`eq:operatornorm` auch durch

$$\|T\|_{X\to Y} = \sup_{\substack{x\in X\\x\not= 0}} \frac{\|T x\|_Y}{\|x\|_X} = \sup_{\|x\|_X \le 1} \|T x\|_Y$$

schreiben.

Zwischen stetigen und beschränkten Operatoren besteht der Zusammenhang

```{admonition} Satz
Sei $T: X\to Y$ ein linearer Operator und $X,Y$ normierte Räume. Dann gilt

$$T\ \text{ist beschränkt} \quad\Leftrightarrow \quad T\ \text{ist stetig.}$$
```

**Beweis**:
a) Sei $T$ beschränkt durch $C>0$ und $x_0\in X$ beliebig. Wir zeigen nun, dass es zu jedem $\varepsilon>0$ ein $\delta=\delta(x_0)$ gibt, so dass

$$\| T x - T x_0\|_Y < \epsilon\quad \forall \|x -x_0\|_X < \delta.$$
Wähle $\delta = \frac{\epsilon}{C}$, dann folgt

$$\| T x - T x_0\|_Y = \|T (x-x_0)\|_Y \le C \underbrace{\|x-x_0\|_X}_{< \delta} < C \frac{\varepsilon}{C} = \varepsilon$$
dh. $T$ ist in $x_0$ stetig und da $x_0\in X$ beliebig ist, in ganz $X$.

b) Sei nun $T$ auf $X$ stetig. Wir zeigen im Widerspruch, dass $T$ beschränkt ist. Daher nehmen wir an: $T$ sei nicht beschränkt. Also gibt es eine Folge $\{x_k\} \subset X$ mit $x_k \not=0$ und $\frac{\|T x_k\|_Y}{\|x_k\|_X} > k$ für alle $k\in \mathbb{N}$. Setze nun $y_k = \frac{x_k}{k \|x_k\|_X}$, so folgt $y_k\in X$ und

$$\|T y_k\|_Y = \left\| T\left(\frac{x_k}{k \|x_k\|_X}\right)\right\|_Y = \frac{\|T x_k\|_Y}{k \|x_k\|_X} > 1$$ (eq:bewstetigbeschraenkt)

für alle $k\in \mathbb{N}$. Andererseits gilt: $\|y_k\|_X = \frac{1}{k} \to 0$ für $k\to \infty$ bzw. $y_k \to 0$. Aus der Stetigkeit von $T$ in $0$ folgt $T y_k \to 0$ für $k\to\infty$, was im Widerspruch zu {eq}`eq:bewstetigbeschraenkt` steht. $\Box$

> Bei linearen Operatoren sind Stetigkeit und Beschränktheit äquivalente Eigenschaften.

**Beispiel**: Sei $X=Y=C[a,b]$, $f\in C[a,b]$ und $\|f\| = \max_{a\le x \le b} |f(x)|$. Betrachte den *Integraloperator* $T$ mit

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

```{admonition} Satz
Seien $X, Y$ normierte Räume und mit $L(X,Y)$ die Menge aller beschränkten linearen Operatoren von $X$ in $Y$ bezeichnet. Dann ist $L(X,Y)$ bezüglich der Operatornorm

$$\|T\| = \sup_{\substack{x\in X\\x\not= 0}} \frac{\|T x\|_Y}{\|x\|_X}$$

wieder ein normierter Raum.

Ist $X$ ein normierter Raum und $Y$ ein *Banachraum*, dann ist $L(X,Y)$ ein Banachraum.
```

## Lineare Funktionale

Wir betrachten nun spezielle lineare Operatoren, welche in den zugrunde liegenden Zahlenkörper abbilden.

```{admonition} Definition: lineares Funktional
Sei $X$ ein normierter Raum. Dann nennt man den Operator $F: X \to \mathbb{K}$ ($\mathbb{R}$ oder $\mathbb{C}$) *lineares Funktional*.
```

Normiert man $\mathbb{K}$ mit

$$\|\alpha\| := |\alpha|,\quad\alpha\in\mathbb{K}$$

so ist $\mathbb{K}$ ein Banachraum und damit die Menge aller beschränkten linearen Funktionale auf $X$ ein Banachraum. Dieser Raum ist insbesondere im Zusammenhang mit partiellen Differentialgleichungen sehr wichtig.

```{admonition} Dualraum von $X$
Der Banachraumm $L(X,\mathbb{K})$ aller beschränkten linearen Funktionale auf $X$ heisst der zu $X$ *konjugierte* oder *duale Raum* und wird mit $X^*$ oder $X'$ bezeichnet.
```

**Beispiel**: Sei $X$ ein Hilbertraum und $y_0$ ein beliebiges (festes) Element aus $X$. Für $x\in X$ wird durch

$$\begin{split}
F : X & \to \mathbb{R}\\
 x & \mapsto y = F x := (x,y_0)\end{split}$$

ein lineares Funktional $F$ definiert. (Die Linearität folgt direkt aus den Eigenschaften des Skalarprodukts.) Mit Hilfe der Schwarzschen Ungleichung folgt

$$\|F x\| = |F x| = |(x,y_0)| \le \|x\|\,\|y_0\|\quad\forall\ x\in X.$$

Damit folgt

$$\frac{\|F x\|}{\|x\|} \le \|y_0\|\quad \forall x\in X,$$

sprich $F$ ist ein beschränktes lineares Funktional, $F\in X^*$ mit $\|F\| \le \|y_0\|$. Da für $x=y_0$

$$\|F y_0\| = |(y_0,y_0)| = \|y_0\|^2 = \|y_0\|\, \|y_0\|$$

gilt, folgt

$$\|F\| = \sup_{\substack{x\in X\\x\not= 0}} \frac{\|F x\|}{\|x\|} = \|y_0\|.$$

Wir kommen nun zum Rieszschen Darstellungssatz: Beschränkte lineare Funktionale eines Hilbertraumes $X$ lassen sich besonders einfach darstellen. Die Darstellung aus obigem Beispiel erfasst **alle** beschränkten linearen Funktionale. Es gilt

```{admonition} Satz: Darstellungssatz von Riesz
Sei $X$ ein Hilbertraum und $F\in X^*$ beliebig. Dann gibt es ein *eindeutig* bestimmtes $y\in X$, so dass $F$ die Darstellung

$$F x = (x,y)\quad\forall x\in X$$

besitzt.
```

**Bemerkung**: Dieser Satz ist das zentrale Ergebnis der Hilbertraum-Theorie. Neben seiner Bedeutung als Darstellungssatz kann er auch als Existenz- und Eindeutigkeitsprinzip aufgefasst werden. Diese Bedeutung des Rieszschen Satzes ist Grundlage für die moderne Theorie der elliptischen partiellen Differentialgleichungen.

