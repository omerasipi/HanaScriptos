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

# Beschränkte lineare Operatoren

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
:label: my-rm-opnorm

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
:label: my-rm-stetigbeschraenkt

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
