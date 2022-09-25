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

## Normierte Räume. Banachräume

Bis jetzt haben wir sehr wenig Eigenschaften eines Raumes benötigt. Was uns noch fehlt, sind abgesehen vom Abstand der Elemente noch *algebraische* Eigenschaften (addieren, multiplizieren, etc.). Dazu definieren wir den *linearen Raum* (oder *Vektorraum*) wie folgt.

```{prf:definition} linearer Raum
:label: my-def-linraum

Ein *linearer Raum* (oder *Vektorraum*) über einem Körper $\mathbb{K}$ besteht aus einer nichtleeren Menge $V$, sowie

* einer Vorschrift, die jedem Paar $(x,y)$ mit $x,y \in V$ genau ein Element $x+y\in V$ zuordnet (*Addition*)
* einer Vorschrift, die jedem Paar $(\lambda,x)$ mit $\lambda\in \mathbb{K}$ und $x \in V$ genau ein element $\lambda x\in V$ zuordnet (*Multiplikation mit Skalaren*), wobei für alle $x,y,z \in V$ und $\lambda, \mu\in\mathbb{K}$ folgende Regeln gelten:

$$\begin{array}{rcll}
x + (y+z) & = & (x+y) + z & \quad\text{Assoziativgesetz} \\
x + y & = & y + x & \quad\text{Kommutativgesetz} \\
x + 0 & = & x & \quad\text{Nullelement} \\
x + x' & = & 0 & \quad\text{Negatives zu $x$}\\
(\lambda + \mu) x & = &  \lambda x + \mu x & \quad\text{1. Distributivgesetz}\\
\lambda (x+y) & = &  \lambda x + \mu x & \quad\text{2. Distributivgesetz}\\
(\lambda \mu) x & = & \lambda (\mu x) & \quad\text{Assoziativgesetz}\\
1 x & = & x & \quad\text{mit $1\in\mathbb{K}$}
\end{array}$$
```

Beispiele für lineare Räume:

* Die Mengen $\mathbb{R}^n$, $\mathbb{C}^n$ sind wohl bekannt.
* Menge $C[a,b]$ aller reellwertigen *stetigen* Funktionen:
  
  $$\begin{split}
  (x+y)(t) & = x(t) + y(t)\\
  (\lambda x)(t) & = \lambda x(t)\\
  0(t) & = 0\\
  (-x)(t) & = -x(t) \end{split}$$

  mit $t \in [a,b]\subset \mathbb{R}$, $\lambda\in\mathbb{R}$.
* $C^k[a,b]$ Menge aller reellwertigen $k$-mal stetig differenzierbare Funktionen.
* $C^{\infty}[a,b]$ Menge aller beliebig oft stetig differenzierbare Funktionen.
* Menge aller Polynome
* Menge $l_p$ aller Zahlenfolgen $x = \{x_k\}$, für die $\sum_{k=1}^{\infty} |x_k|^p < \infty$ konvergiert:

  $$\begin{split}
  x+y & = \{x_k\} + \{y_k\} = \{x_k + y_k\}\\
  \lambda x & = \lambda \{x_k\} = \{\lambda x_k\},\quad \lambda\in\mathbb{R}\end{split}$$
  
Wie in der linearen Algebra sind folgende Begriffe analog definiert

```{prf:definition} Unterraum, lineare Mannigfaltigkeit, lineare Hülle / Span, linear unabhängig, Dimension, Basis
:label: my-def-URMannigf

* Eine nicht leere Teilmenge $S$  von $V$ heisst *Unterraum* von $V$, wenn für bel. $x,y \in S$  und $\lambda \in \mathbb{K}$ stets

  $$x+y \in S\quad\text{und}\quad \lambda x \in S$$

  folgt. Insbesondere ist $S$ selbst ein linearer Raum über $\mathbb{K}$.
* Ist $S$ ein Unterraum von $V$ und $x_0\in V$ beliebig, so nennt man

  $$M = x_0 + S := \{x_0+y\ |\ y\in S\}$$
  eine *lineare Mannigfaltigkeit* von $V$.
* Ist $A$ eine beliebige nichtleere Teilmenge von $V$, so bilden alle Linearkombinationen $\sum_{k=1}^m \lambda_k x_k$ mit beliebigem $m \in \mathbb{N}$, $\lambda_k\in\mathbb{K}$, $x_k \in A$ einen Unterraum $S\subset V$. Er wird *lineare Hülle von* $A$ oder $\mathop{span} A := S$ genannt.

  Man sagt: $A$ spannt $S$ auf oder $A$ ist ein *Erzeugendensystem* von $S$. Im Falle $S=V$ spannt $A$ den ganzen Raum $V$ auf.
* Sind $S$ und $T$ Unterräume von $V$, dann ist die *Summe* $S+T$ definiert durch

  $$S+T := \mathop{span} S \cup T.$$
* Die (endlich vielen) Elemente $x_1, \ldots, x_n\in V$ heissen *linear unabhängig*, wenn aus

  $$\alpha_1 x_1 + \ldots + \alpha_n x_n = 0\quad\text{stets}\quad \alpha_1 = \ldots = \alpha_n = 0$$

  folgt.
* Sei $S$ ein Unterraum von $V$. Wir sagen, $S$ besitzt die *Dimension* $n$, $\mathop{dim} S = n$, wenn es $n$ linear unabhängige Elemente von $S$ gibt, aber $n+1$ Elemente von $S$ stets linear abhängig sind.

  $S$ heisst *Basis* von $V$, wenn die Elemente von $S$ linear unabhängigsind und $\mathop{span} S = V$ gilt.

  Besitzt $V$ keine endlich dimensionale Basis, nennt man $V$ *unendlich-dimeansional* ($\mathop{dim} V = \infty$).
```

```{prf:remark}
:label: my-rm-unendlichdimfuncraeume

Die oben erwähnten Funktionenräume $C[a,b]$, $C^k[a,b]$, Polynome sind unendlich-dimensional. Ebenso ist der Folgenraum $l_p$ unendlich-dimensional:

> Man betrachte
>
> $$x^{(k)} = \{0, \ldots, 0, 1, 0, \ldots \}\in l_p$$
> mit 1 an der Stelle $k$.
```

Im folgenden sind wir an Räumen interessiert, für welche eine lineare Struktur und eine Metrik gegeben ist.

```{prf:definition} normierter Raum
:label: my-def-normRaum

Sei $V$ ein metrischer und linearer Raum. Zu dem sei die Metrik $d$ von $V$ *translationsinvariant*

$$d(x+z, y+z) = d(x,y)\quad\forall\ x,y,z\in V$$
und *homogen*

$$d(\alpha x, \alpha y) = |\alpha| d(x,y)\quad \forall\ \alpha\in\mathbb{K}, x,y\in V.$$
Dann nennt man $V$ einen *normierten Raum*. Der durch

$$\|x\| := d(x,0)\quad \forall\ x\in V$$
erklärte Ausdruck heisst *Norm* von $x$.
```

```{prf:remark}
:label: my-rm-grmetrik

* Neben der kurzen Schreibweise $V$, verwendet man häufig auch die Bezeichnung $(V, \|\cdot\|)$. Der Punkt in $\|\cdot\|$ ist als Platzhalter zu verstehen.
* Führt man den normierten Raum $V$ mit Hilfe einer Norm ein, so ist durch

  $$d(x,y) = \|x-y\|\quad\text{für alle}\ x,y\in V$$
  eine Metrik in $V$ gegeben.
```

```{admonition} Folgerung
Ein normierter Raum $(V, \|\cdot\|)$ ist ein linearer Raum, auf dem eine Norm $\|\cdot\|$ erklärt ist, die für alle $x,y\in V$ und $\alpha \in \mathbb{K}$

$$\begin{split}
\|x\| & \ge 0,\quad \|x\|=0\quad \text{genau dann, wenn $x=0$ ist},\\
\|\alpha x\| & = |\alpha| \|x\|\\
\|x + y\| & \le \|x\| + \|y\|\end{split}$$
erfüllt.
```

Damit können wir einen wichtigen Begriff der Funktionalanalysis einführen, den Banachraum:

```{prf:definition} Banachraum
:label: my-def-Banachraum

*Vollständig normierte* Räume $V$ sind diejenigen, für die jede Cauchy-Folge in $V$ gegen ein Element in $V$ konvergiert.

Ein vollständiger normierter Raum heisst *Banachraum*.
```

**Beispiele**: Folgende Räume sind Banachräume

* $\mathbb{R}^n$ mit der Norm $\|x\| = \sqrt{\sum_{k=1}^n |x_k|^2}$.
* $C[a,b]$ mit der Norm $\|x\| := \max_{a\le t \le b} |x(t)|$.
* $C^k[a,b]$ mit der Norm $\|x\| := \max_{a\le t \le b} |x(t)| + \max_{a\le t \le b} |x'(t)| + \ldots + \max_{a\le t \le b} |x^{(k)}(t)|.$  
* $l_p\ (1\le p < \infty)$ mit der Norm $\|x\| = \left(\sum_{k=1}^{\infty} |x_k|^p \right)^{1/p}$

```{prf:definition} abzählbare Basis
:label: my-def-abzBasis

Man sagt, dass $V$ eine *abzählbare Basis* $\{x_k\}_{k=1}^{\infty}$ mit $x_k \in V$ besitzt, falls jedes $x\in V$ eindeutig in der Form $x = \sum_{k=1}^{\infty} \alpha_k\,x_k$ darstellbar ist, wobei die Konvergenz bezüglich der Norm von $V$ zu verstehen ist.
````

Lineare Räume können durchaus verschieden normiert werden. Als Beispiel betrachte dazu den Raum $V=\mathbb{R}^n$ mit den Normen

$$\begin{split}
\|x\|_1 & = \sum_{k=1}^n |x_k|\quad\text{(Betragsnorm)}\\
\|x\|_2 & = \sqrt{\sum_{k=1}^n |x_k|^2}\quad\text{(Euklidische Norm, Quadratnorm)}\\
\|x\|_{\infty} & = \max_{1\le k \le n} |x_k|\quad\text{(Maximumsnorm)}\end{split}$$

```{admonition} Aufgabe
Stelle den Einheitskreis $K_{*} = \{x\in\mathbb{R}^n\, \big|\, \|x\|_{*} = 1\}$ bezüglich den drei verschiedenen Normen $*$ graphisch dar.
```

```{code-cell} ipython3
:tags: [hide-cell]

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

x,y = np.meshgrid(np.linspace(-1.2,1.2,81),np.linspace(-1.2,1.2,81))

z1 = np.abs(x)+np.abs(y)
z2 = np.sqrt(x**2+y**2)
zinf = np.max([np.abs(x),np.abs(y)],axis=0)

fig, (ax1, ax2, ax3)  = plt.subplots(1, 3)
ax1.contour(x,y,z1-1,0,colors='tab:blue')
ax2.contour(x,y,z2-1,0,colors='tab:orange')
ax3.contour(x,y,zinf-1,0,colors='tab:green')
for a,t in zip((ax1, ax2, ax3),('$K_1$','$K_2$','$K_\infty$')):
    a.set_aspect(1)
    a.grid()
    a.set_title(t)
plt.tight_layout()
plt.show()
```

Die drei verschiedenen Normen führen zur Frage, wie die Normen zusammenhängen.

```{admonition} Defintion: äquivalente Normen
Zwei Normen $\|\cdot\|_a$ und $\|\cdot\|_b$ heissen *äquivalent*, wenn jede bezüglich der Norm $\|\cdot\|_a$ konvergente Folge auch bezüglich der Norm $\|\cdot\|_b$ konvergent ist und umgekehrt.
```

Im endlich dimensionalen gilt ein pauschaler Satz:

```{prf:theorem}
:label: my-thm-normaequiv

Alle Normen in einem **endlich**-dimensionalen Raum $V$ sind äquivalent.
```

Dieses Resultat gilt für unendlich-dimensionale Räume **nicht**. Als Beispiel sei der lineare Raum $C[a,b]$ mit der Maximumsnorm und der Quadratnorm erwähnt. Die beiden Normen sind nicht äquivalent, vgl. das Gegenbeispiel {eq}`eq:Gegenbeispiel`. 

Mit diesem Satz folgt

```{prf:theorem}
:label: my-thm-enddimbanach

Jeder endlich-dimensionale normierte Raum $V$ ist vollständig, also ein Banachraum.
```