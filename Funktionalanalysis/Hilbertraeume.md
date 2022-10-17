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

# Skalarprodukträume. Hilberträume

Das aus der linearen Algebra bekannte Skalarprodukt lässt sich auch auf unendlich-dimensionale Räume übertragen. Wir definieren ganz allgemein

```{prf:definition} Skalarprodukt, Skalarproduktraum
:label: my-def-skalarprodukt

Unter einem *Skalarproduktraum* versteht man einen linearen Raum $V$ über $\mathbb{K}$, in dem ein *Skalarprodukt* $(x,y)$ mit folgenden Eigenschaften definiert ist: Für beliebige $x,y,z \in V$ und $\alpha\in\mathbb{K}$ ist

$$(\cdot, \cdot) : V \times V \to \mathbb{K}$$

und es gilt

$$\begin{split}
(x,x) & \ge 0, \quad (x,x) = 0\ \Leftrightarrow\ x=0\\
(x,y) & = \overline{(y,x)}\\
(\alpha x, y) & = \alpha (x,y)\\
(x+y,z) & = (x,z) + (y,z).
\end{split}$$ (eq:eigenschaftenskalarprodukt)
```

Beispielsweise lässt sich auf $V=C[a,b]$ Menge der reellwertigen stetigen Funktionen auf dem Intervall $[a,b]$ durch

$$(x,y) := \int_a^b x(t)\,y(t) dt\quad x,y \in V$$ (eq:L2SkalarProdukt)

ein Skalarprodukt definieren.

```{admonition} Aufgabe
Weise die Eigenschaften {eq}`eq:eigenschaftenskalarprodukt` für das Skalarprodukt {eq}`eq:L2SkalarProdukt` nach.
```

```{prf:theorem}
:label: my-thm-skal

Es sei $V$ ein Skalarproduktraum. Für beliebige $x,y,z\in V$ und $\alpha \in \mathbb{K}$ gilt

$$\begin{split}(x, \alpha y) & = \overline{\alpha} (x,y)\\
(x, y+z) & = (x,y) + (x,z)\end{split}$$
```

````{prf:proof}
Selber durchführen.
````

```{prf:theorem} induzierte Norm
:label: my-thm-indnorm

In jedem Skalarproduktraum $V$ lässt sich durch

$$\|x\| := \sqrt{(x,x)}$$

eine Norm definieren. Man bezeichnet sie als *die durch das Skalarprodukt* $(x,y)$ *induzierte Norm*.
```

Der Beweis lässt sich einfach mit Hilfe der *Schwarz'schen Ungleichung* durchführen.

```{prf:theorem} Schwarzsche Ungleichung
:label: my-thm-schwarz

Es sei $V$ ein Skalarproduktraum. Dann gilt für alle $x,y\in V$

$$|(x,y)| \le \|x\|\,\|y\|.$$
```

Die Beweise sind dem Leser überlassen (vgl. {cite:p}`burg_wille_haf_meister_PDE` S. 43, 44).

```{prf:theorem}
:label: my-thm-skalstetig

Sei $V$ ein Skalarproduktraum. Ferner seien $\{x_n\}$ und $\{y_n\}$ Folgen aus $V$ mit $x_n \to x$ und $y_n \to y$ für $n\to\infty$, wobei die Konvergenz im Sinne der induzierten Norm zu verstehen ist. Dann gilt

$$(x_n, y_n) \to (x,y)\quad\text{für}\quad n\to\infty,$$

dh. das Skalarprodukt ist eine stetige Funktion bezüglich der Normkonvergenz.
```

```{prf:definition} Hilbertraum
:label: my-def-Hilbertraum

Ein Skalarproduktraum $V$, der bezüglich der durch das Skalarprodukt induzierten Norm

$$\|x\| = \sqrt{(x,x)}\quad\text{für}\ x\in V$$
vollständig ist, heisst *Hilbertraum*.
```

**Beispiele**:

* $\mathbb{R}^n$ bzw. $\mathbb{C}^n$ sind mit den Skalarprodukte 

  $$(x,y) = \sum_{k=1}^n x_k y_k\quad\text{bzw.}\quad (x,y) = \sum_{k=1}^n x_k \overline{y_k}$$

  Hilberträume.

* $l_2$ ist mit dem Skalarprodukt

  $$(x,y) = \sum_{k=1}^\infty x_k \overline{y_k}$$

  ein Hilbertraum.

* $C[a,b]$ ist bezüglich der Quadratnorm

  $$\|x\|_2 = \sqrt{(x,x)} = \left(\int_a^b |x(t)|^2 dt\right)^{1/2}$$

  **kein** Hilbertraum, da $(C[a,b], \|\cdot\|_2)$ nicht vollständig ist. Der Raum $L_2[a,b]$ erweist sich als Vervollständigung dieses Raumes, welcher jedoch eine Verallgemeinerung des Riemannschen Integralbegriffs (dem Lebesgues Mass) erfordert.

Viele Eigenschaften des euklischen Raumes $\mathbb{R}^n$, die mit dem Skalarprodukt zusammenhängen, können auf einen beliebigen Hilbertraum übertragen werden. Der Begriff der Ortohogonalität ist dabei sehr zentral.

```{prf:definition} orthogonal
:label: my-def-orthogonal

Es sei $V$ ein Hilbertraum.

* Zwei Elemente $x,y \in V$ heissen *orthogonal* ($x\perp y$), wenn

  $$(x,y) = 0$$
  
  gilt.
* Zwei Teilmengen $A,B \subset V$ heissen *orthogonal* ($A\perp B$), wenn
  
  $$(x,y) = 0 \quad \forall\ x\in A, y\in B$$

  gilt.
* Sei $M\subset V$, dann heisst

  $$M^\perp := \{y\in V | (x,y)=0\quad\forall\,x\in M\}$$

  *Orthogonalraum* von $M$.
* Sei $V'$ ein abgeschlossener Unterraum von $V$. $V''$ wird *orthogonales Komplement* von $V'$ genannt, wenn

  $$V''\perp V' \quad \text{und}\quad V' \oplus V'' = V$$
  gilt. Mit $\oplus$ ist die direkte Summe bezeichnet.
```

Es gilt

```{prf:theorem} Pythagoras
:label: my-thm-pythagoras

Es sei $V$ ein Hilbertraum und seien $x,y \in V$ mit $x\perp y$. Dann gilt

$$\|x+y\|^2 = \|x\|^2 + \|y\|^2.$$
```
````{prf:proof}
Einfaches Nachrechnen.
````

```{prf:theorem}
:label: my-thm-unterraum

Es sei $V$ ein Hilbertraum und $M$ eine beliebige Teilmenge von $V$. Dann ist der Orthogonalraum $M^\perp$ von $M$ ein abgeschlossener Unterraum von $V$.
```

Für den Beweis sei auf {cite:p}`burg_wille_haf_meister_PDE` S. 50 verwiesen.

Wir kommen nun auf das Approximationsproblem aus {ref}`BestapproximationMetrisch` zurück und können die Struktureigenschaften des Hilbertraumes nutzen:

```{prf:theorem}
:label: my-thm-eigabUR

Es sei $V'$ ein abgeschlossener Unterraum von $V$ und $x\in V$ beliebig. 
* Dann existiert genau ein $x_0 \in V'$ mit

  $$\|x-x_0\| = \min_{x'\in V'} \|x-x'\|,$$

  dh. zu jedem $x\in V$ gibt es genau ein bestapproximierendes Element bezüglich $V'$.
* Es gilt $x_0\in V$ ist genau dann bestapproximierend an $x\in V$, wenn

  $$(x-x_0,y) = 0\quad \forall y\in V'$$
  gilt.
```

Für den Fall, dass $V'$ endlich-dimensional ist, lässt sich das bestapproximierende Element konstruieren. Es gilt

```{prf:theorem}
:label: my-thm-bestapprox

Sei $V'\subset V$ mit $\dim V' < \infty$ ein Unterraum des Hilbertraumes $V$ und sei $x_1, \ldots, x_n$ eine Basis von $V'$. Dann lässt sich das eindeutig bestimmte bestapproximierende Element $x_0\in V'$ an $x\in V$ in der Form

$$x_0 = \sum_{k=1}^n \lambda_k x_k$$

darstellen, wobei die Koeffizienten $\lambda_1, \ldots, \lambda_n$ durch das lineare Gleichungssystem

$$(x-x_0,x_i) = 0 \quad i = 1, \ldots, n$$

und mit $x_0$ eingesetzt

$$(x,x_i) - \sum_{k=1}^n \lambda_k (x_k,x_i) = 0 \quad i = 1, \ldots, n$$ (eq:bestapproxHilbert)

gegeben sind.
```

```{prf:remark}
:label: my-rm-orthogonalsystem

Bildet $x_k$, $k=1,\ldots, n$ ein *Orthonormalsystem*

$$(x_i,x_k) = \delta_{i k} = \begin{cases}
0\quad \text{für}\ i\not= k\\
1\quad \text{für}\ i = k\end{cases},$$
so folgt aus {eq}`eq:bestapproxHilbert` sofort

$$\lambda_i = (x,x_i)\quad i = 1, \ldots,n.$$

Die Koeffzienten $\lambda_i$ nennt man auch Fourierkoeffizienten!
```

Wir werden zeigen, dass sich mit Hilfe des Schmidtschen Orthogonalisierungsverfahrens aus $n$ linear unabhängigen Elementen stets ein Orthogonalsystem konstruieren lässt.

```{prf:definition} Orthonormalsystem (ONS)
:label: my-def-ONS

Es sei $V$ ein Hilbertraum. Man nennt die Folge $\{x_k\}_{k\in\mathbb{N}}$ ein (abzählbares) *Orthonormalsystem* (kurz ONS) von $V$, wenn 

$$(x_j, x_k) = \delta_{i k} = \begin{cases}
0\quad \text{für}\ i\not= k\\
1\quad \text{für}\ i = k\end{cases}\quad \text{für alle}\ j,k\in\mathbb{N}$$

erfüllt ist.
```

**Beispiele**
* Im $l_2$ bilden die Folgen $\{1,0,0,\ldots\}$, $\{0,1,0,\ldots\},\ldots $ ein ONS.
* Für den (nicht vollständigen) reellen Skalarproduktraum $C[0,2\pi]$ bilden die trigonometrischen Funktionen

  $$\frac{1}{\sqrt{2\pi}}, \frac{1}{\sqrt{\pi}} \cos t, \frac{1}{\sqrt{\pi}} \sin t, \frac{1}{\sqrt{\pi}} \cos 2t, \frac{1}{\sqrt{\pi}} \sin 2t, \ldots $$
  ein ONS.
* Weitere Beispiele sind Hermitesche und Legendresche Polynome.

Es gilt ganz allgemein:

```{prf:theorem} Orthogonalisierungsverfahren nach Erhard Schmidt
:label: my-thm-ErhardSchmidt

Gegeben sei eine abzählbare (nicht zwingend endliche) linear unabhängige Folge $\{y_k\} \subset V$ aus dem Hilbertraum $V$. Dann gibt es ein ONS aus $n$ bzw. abzählbar unendlich viele Elementen $\{x_k\}$ so, dass der von der Folge $\{y_k\}$ aufgespannte (abgeschlossene) Unterraum mit dem der Folge $\{x_k\}$ aufgespannten (abgeschlossene) Unterraum übereinstimmt.
```

Der Beweis beruht auf der Konstruktion: sei

$$x_1 = \frac{y_1}{\|y_1\|}$$

so ist $\mathop{span}\{y_1\} = \mathop{span}\{x_1\}$. Wir nehmen nun an, es seien bereits $k$ orthonormierte Elemente $x_1, \ldots, x_k$ mit $\mathop{span}\{y_1, \ldots, y_k\} = \mathop{span}\{x_1,\ldots, x_k\}$ konstruiert. Dann setze

$$z_{k+1} = y_{k+1} - \sum_{j=1}^k (y_{k+1},x_j) x_j.$$

In dem Fall gilt $(z_{k+1},x_i) = 0$ für alle $i = 1, \ldots, k$. Mit

$$x_{k+1} = \frac{z_{k+1}}{\|z_{k+1}\|}$$

folgt $\mathop{span}\{y_1, \ldots, y_{k+1}\} = \mathop{span}\{x_1,\ldots, x_{k+1}\}$.

```{prf:remark}
:label: my-rm-numschmidt

Das Verfahren wird auch in der Numerik angewandt als wichtiges Beispiel sei die Arnoldi Iteration für die Berechnung von **Eigenwerten** erwähnt. Ebenso kann das Verfahren für die Berechnung der **QR-Zerlegung** benutzt werden. Wobei zu erwähnen ist, dass die Householder-Transformation numerisch stabiler ist.
```

:::{seealso}
[Beispiel für das Orthogonalisierungsverfahren.](BeispielOrthogonalisierungSchmidt.ipynb)
:::

Die Elemente eines Hilbertraumes können mit Hilfe eines vollständigen Orthogonalsystems dargestellt werden. Dies gelingt mit Hilfe der _verallgemeinerten Fourierreihen_:

```{prf:theorem} Fourierentwicklung
:label: my-thm-Fourierentwicklung

Sei $V$ ein Hilbertraum mit einem vollständigen ONS $\{x_k\}_{k\in\mathbb{N}}$.
* Dann lässt sich jedes $x\in V$ in der Summenform

  $$x = \sum_{k=1}^\infty a_k\,x_k\quad\text{Fourierentwicklung von $x$}$$

  mit eindeutig bestimmten Koeffizienten

  $$a_k = (x,x_k)\in\mathbb{K}$$

  darstellen und die Reihe $\sum_{k=1}^\infty |a_k|^2$ ist konvergent.
* Umgekehrt gibt es zu jeder Zahlenfolge $\{a_k\}_{k\in\mathbb{N}}$ in $\mathbb{K}$, für die $\sum_{k=1}^\infty |a_k|^2$ konvergiert, genau ein $x\in V$ mit $x=\sum_{k=1}^\infty a_k x_k$.
```

```{prf:remark}
:label: my-rm-hilbertbasis

Aufgrund der Darstellung $x = \sum_{k=1}^\infty a_k\,x_k$ nennt man ein vollständiges ONS auch eine *Hilbertraumbasis*.
```

```{prf:theorem} Struktur von Hilberträumen
:label: my-thm-StrukturHR

Es sei $V$ ein Hilbertraum und $\{x_k\}_{k\in\mathbb{N}}$ ein (abzählbares) ONS in $V$. Dann sind die folgenden Aussagen äquivalent:
* $V = \overline{\underset{k\in\mathbb{N}}{\bigoplus} \mathop{span}(x_k)}$.
* Das ONS $\{x_k\}_{k\in\mathbb{N}}$ ist abgeschlossen.
* Für alle $x\in V$ gilt die *Parsevalsche Gleichung*

  $$\sum_{k=1}^\infty |(x,x_k)|^2 = \|x\|^2\quad\text{(Vollständigkeitsrelation)}$$

* Jedes Element $x\in V$ besitzt die Fourierentwicklung

  $$x = \sum_{k=1}^\infty (x,x_k)\,x_k.$$
```

:::{seealso}
[Beispiel zur Fourierentwicklung.](BeispielFourierEntwicklung.ipynb)
:::
