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

# Grundlegende Räume

## Vektorräume

In den Modulen linearen Algebra, Analysis und Numerik wurden Vektorräume und Normen für den $\mathbb{R}^n$ eingeführt und benutzt. Wir definieren hier der Vollständigkeit halber die Begriffe nochmals und erweitern die Anwendung auf allgemeinere Räume, insbesondere müssen diese nicht endlich dimensional sein.

Beginnen wir mit dem Begriff der Vektorräume. Um diesen definieren zu können, benötigen wir einen Zahlenkörper. In aller Regel benutzen wir die reellen Zahlen $\mathbb{R}$. Wir legen mit dem folgenden Axiom fest, was die reellen Zahlen sind:

```{admonition} Definition: reelle Zahlen

Es existiere eine Menge $\mathbb{R}$ mit folgenden Eigenschaften:
- Es existieren Operationen

  $$\begin{split}  + : \mathbb{R} \times \mathbb{R} \ & \to \ \mathbb{R}\quad\text{Addition}\\
  \cdot : \mathbb{R} \times \mathbb{R} \ & \to \ \mathbb{R}\quad\text{Multiplikation}\end{split}$$
  mit den Eigenschaften:
  - Assoziativgesetze

    $$\begin{split}
    (a+b) + c & = a + (b+c)\quad \forall\ a,b,c \in \mathbb{R}\\
    (a\cdot b) \cdot c & = a \cdot (b\cdot c)\quad \forall\ a,b,c \in \mathbb{R}
    \end{split}$$
 
  - Kommutativgesetze
    
    $$\begin{split}
    a + b & = b + a \quad \forall\ a,b \in \mathbb{R}\\
    a\cdot b & = b \cdot a\quad \forall\ a,b \in \mathbb{R}
    \end{split}$$

  - Neutrale Elemente
    
    $$\begin{split}
    \exists 0: a + 0 = a\quad \forall a \in \mathbb{R}\\
    \exists 1: a \cdot 1 = a\quad \forall a \in \mathbb{R}
    \end{split}$$
  
  - Inverse Elemente
    
    $$\begin{split}
    \forall a\in \mathbb{R}\ \exists a'\in\mathbb{R} : a + a' = 0\\
    \forall a\in \mathbb{R}\setminus \{ 0\}\ \exists \tilde{a}\in\mathbb{R} : a \cdot \tilde{a} = 1
    \end{split}$$
    Schreibweise: $-a := a', a - a := a + (-a), \frac{1}{a} := \tilde{a}, \frac{a}{a} := a\cdot \frac{1}{a}$
  
  - Distributivgesetz

    $$a\cdot (b+c) = a\cdot b + a\cdot c$$

- $\mathbb{R}$ ist geordnet, d.h. es existiert eine Relation $\le$ so, dass gilt
  - $\mathbb{R}$ ist totalgeordnet, d.h.
    1. $\mathbb{R}$ ist teilgeordnet, d.h.
       - Für alle $x\in\mathbb{R}$ gilt: $x \le x$.
       - Ist $x \le y$ und $y\le z$ für $x,y,z\in\mathbb{R}$, so ist $x \le z$.
       - Ist $x \le y$ und $y \le x$ für $x,y \in \mathbb{R}$, so ist $x = y$.
       
       (Mit $x<y$ bezeichnen wir den Fall $x \le y$ und $x \not= y$.)
    2. Für je zwei Elemente $x,y \in\mathbb{R}$ gilt
       
       $$x \le y\quad \text{oder}\quad y\le x.$$
  - Die Ordnung ist verträglich mit Addition und Multiplikation, d.h. für $a,b,c\in\mathbb{R}$ gilt
  
  $$\begin{array}{c}
  a \le b \Rightarrow a+c \le b+c\\
  a \le b, 0 \le c \Rightarrow a\cdot c \le b\cdot c
  \end{array}$$

- $\mathbb{R}$ ist vollständig. Das heisst, dass jede nicht leere nach oben beschränkte Menge reeller Zahlen eine kleinste obere Schranke besitzt.
```

Im Allgemeinen ist ein (Zahlen)-Körper über eine Gruppe wie folgt definiert:

```{admonition} Definition: Gruppe
Ein Tupel $(G, \cdot)$ bestehend aus einer Menge $G$ und einer Verknüpfung $\cdot : G \to G$ heisst *Gruppe*, falls die Verknüpfung assoziativ ist, ein neutrales Element $e\in G$ existiert und für alle $a \in G$ ein $b\in G$ exisitert so dass $a\cdot b = b\cdot a = e$ gilt. Ist die Verknüpfung kommutativ, nennt man die Gruppe kommutativ.
```

Ein Körper lässt sich somit wie folgt allgemein definieren:

```{admonition} Definition: Körper
Ein *Körper* ist eine Tripel $(K, +, \cdot)$ mit folgenden Eigenschaften:

- $(K, +)$ ist eine kommutative Gruppe mit neutralem Element $0_K$.
- $(K, \cdot)$ ist eine kommutative Gruppe mit neutralem Element $1_K$.
- (Distributivgesetz) Für alle $a,b,c\in K$ gilt

$$\begin{split}a\cdot (b+c) & = a\cdot b + a\cdot c\\
(a+b)\cdot c & = (a\cdot c) + (b\cdot c)\end{split}$$
```

Beispiele für Körper sind folgende Tupel

- $(\mathbb{R}, +, \cdot)$ und $(\mathbb{Q}, +, \cdot)$  mit der üblichen Addition und Multiplikation,
- $(\mathbb{Z}, +, \cdot)$ bildet **keine** Gruppe.

Mit Hilfe eines Körpers können wir nun einen $\mathbb{K}$-Vektorraum wie folgt definieren:

```{admonition} Definition: $\mathbb{K}-Vektorraum$
Ein $\mathbb{K}$-Vektorraum ist ein Tripel $(V, +, \cdot)$ mit den Eigenschaften

1. $(V,+)$ ist eine kommutative Gruppe
2. Die Abbildung $\cdot : \mathbb{K} \times V \to V$ genügt den Eigenschaften

    - $\alpha\cdot (\beta\cdot x) = (\alpha\cdot \beta)\cdot x$
    - $\alpha\cdot (x+y) = (\alpha\cdot x) + (\alpha\cdot y)$
    - $(\alpha+\beta)\cdot x = (\alpha\cdot x) + (\beta\cdot x)$
    - $1_K \cdot x = x$

    wobei $\alpha, \beta\in\mathbb{K},\ x, y\in V$. 
```

```{admonition} Bemerkung
Wenn man in einem allgemeinen Kontext von Vektoren spricht, so meint man damit die Elemente eines Vektorraumes. Spricht man von Skalaren, so sind die Elemente des zugrundeliegenden Körpers gemeint.
```

**Beispiele**:

*  Vektorraum $\mathbb{R}^n$

   Anwendung in python:

````{code-cell} ipython3
import numpy as np
x = np.array([1,2,3,4])
y = np.array([-4,-3,-2,-1])

print('x+y=',x+y)
print('5*(x+y)=',5*(x+y),'= 5*x+5*y = ',5*x+5*y)
````
*  Vektorraum der stetigen Funktionen.

   Sei $\alpha \in \mathbb{R}$ und $u,v : [a,b]\subset \mathbb{R} \to \mathbb{R}$ zwei stetige Funktionen. Zeige, dass die Summe zweier stetiger Funktionen $(u+v)$ wieder stetig ist und dass das Vielfache einer stetigen Funktion $(\alpha u)$ ebenso stetig ist.

   ```{admonition} Aufgabe
   Beweise die Aussage.
   ```

## Metrische Räume

### Metrik

In der Analysis will man oft eine Distanz zwischen zwei Elemente eines Vektorraumes angeben. Insbesondere bei Konvergenzbetrachtungen ist der Abstand zweier Elemente existentiell wichtig. Der Konvergenzbegriff in $\mathbb{R}$ unter Hinzunahme der folgenden Distanzfunktion

$$\begin{split} d : \mathbb{R}\times\mathbb{R} & \to \mathbb{R}^+\\
(x,y) & \mapsto  d = d(x,y) := |x-y|\end{split}$$

lautet wie folgt: Die Folge $\{x_n\}$ aus $\mathbb{R}$ heisst konvergent gegen $x_0\in\mathbb{R}$, wenn es zu jedem $\varepsilon > 0$ eine natürliche Zahl $n_0 = n_0(\varepsilon)$ gibt, so dass

$$d(x_n,x_0) = |x_n-x_0| < \varepsilon$$

für alle $n \ge n_0$ gilt.

```{admonition} Definition: Metrik

Eine nichtleere Menge $X$ mit *Elemente* $x, y, z, \ldots$ heisst ein *metrischer Raum*, wenn jedem Paar $x, y \in X$ eine reelle Zahl $d(x,y)$, genannt *Abstand* oder *Metrik*, zugeordnet ist, mit den Eigenschaften: Für alle $x,y,z\in X$ gilt

1. $d(x,y) \ge 0,\ d(x,y) = 0$ genau dann, wenn $x=y$ ist
2. $d(x,y) = d(y,x)$ *Symmetrieeigenschaft*
3. $d(x,y) \le d(x,z) + d(z,y)$ *Dreiecksungleichung*.
```

Für metrische Räume verwenden wir wieder die Schreibweisen: $(X, d)$ oder kurz $X$, falls der Kontext klar ist.

```{admonition} Aufgabe
Zeige, dass für $X = C[0,1]$, den stetigen Funktionen auf dem Intervall $[0,1]$ die Abbildung

$$\begin{split}d_{\text{max}} : X \times X & \to \mathbb{R}^+\\
(x,y) & \mapsto d_{\text{max}}(x,y) := \max_{t\in[0,1]} |x(t)-y(t)|\end{split}$$ (eq:maxnorm)

eine Metrik definiert (*Maximumsmetrik*).
```

**Beispiel**:

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

x = lambda t: t**2
y = lambda t: t*(1-t)+2/(1+((t-0.5)/0.02)**2)

t = np.linspace(0,1,400)

plt.plot(t,x(t), label='$x(t)$')
plt.plot(t,y(t), label='$y(t)$')
plt.legend()
plt.show()
```

Die Maximummetrik berechnet die maximale Differenz der beiden Funktionen:

```{code-cell} ipython3
:tags: [hide-input]

plt.plot(t,np.abs(x(t)-y(t)), label='$|x(t)-y(t)|$')
plt.legend()
plt.show()
```

```{admonition} Aufgabe

Berechne analytisch den exakten Wert der Maximumsmetrik für die beiden Funktionen $x,y$ aus obigem Beispiel.
```

Eine weitere wichtige Metrik ist die **Integralmetrik**. Es sei $X$ die Menge aller reellwertigen Funktionen, die auf einem (nicht notwendig beschränkten) Intervall $(a,b)$ stetig sind und für die das Integral

$$\int_a^b |x(t)|^p dt, \quad 1\le p < \infty$$

im Riemannschen Sinne existiert. Setzen wir für $x(t), y(t)\in X$

```{admonition} Definition: Integralmetrik
$$d_p(x,y) := \left(\int_a^b |x(t)-y(t)|^p dt\right)^{1/p},\quad 1\le p < \infty$$
```

so wird $X$ mit dieser *Integralmetrik* zu einem metrischen Raum. Der Beweis nutzt die Minkowski-Ungleichung für Integrale

$$\left(\int_a^b |u(t)-v(t)|^p dt\right)^{1/p} \le \left(\int_a^b |u(t)|^p dt\right)^{1/p} + \left(\int_a^b |v(t)|^p dt\right)^{1/p}.$$

**Beispiel**:

```{code-cell} ipython3
:tags: [hide-input]

plt.plot(t,x(t), label='$x(t)$')
plt.plot(t,y(t), label='$y(t)$')
plt.fill_between(t,x(t),y(t),alpha=0.3, label='$x(t)-y(t)$')
plt.legend()
plt.show()
```

```{admonition} Aufgabe
Berechne die Integralmetrik für $p=2$ und $p=1/2$ für das Beispiel oben.
```

```{admonition} Aufgabe
In der Codierungstheorie ist ein $n$-stelliges Binärwort ein $n$-Tuppel $(\xi_1, \ldots, \xi_n)$, wobei $\xi_k \in \{0,1\}$ für $k=1,\ldots, n$. Bezeichne $X$ die Menge aller dieser Binärwörter. Für $x = \xi_1 \xi_2 \ldots \xi_n$, $y = \eta_1 \eta_2 \ldots \eta_n$ ist die *Hamming-Distanz* zwischen $x$ und $y$ durch

$$d_H(x,y) := \text{Anzahl der Stellen an denen sich $x$ und $y$ unterscheiden}$$

definiert. Zeige:

1. $d_H(x,y)$ lässt sich durch

   $$d_H(x,y) := \sum_{k=1}^n [(\xi_k+\eta_k) \mod 2]$$

   darstellen.

2. $(X,d_H)$ ist ein metrischer Raum.
```

Die topologischen Begriffe wie *offene Kugel*, *innerer Punkt*, *Häufungspunkt*, *abgeschlossen*, *beschränkt* lassen sich mit Hilfe der Metrix $d$ für einen metrischen Raum $(X,d)$ direkt aus dem aus der Analysis bekannten $\mathbb{R}^n$ übertragen.

```{admonition} Definition: Konvergenz
Eine Folge $\{x_n\}\subset X$ von Elemente aus $X$ heisst *konvergent*, wenn es ein $x_0\in X$ gibt mit

$$d(x_n,x_0) \to 0\quad \text{für}\quad n\to\infty,$$

dh. wenn es zu jedem $\varepsilon > 0$ ein $n_0 = n(\varepsilon)\in\mathbb{N}$ gibt, mit

$$d(x_n,x_0) < \varepsilon\quad\text{für alle}\quad n\ge n_0.$$

Schreibweise: $x_n \to x_0$ für $n \to \infty$ oder $\lim_{n\to\infty} x_n = x_0$, $x_0$ heisst *Grenzwert* der Folge $\{x_n\}$
```

**Der Grenzwert einer konvergenten Folge ist eindeutig bestimmt.** Dies lässt sich per Widerspruch wie folgt leicht zeigen. Seien $x_0$ und $y_0$ zwei verschiedene Grenzwerte. Daher gilt

$$\begin{split}
0 < d(x_0,y_0) & \le d(x_0, x_n) + d(x_n,y_0)\\
& = \underbrace{d(x_n,x_0)}_{\to 0\ \text{für}\, n\to 0} + \underbrace{d(x_n,y_0)}_{\to 0\ \text{für}\, n\to 0} \to 0\ \text{für}\, n\to 0.\end{split}$$

Es folgt damit $d(x_0,y_0)=0$ und damit $x_0 = y_0$.

### Funktionenfolgen

Wir starten mit einem intuitiven Begriff der Konvergenz für Funktionenfolgen:

```{admonition} Definition: Punktweise Konvergenz
Man nennt eine Funktionenfolge $\{x_n(t)\} \subset C[a,b]$ *punktweise konvergent*, wenn für jedes $t\in [a,b]$ die Zahlenfolge $x_1(t), x_2(t), \ldots $ konvergiert. Die *Grenzfunktion* $x$ ist dabei durch

$$\lim_{n\to\infty} x_n(t) = x(t)\quad\text{für jedes}\quad t\in [a,b]$$

gegeben.
```

Die punktweise Konvergenz erweist sich für die Analysis als zu schwach. Als Beispiel dazu betrachten wir die Folge der Funktionen $x_n \subset C(\mathbb{R})$

$$x_n(t) = \frac{1}{1+x^{2n}},\quad n=1,2,3,\ldots$$

Wie im Python Code unten leicht zu sehen ist, konvergiert die Folge punktweise gegen

$$x(t) = \begin{cases}
1,\quad \text{für}\ |t| < 1,\\
1/2, \quad \text{für}\ |t| = 1,\\
0, \quad \text{für}\ |t| > 1.\end{cases}$$

Obwohl alle Funktionen $x_n$ stetig sind, ist die Grenzfunktion $x$ unstetig und damit nicht in unserem Funktionenraum (oder Vektorraum) $x\not\in C(\mathbb{R})$! Das ist für uns unbrauchbar.

```{code-cell} ipython3
:tags: [hide-input]

def x(t):
    y = np.zeros_like(t)
    y[np.abs(t)<1] = 1
    y[np.abs(t)==1] = 0.5
    return y

xn = lambda t,n : 1/(1+t**(2*n))
t = np.linspace(-3,3,400)

plt.plot(t,xn(t,1), label='$xn=1$')
plt.plot(t,xn(t,4), label='$xn=4$')
plt.plot(t,xn(t,8), label='$xn=8$')
plt.plot(t,x(t),'--', label='Grenzfunktion')
plt.legend()
plt.show()
```

Das Problem der punktweisen Konvergenz besteht darin, dass für jedes $t\in\mathbb{R}$ eine eigene Schranke $\varepsilon = \varepsilon(t)$ gewählt werden kann. Dies ist zu schwach, um garantieren zu können, dass eine konvergente Funktionenfolge wieder stetig ist. Wir benutzen nun die Maximumsmetrik {eq}`eq:maxnorm`. Daraus folgt, dass zu jedem $\varepsilon > 0$ es ein $n_0 = n_0(\varepsilon) \in \mathbb{N}$ mit

$$\max_{t\in\mathbb{R}} |x_n(t) - x(t)| < \varepsilon\quad \text{für alle}\ n \ge n_0$$

geben muss. Zu jedem $\varepsilon > 0$ ist hier im Sinne von "beliebig klein" zu verstehen. Der grosse Unterschied zur punktweisen Konvergenz ist, dass hier das $\varepsilon$ **für alle** $t\in\mathbb{R}$ das selbe ist. In diesem Sinne ist jedoch die Funktionenfolge $x_n$ nicht mehr konvergent:

```{code-cell} ipython3
:tags: [hide-input]

epsilon = 0.1

plt.plot(t,xn(t,1), label='$xn=1$')
plt.plot(t,xn(t,4), label='$xn=4$')
plt.plot(t,xn(t,8), label='$xn=8$')
plt.plot(t,x(t),'--', label='Grenzfunktion')
plt.fill_between(t,x(t)-epsilon,x(t)+epsilon,label=r'$\varepsilon$-Schranke',alpha=0.3)
plt.legend()
plt.show()
```

Für ein $\varepsilon < 1$ (Sprunghöhe) finden wir kein $n_0$ so, dass der Abstand zur Grenzfunktion der Bedingung genügt. Die Funktionenfolge ist daher nicht konvergent bezüglich der Maximumsmetrik.

```{admonition} Definition: Gleichmässige Konvergenz
Die Konvergenz in $(C[a,b],d_{\max})$ nennt man **gleichmässige** Konvergenz auf dem Intervall $[a,b]$.
```

```{note}
Ist eine stetige Funktionenfolge gleichmässig konvergent, so ist die Grenzfunktion wiederum stetig.
```

### Cauchy-Folge

Der Begriff der Cauchy-Folge folge kann direkt auf metrische Räume übertragen werden:

```{admonition} Definition: Cauchy-Folge
Eine Folge $\|x_n\}$ aus dem metrischen Raum $X$ heisst *Cauchy-Folge*, wenn

$$\lim_{m,n\to\infty} d(x_n,x_m) = 0$$

ist, dh. wenn es zu jedem $\varepsilon > 0$ eine natürliche Zahl $n_0 = n_0(\varepsilon)\in\mathbb{N}$ git mit

$$d(x_n,x_m) < \varepsilon\quad\text{für alle}\quad n,m \ge n_0.$$
```

```{admonition} Satz
Jede konvergente Folge im metrischen Raum $X$ ist auch eine Cauchy-Folge.
```

**Beweis**: Aus der Konvergenz der Folge $\{x_n\}$ folgt: Zu jedem $\varepsilon>0$ gibt es ein $n_0 = n_0(\varepsilon)\in\mathbb{N}$ und ein $x_0\in X$ mit

$$d(x_n,x_0) < \varepsilon\quad\text{und}\quad d(x_m,x_0) < \varepsilon$$

für alle $n,m \ge n_0$. Mit Hilfe der Dreiecksungleichung folgt

$$d(x_n,x_m) \le d(x_n,x_0) + d(x_0,x_m) = d(x_n,x_0) + d(x_m,x_0) < 2 \varepsilon$$

für alle $n,m \ge n_0$. $\Box$

Die Umkehrung gilt im allgemeinen nicht.

**Beispiel**: Betrachte $X = (0,1)$ mit der Metrik $d(x,y) := |x-y|$ und der Folge $\{x_n\}$ mit $x_n = \frac{1}{1+n}$. Die Folge ist eine Cauchy-Folge im metrischen Raum $(X,d)$, besitzt jedoch keinen Grenzwert in diesem $(0\not\in X)$.

Das führt uns zu einem neuen Begriff, der **Vollständigkeit**. Wir interessieren uns insbesondere für diejenigen metrischen Räume, in denen Cauchy-Folgen konvergieren.

```{admonition} Definition: Vollständig
Ein metrischer Raum $X$ heisst *vollständig*, wenn jede Cauchy-Folge in $X$ gegen ein Element von $X$ konvergiert.
```

Betrachten wir ein paar Beispiele:

1. $\mathbb{R}^n$ mit der euklidischen Metrik

   $$d(x,y) = \sqrt{\sum_{k=1}^n |x_k-y_k|^2}$$

   versehen, ist ein vollständiger metrischer Raum. Dies folgt aus dem Cauchyschen Konvergenzkriterium für Punktfolgen.

2. $C[a,b]$ mit der Maximumsmetrik

   $$d(x,y) = \max{a\le t\le b} |x(t)-y(t)|$$

   versehen ist vollständig.

Dagegen ist $C[a,b]$ bezüglich der Integralmetrik

$$d(x,y) = \left(\int_a^b |x(t)-y(t)|^p dt \right)^{1/p},\quad 1\le p < \infty$$ (eq:LpMetrik)

**nicht** vollständig. Wir betrachten dazu für $p=2$ auf $C[0,1]$ folgendes Gegenbeispiel:

Sei $t\in[0,1]$

$$x_n(t) := \begin{cases}
n^\alpha\quad\text{für}\ t \le 1/n\\
\frac{1}{t^\alpha}\quad\text{für}\ t > 1/n\end{cases}$$ (eq:Gegenbeispiel)

und $x(t) = \frac{1}{t^\alpha}$ für $0<\alpha<1/2$.

```{code-cell} ipython3
:tags: [hide-input]

def xn(t,n,alpha):
    y = np.zeros_like(t)
    y[t<=1/n] = n**alpha
    y[t>1/n] = 1/(t[t>1/n]**alpha)
    return y
x = lambda t, alpha: 1/t**alpha
t = np.linspace(0,1,400)
plt.plot(t,xn(t,3,1/3), label='$n=3$')
plt.plot(t,xn(t,10,1/3), label='$n=10$')
plt.plot(t,xn(t,20,1/3), label='$n=20$')
plt.plot(t[1:],x(t[1:],1/3),'--', label='Grenzfunktion')
plt.title(r'$\alpha = 1/3$')
plt.legend()
plt.ylim(0,4)
plt.show()
```

Es gilt $x_n \in C[0,1]$ für alle $n\in\mathbb{N}$ und

$$d(x_n,x)^2 = \int_0^1 |x_n(t)-x(t)|^2 dt = \int_0^{1/n} (n^\alpha - t^{-\alpha})^2 dt.$$

Mit Hilfe der Ungleichung $(a-b)^2 \le 2 (a^2+b^2)$ folgt

$$d(x_n,x)^2 \le 2 \int_0^{1/n} (n^{2\alpha} + t^{-2\alpha})^2 dt = \frac{2}{n^{1-2\alpha}} + \frac{2}{1-2\alpha} \frac{1}{n^{1-2\alpha}} \to 0 \quad \text{für}\ n\to\infty$$

da $1-2\alpha > 0$ gilt. Mit dem  konviergiert die Folge $\{x_n\}$ gegen $x$ in der Integralmetrik ($p=2$). Wir zeigen, dass $\{x_n(t)\}$ keine auf $[0,1]$ stetige Grenzfunktion besitzt. Dazu nehmen wir an, dass $y(t)\in C[0,1]$ sei Grenzfunktion der Folge $\{x_n(t)\}$. Wir setzen

$$M = \max_{0\le t\le 1} |y(t)|.$$

(Muss für eine wohl definierte Obersumme endlich sein!)
Für $t \le (2M)^{-1/\alpha}$ und $n > (2M)^{1/\alpha}$ gilt

$$x_n(t) - y(t) \ge M\quad \text{für}\quad t\le (2M)^{-1/\alpha}$$

und somit

$$d^2(x_n,y) = \int_0^1 |x_n(t)-y(t)|^2 dt \ge \int_0^{(2M)^{-1/\alpha}} |x_n(t)-y(t)|^2 dt \ge (2M)^{-1/\alpha}\cdot M^2 > 0, \quad \text{für}\quad n > (2M)^{1/\alpha}$$

im Widerspruch zur Annahme, dass $\{x_n(t)\}$ gegen $y(t)$ konvergiert. Damit ist die Behauptung bewiesen. $\Box$

Die Tatsache, daß $C[a,b]$, versehen mit einer Integralmetrik, kein vollständiger metrischer Raum ist, bedeutet einen schwerwiegenden Mangel dieses Raumes. Jedoch gibt es mehrere Wege der Vervollständigung:

1. Man erweitert die Klasse $C[a,b]$ zur Klasse der auf $[a,b]$ Lebesgue-integrierbaren Funktionen und interpretiert das in {eq}`eq:LpMetrik` auftretende Integral im Lebesgueschen Sinn. Dadurch gelangt man zum vollständigen metrischen Raum $L_p[a, b]$ (vgl. {cite:p}`heuser_Analysis2`, S. 109).
2. Ein anderer Weg besteht darin, dass der vollständige Erweiterungsraum als Menge von linearen Funktionalen auf einem geeigneten Grundraum nach dem Vorbild der Distributionentheorie aufgefasst wird.
3. Ein weiterer Weg besteht in der abstrakten Konstruktion eines vollständigen Erweiterungsraums mit Hilfe von Cauchy-Folgen. Auf diese Weise lässt sich **jeder** nichtvollständige metrische Raum vervollständigen (vgl. {cite:p}`heuser_FA`, S. 251)

Wir definieren kompakt für metrische Räume wie folgt

```{admonition} Definition: Kompakt
Sein $X$ ein metrischer Raum. $A \subset X$ heisst *kompakt*, wenn jede Folge $\{x_n\}$ aus $A$ eine Teilfolge enthält, die gegen ein Grenzelement $x\in A$ konvergiert.
```

```{admonition} Satz
Jede kompakte Teilmenge $A$ eines metrischen Raumes $X$ ist beschränkt und abgeschlossen.
````

Die Umkehrung gilt im allgemeinen nicht. 

(BestapproximationMetrisch)=
### Bestapproximation in metrischen Räumen

In der Approximationstheorie stellt sich das Grundproblem: In einem metrischen Raum $X$ sei eine Teilmenge $A$ und ein fester Punkt $y\in X$ vorgegeben. Zu bestimmen ist ein Punkt $x_0 \in A$, der von $y$ minimalen Abstand hat. Das Problem beginnt schon damit, dass es nicht klar, ist, dass es einen solchen Punkt überhaupt gibt:

> Betrachte den Raum $(\mathbb{R}, d)$ mit $d(x_1,x_2) = |x_1-x_2|$. Die Teilmenge $A$ sei gegeben durch $A = (0,1)$ und $y=2$. In $A$ gibt es keinen Punkt $x_0$, für den $d(x_0,2)$ minimal ist ($1\not\in A$). 

Zur Erinnerung:

```{admonition} Definition: Supremum, Infimum
Sei $A\subset \mathbb{R}$, dann bezeichnet man mit dem *Supremum* von $A$ die kleinste reelle Zahl $\lambda$ mit $x\le \lambda$ für alle $x\in A$. Analog bezeichnet man mit dem *Infimum* die grösste reelle Zahl $\mu$ mit $x\ge \mu$ für alle $x\in A$.
```

```{admonition} Satz
Es sei $X$ ein metrischer Raum und $A$ eine **kompakte** Teilmenge von $X$. Dann gibt es zu jedem festen Punkt $y \in X$ einen Punkt $x_0 \in A$, der von $y$ kleinsten Abstand hat.
````

Betrachten wir das obige Beispiel angepasst auf die Voraussetzung im Satz: Sei $A = [0,1] \subset \mathbb{R}$ ein kompaktes Intervall, dann ist der Punkt $x_0 = 1$ bestapproximierendes Element.

% Anwendung!

% Macht Probleme beim LaTeX Build...
```{figure} Bestapproximation.png
---
align: left
height: 250px
name: fig-Bestapproximation
---
Bestapproximation
```

## Normierte Räume. Banachräume

Bis jetzt haben wir sehr wenig Eigenschaften eines Raumes benötigt. Was uns noch fehlt, sind abgesehen vom Abstand der Elemente noch *algebraische* Eigenschaften (addieren, multiplizieren, etc.). Dazu definieren wir den *linearen Raum* (oder *Vektorraum*) wie folgt.

```{admonition} Definition: linearer Raum
Ein *linearer Raum* (oder *Vektorraum*) über einem Körper $\mathbb{K}$ besteht aus einer nichtleeren Menge $X$, sowie

* einer Vorschrift, die jedem Paar $(x,y)$ mit $x,y \in X$ genau ein Element $x+y\in X$ zuordnet (*Addition*)
* einer Vorschrift, die jedem Paar $(\lambda,x)$ mit $\lambda\in \mathbb{K}$ und $x \in X$ genau ein element $\lambda x\in X$ zuordnet (*Multiplikation mit Skalaren*), wobei für alle $x,y,z \in X$ und $\lambda, \mu\in\mathbb{K}$ folgende Regeln gelten:

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

* Die Mengen $\mathbb{R}^n$, \mathbb{C}^n$ sind wohl bekannt.
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

```{admonition} Definition: Unterraum, lineare Mannigfaltigkeit, lineare Hülle / Span, linear unabhängig, Dimension, Basis
* Eine nicht leere Teilmenge $S$  von $X$ heisst *Unterraum* von $X$, wenn für bel. $x,y \in S$  und $\lambda \in \mathbb{K}$ stets

  $$x+y \in S\quad\text{und}\quad \lambda x \in S$$

  folgt. Insbesondere ist $S$ selbst ein linearer Raum über $\mathbb{K}$.
* Ist $S$ ein Unterraum von $X$ und $x_0\in X$ beliebig, so nennt man

  $$M = x_0 + S := \{x_0+y\ |\ y\in S\}$$
  eine *lineare Mannigfaltigkeit* von $X$.
* Ist $A$ eine beliebige nichtleere Teilmenge von $X$, so bilden alle Linearkombinationen $\sum_{k=1}^m \lambda_k x_k$ mit beliebigem $m \in \mathbb{N}$, $\lambda_k\in\mathbb{K}$, $x_k \in A$ einen Unterraum $S\subset X$. Er wird *lineare Hülle von* $A$ oder $\mathop{span} A := S$ genannt.

  Man sagt: $A$ spannt $S$ auf oder $A$ ist ein *Erzeugendensystem* von $S$. Im Falle $S=X$ spannt $A$ den ganzen Raum $X$ auf.
* Sind $S$ und $T$ Unterräume von $X$, dann ist die *Summe* $S+T$ definiert durch

  $$S+T := \mathop{span} S \cup T.$$
* Die (endlich vielen) Elemente $x_1, \ldots, x_n\in X$ heissen *linear unabhängig*, wenn aus

  $$\alpha_1 x_1 + \ldots + \alpha_n x_n = 0\quad\text{stets}\quad \alpha_1 = \ldots = \alpha_n = 0$$

  folgt.
* Sei $S$ ein Unterraum von $X$. Wir sagen, $S$ besitzt die *Dimension* $n$, $\mathop{dim} S = n$, wenn es $n$ linear unabhängige Elemente von $S$ gibt, aber $n+1$ Elemente von $S$ stets linear abhängig sind.

  $S$ heisst *Basis* von $X$, wenn die Elemente von $S$ linear unabhängigsind und $\mathop{span} S = X$ gilt.

  Besitzt $X$ keine endlich dimensionale Basis, nennt man $X$ *unendlich-dimeansional* ($\mathop{dim} X = \infty$).
```

**Bemerkungen**: Die oben erwähnten Funktionenräume $C[a,b]$, $C^k[a,b]$, Polynome sind unendlich-dimensional. Ebenso ist der Folgenraum $l_p$ unendlich-dimensional:

> Man betrachte
>
> $$x^{(k)} = \{0, \ldots, 0, 1, 0, \ldots \}\in l_p$$
> mit 1 an der Stelle $k$.

Im folgenden sind wir an Räumen interessiert, für welche eine lineare Struktur und eine Metrik gegeben ist.

```{admonition} Definition: normierter Raum
Sei $X$ ein metrischer und linearer Raum. Zu dem sei die Metrik $d$ von $X$ *translationsinvariant*

$$d(x+z, y+z) = d(x,y)\quad\forall\ x,y,z\in X$$
und *homogen*

$$d(\alpha x, \alpha y) = |\alpha| d(x,y)\quad \forall\ \alpha\in\mathbb{K}, x,y\in X.$$
Dann nennt man $X$ einen *normierten Raum*. Der durch

$$\|x\| := d(x,0)\quad \forall\ x\in X$$
erklärte Ausdruck heisst *Norm* von $x$.
```

**Bemerkungen**: 
* Neben der kurzen Schreibweise $X$, verwendet man häufig auch die Bezeichnung $(X, \|\cdot\|)$. Der Punkt in $\|\cdot\|$ ist als Platzhalter zu verstehen.
* Führt man den normierten Raum $X$ mit Hilfe einer Norm ein, so ist durch

  $$d(x,y) = \|x-y\|\quad\text{für alle}\ x,y\in X$$
  eine Metrik in $X$ gegeben.


```{admonition} Folgerung
Ein normierter Raum $(X, \|\cdot\|)$ ist ein linearer Raum, auf dem eine Norm $\|\cdot\|$ erklärt ist, die für alle $x,y\in X$ und $\alpha \in \mathbb{K}$

$$\begin{split}
\|x\| & \ge 0,\quad \|x\|=0\quad \text{genau dann, wenn $x=0$ ist},\\
\|\alpha x\| & = |\alpha| \|x\|\\
\|x + y\| & \le \|x\| + \|y\|\end{split}$$
erfüllt.
```

Damit können wir einen wichtigen Begriff der Funktionalanalysis einführen, den Banachraum:

```{admonition} Definition: Banachraum
*Vollständig normierte* Räume $X$ sind diejenigen, für die jede Cauchy-Folge in $X$ gegen ein Element in $X$ konvergiert.

Ein vollständiger normierter Raum heisst *Banachraum*.
```

**Beispiele**: Folgende Räume sind Banachräume

* $\mathbb{R}^n$ mit der Norm $\|x\| = \sqrt{\sum_{k=1}^n |x_k|^2}$.
* $C[a,b]$ mit der Norm $\|x\| := \max_{a\le t \le b} |x(t)|$.
* $C^k[a,b]$ mit der Norm $\|x\| := \max_{a\le t \le b} |x(t)| + \max_{a\le t \le b} |x'(t)| + \ldots + \max_{a\le t \le b} |x^{(k)}(t)|.$  
* $l_p\ (1\le p < \infty)$ mit der Norm $\|x\| = \left(\sum_{k=1}^{\infty} |x_k|^p \right)^{1/p}$

```{admonition} Definition: abzählbare Basis
Man sagt, dass $X$ eine *abzählbare Basis* $\{x_k\}_{k=1}^{\infty}$ mit $x_k \in X$ besitzt, falls jedes $x\in X$ eindeutig in der Form $x = \sum_{k=1}^{\infty} \alpha_k\,x_k$ darstellbar ist, wobei die Konvergenz bezüglich der Norm von $X$ zu verstehen ist.
````

Lineare Räume können durchaus verschieden normiert werden. Als Beispiel betrachte dazu den Raum $X=\mathbb{R}^n$ mit den Normen

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

```{admonition} Satz
Alle Normen in einem **endlich**-dimensionalen Raum $X$ sind äquivalent.
```

Dieses Resultat gilt für unendlich-dimensionale Räume **nicht**. Als Beispiel sei der lineare Raum $C[a,b]$ mit der Maximumsnorm und der Quadratnorm erwähnt. Die beiden Normen sind nicht äquivalent, vgl. das Gegenbeispiel {eq}`eq:Gegenbeispiel`. 

Mit diesem Satz folgt

```{admonition} Satz
Jeder endlich-dimensionale normierte Raum $X$ ist vollständig, also ein Banachraum.
```

## Skalarprodukträume. Hilberträume

Das aus der linearen Algebra bekannte Skalarprodukt lässt sich auch auf unendlich-dimensionale Räume übertragen. Wir definieren ganz allgemein

```{admonition} Definition: Skalarprodukt, Skalarproduktraum
Unter einem *Skalarproduktraum* versteht man einen linearen Raum $X$ über $\mathbb{K}$, in dem ein *Skalarprodukt* $(x,y)$ mit folgenden Eigenschaften definiert ist: Für beliebige $x,y,z \in X$ und $\alpha\in\mathbb{K}$ ist

$$(\cdot, \cdot) : X \times X \to \mathbb{K}$$

und es gilt

$$\begin{split}
(x,x) & \ge 0, \quad (x,x) = 0\ \Leftrightarrow\ x=0\\
(x,y) & = \overline{(y,x)}\\
(\alpha x, y) & = \alpha (x,y)\\
(x+y,z) & = (x,z) + (y,z).
\end{split}$$ (eq:eigenschaftenskalarprodukt)
```

Beispielsweise lässt sich auf $X=C[a,b]$ Menge der reellwertigen stetigen Funktionen auf dem Intervall $[a,b]$ durch

$$(x,y) := \int_a^b x(t)\,y(t) dt\quad x,y \in X$$ (eq:L2SkalarProdukt)

ein Skalarprodukt definieren.

```{admonition} Aufgabe
Weise die Eigenschaften {eq}`eq:eigenschaftenskalarprodukt` für das Skalarprodukt {eq}`eq:L2SkalarProdukt` nach.
```

```{admonition} Satz
Es sei $X$ ein Skalarproduktraum. Für beliebige $x,y,z\in X$ und $\alpha \in \mathbb{K}$ gilt

$$\begin{split}(x, \alpha y) & = \overline{\alpha} (x,y)\\
(x, y+z) & = (x,y) + (x,z)\end{split}$$
```

**Beweis**: selber durchführen.

```{admonition} Satz: induzierte Norm
In jedem Skalarproduktraum $X$ lässt sich durch

$$\|x\| := \sqrt{(x,x)}$$

eine Norm definieren. Man bezeichnet sie als *die durch das Skalarprodukt* $(x,y)$ *induzierte Norm*.
```

Der Beweis lässt sich einfach mit Hilfe der *Schwarz'schen Ungleichung* durchführen.

```{admonition} Satz: Schwarzsche Ungleichung
Es sei $X$ ein Skalarproduktraum. Dann gilt für alle $x,y\in X$

$$|(x,y)| \le \|x\|\,\|y\|.$$
```

Die Beweise sind dem Leser überlassen (vgl. {cite:p}`burg_wille_haf_meister_PDE` S. 43, 44).

```{admonition} Satz
Sei $X$ ein Skalarproduktraum. Ferner seien $\{x_n\}$ und $\{y_n\}$ Folgen aus $X$ mit $x_n \to x$ und $y_n \to y$ für $n\to\infty$, wobei die Konvergenz im Sinne der induzierten Norm zu verstehen ist. Dann gilt

$$(x_n, y_n) \to (x,y)\quad\text{für}\quad n\to\infty,$$

dh. das Skalarprodukt ist eine stetige Funktion bezüglich der Normkonvergenz.
```

```{admonition} Definition: Hilbertraum
Ein Skalarproduktraum $X$, der bezüglich der durch das Skalarprodukt induzierten Norm

$$\|x\| = \sqrt{(x,x)}\quad\text{für}\ x\in X$$
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

```{admonition} Definition: orthogonal
Es sei $X$ ein Hilbertraum.

* Zwei Elemente $x,y \in X$ heissen *orthogonal* ($x\perp y$), wenn

  $$(x,y) = 0$$
  
  gilt.
* Zwei Teilmengen $A,B \subset X$ heissen *orthogonal* ($A\perp B$), wenn
  
  $$(x,y) = 0 \quad \forall\ x\in A, y\in B$$

  gilt.
* Sei $M\subset X$, dann heisst

  $$M^\perp := \{y\in X | (x,y)=0\quad\forall\,x\in M\}$$

  *Orthogonalraum* von $M$.
* Sei $X'$ ein abgeschlossener Unterraum von $X$. $X''$ wird *orthogonales Komplement* von $X'$ genannt, wenn

  $$X''\perp X' \quad \text{und}\quad X' \oplus X'' = X$$
  gilt. Mit $\oplus$ ist die direkte Summe bezeichnet.
```

Es gilt

```{admonition} Satz: Pythagoras
Es sei $X$ ein Hilbertraum und seien $x,y \in X$ mit $x\perp y$. Dann gilt

$$\|x+y\|^2 = \|x\|^2 + \|y\|^2.$$
```

**Beweis**: Einfaches Nachrechnen.

```{admonition} Satz
Es sei $X$ ein Hilbertraum und $M$ eine beliebige Teilmenge von $X$. Dann ist der Orthogonalraum $M^\perp$ von $M$ ein abgeschlossener Unterraum von $X$.
```

Für den Beweis sei auf {cite:p}`burg_wille_haf_meister_PDE` S. 50 verwiesen.

Wir kommen nun auf das Approximationsproblem aus {ref}`BestapproximationMetrisch` zurück und können die Struktureigenschaften des Hilbertraumes nutzen:

```{admonition} Satz
Es sei $X'$ ein abgeschlossener Unterraum von $X$ und $x\in X$ beliebig. 
* Dann existiert genau ein $x_0 \in X'$ mit

  $$\|x-x_0\| = \min_{x'\in X'} \|x-x'\|,$$

  dh. zu jedem $x\in X$ gibt es genau ein bestapproximierendes Element bezüglich $X'$.
* Es gilt $x_0\in X$ ist genau dann bestapproximierend an $x\in X$, wenn

  $$(x-x_0,y) = 0\quad \forall y\in X'$$
  gilt.
```

Für den Fall, dass $X'$ endlich-dimensional ist, lässt sich das bestapproximierende Element konstruieren. Es gilt

```{admonition} Satz
Sei $X'\subset X$ mit $\dim X' < \infty$ ein Unterraum des Hilbertraumes $X$ und sei $x_1, \ldots, x_n$ eine Basis von $X'$. Dann lässt sich das eindeutig bestimmte bestapproximierende Element $x_0\in X'$ an $x\in X$ in der Form

$$x_0 = \sum_{k=1}^n \lambda_k x_k$$

darstellen, wobei die Koeffizienten $\lambda_1, \ldots, \lambda_n$ durch das lineare Gleichungssystem

$$(x-x_0,x_i) = (x,x_i) - \sum_{k=1}^n \lambda_k (x_k,x_i) = 0 \quad i = 1, \ldots, n$$ (eq:bestapproxHilbert)

gegeben sind.
```

**Bemerkung**: Bildet $x_k$, $k=1,\ldots, n$ ein *Orthonormalsystem*

$$(x_i,x_k) = \delta_{i k} = \begin{cases}
0\quad \text{für}\ i\not= k\\
1\quad \text{für}\ i = k\end{cases},$$
so folgt aus {eq}`eq:bestapproxHilbert` sofort

$$\lambda_i = (x,x_i)\quad i = 1, \ldots,n.$$

Die Koeffzienten $\lambda_i$ nennt man auch Fourierkoeffizienten!

> Wir werden zeigen, dass sich mit Hilfe des Schmidtschen Orthogonalisierungsverfahrens aus $n$ linear unabhängigen Elementen stets ein Orthogonalsystem konstruieren lässt.

```{admonition} Definition: Orthonormalsystem (ONS)
Es sei $X$ ein Hilbertraum. Man nennt die Folge $\{x_k\}_{k\in\mathbb{N}}$ ein (abzählbares) *Orthonormalsystem* (kurz ONS) von $X$, wenn 

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

```{admonition} Satz: Orthogonalisierungsverfahren nach Erhard Schmidt
Gegeben sei eine abzählbare (nicht zwingend endliche) linear unabhängige Folge $\{y_k\} \subset X$ aus dem Hilbertraum $X$. Dann gibt es ein ONS aus $n$ bzw. abzählbar unendlich viele Elementen $\{x_k\}$ so, dass der von der Folge $\{y_k\}$ aufgespannte (abgeschlossene) Unterraum mit dem der Folge $\{x_k\}$ aufgespannten (abgeschlossene) Unterraum übereinstimmt.
```

Der Beweis beruht auf der Konstruktion: sei

$$x_1 = \frac{y_1}{\|y_1\|}$$

so ist $\mathop{span}\{y_1\} = \mathop{span}\{x_1\}$. Wir nehmen nun an, es seien bereits $k$ orthonormierte Elemente $x_1, \ldots, x_k$ mit $\mathop{span}\{y_1, \ldots, y_k\} = \mathop{span}\{x_1,\ldots, x_k\}$ konstruiert. Dann setze

$$z_{k+1} = y_{k+1} - \sum_{j=1}^k (y_{k+1},x_j) x_j.$$

In dem Fall gilt $(z_{k+1},x_i) = 0$ für alle $i = 1, \ldots, k$. Mit

$$x_{k+1} = \frac{z_{k+1}}{\|z_{k+1}\|}$$

folgt $\mathop{span}\{y_1, \ldots, y_{k+1}\} = \mathop{span}\{x_1,\ldots, x_{k+1}\}$.

**Bemerkung**: Das Verfahren wird auch in der Numerik angewandt.

:::{seealso}
[Beispiel für das Orthogonalisierungsverfahren.](BeispielOrthogonalisierungSchmidt.ipynb)
:::

Die Elemente eines Hilbertraumes können mit Hilfe eines vollständigen Orthogonalsystems dargestellt werden. Dies gelingt mit Hilfe der _verallgemeinerten Fourierreihen_:

```{admonition} Satz: Fourierentwicklung

Sei $X$ ein Hilbertraum mit einem vollständigen ONS $\{x_k\}_{k\in\mathbb{N}}$.
* Dann lässt sich jedes $x\in X$ in der Summenform

  $$x = \sum_{k=1}^\infty a_k\,x_k\quad\text{Fourierentwicklung von $x$}$$

  mit eindeutig bestimmten Koeffizienten

  $$a_k = (x,x_k)\in\mathbb{K}$$

  darstellen und die Reihe $\sum_{k=1}^\infty |a_k|^2$ ist konvergent.
* Umgekehrt gibt es zu jeder Zahlenfolge $\{a_k\}_{k\in\mathbb{N}}$ in $\mathbb{K}$, für die $\sum_{k=1}^\infty |a_k|^2$ konvergiert, genau ein $x\in X$ mit $x=\sum_{k=1}^\infty a_k x_k$.
```

**Bemerkung**: Aufgrund der Darstellung $x = \sum_{k=1}^\infty a_k\,x_k$ nennt man ein vollständiges ONS auch eine *Hilbertraumbasis*.

```{admonition} Satz: Struktur von Hilberträumen
Es sei $X$ ein Hilbertraum und $\{x_k\}_{k\in\mathbb{N}}$ ein (abzählbares) ONS in $X$. Dann sind die folgenden Aussagen äquivalent:
* $X = \overline{\underset{k\in\mathbb{N}}{\bigoplus} \mathop{span}(x_k)}$.
* Das ONS $\{x_k\}_{k\in\mathbb{N}}$ ist abgeschlossen.
* Für alle $x\in X$ gilt die *Parsevalsche Gleichung*

  $$\sum_{k=1}^\infty |(x,x_k)|^2 = \|x\|^2\quad\text{(Vollständigkeitsrelation)}$$

* Jedes Element $x\in X$ besitzt die Fourierentwicklung

  $$x = \sum_{k=1}^\infty (x,x_k)\,x_k.$$
```

:::{seealso}
[Beispiel zur Fourierentwicklung.](BeispielFourierEntwicklung.ipynb)
:::
