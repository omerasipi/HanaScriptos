# Variationsrechnung

## Einführendes Beispiel

Wir beginnen mit einem eindimensionalen Problem: gesucht ist $u\in C^1(a,b)$ so, dass

$$J(u) = \int_a^b \left(\frac{1}{2} u'(x)^2 - f(x) u(x)\right) dx \to \min$$

gilt. 

Wie in der Analysis reeller Funktionen soll für das Beispiel ein notwendiges Kriterium für eine Extremalstelle des Funktionals berechnet werden. Für die Herleitung gehen wir davon aus, dass $u(x)$ ein Minimierer des Funktionals $J$ ist. Diesen stören wir mit einer beliebigen Funktion $v(x) \in C_0^\infty(a,b)$. Das bedeutet, dass $v$ beliebig oft stetig differenzierbar ist und $v(a) = v(b) = 0$ gilt. Wir betrachten daher $J(u(x) + \varepsilon v(x))$. Falls $u(x)$ eine Extremalstelle ist, so muss

$$\frac{d}{d\varepsilon} J(u(x) + \varepsilon v(x))\big|_{\epsilon = 0} = 0$$

gelten.

```{figure} GateauxAbleitung.png
---
align: left
height: 250px
name: fig-GateauxAbleitung
---
Störung der Funktion $u(x)$ durch $\varepsilon v(x)$
```

Wir berechnen daher

$$\begin{split}
\frac{d}{d\varepsilon} J(u(x) + \varepsilon v(x))\big|_{\epsilon = 0} & = \frac{d}{d\varepsilon} \int_a^b \left( \frac{1}{2} (u'(x) + \varepsilon v'(x))^2 - f(x) (u(x) + \varepsilon v(x)) \right) dx \Big|_{\epsilon = 0}\\
& = \int_a^b \left((u'(x) + \varepsilon v'(x)) \cdot v'(x) - f(x) v(x)\right) dx  \Big|_{\epsilon = 0}\\
& = \int_a^b \left(u'(x) \cdot v'(x) - f(x) v(x) \right) dx
\end{split}$$

Es folgt die **schwache Gleichung**

$$ \int_a^b \left(u'(x) \cdot v'(x) - f(x) v(x) \right) dx = 0\quad \forall v \in C_0^\infty(a,b).$$

Unter der Voraussetzung, dass $u \in C^2(a,b)$ gilt, folgt mit Hilfe der partiellen Integration

$$ \underbrace{\left[ u'(x) v(x)\right]_a^b}_{=0\ \text{da}\ v(a)=v(b)=0} + \int_a^b \left(- u''(x) v(x) - f(x) v(x) \right) dx = 0\quad \forall v \in C_0^\infty(a,b).$$

und somit

$$ \int_a^b \left(- u''(x) - f(x)\right) v(x) dx = 0\quad \forall v \in C_0^\infty(a,b).$$

Der **Fundamentalsatz der Variationsrechnung** besagt, dass die Gleichung nur dann erfüllt sein kann, wenn die **starke Gleichung**

$$- u''(x) - f(x) = 0\quad\text{oder}\quad - u''(x) = f(x)$$

gilt.

Wie aus dem Beispiel zu sehen ist, sind die Ansprüche an die Regularität in der starken Gleichung in der Regel grösser als die der schwachen Gleichung. Für skalare Funktionen ist das in der Regel kein grösseres Problem. Hingegen im Mehrdimensionalen führt dies zur Regularitätstheorie. Es muss studiert werden, wann eine schwache Lösung auch eine starke ist. Dies hängt im Mehrdimensionalen mit unter auch vom Gebiet $\Omega \subset \mathbb{R}^n$ ab.

Im Beispiel sind ein paar Grundkonzepte enthalten, welche wir wie folgt festhalten

```{admonition} Definition: Gâteaux-differenzierbar
Sei $X$ ein Banachraum. Existiert für $J: X \to \mathbb{R}$ ein gegebenes Funktional

$$G'(\varepsilon=0) = \lim_{\varepsilon \to 0} \frac{J(u + \varepsilon v) - J(u)}{\varepsilon},$$

so ist $J$ in $u$ in Richtung $v$ Gâteaux-differenzierbar und die Ableitung wird als $dJ(u,v)$ bezeichnet.
```

Das Gâteaux-Differential muss nicht zwingend linear noch stetig in $v$ sein. 

```{admonition} Definition: Fréchet-Ableitung
Existiert $dJ(u,v)$ in $u\in D \subset X$ für $v\in X$ und ist $dJ(u,v)$ linear in $v$, so heisst $dJ(u,v)$ die **erste Variation** oder **Fréchet-Ableitung** von $J$ in Richtung $v$. Wir schreiben in dem Fall 

$$dJ(u,v) = \delta J(u) v = J'(u) v.$$

Ist dies für alle $v\in D\subset X$ richtig, wobei $D$ ein Unterraum von $X$ ist, so ist

$$J'(u) : D \to \mathbb{R}$$

ein lineares Funktional.
```

Betrachten wir obiges Beispiel noch im mehrdimensionalen: gesucht ist $u\in H_0^1(\Omega)$ so, dass

$$J(u) = \int_\Omega \left(\frac{1}{2} |\nabla u|^2 - f(x) u(x)\right) dx \to \min$$

gilt. 

Wir berechnen die Gâteaux-Ableitung

$$\frac{d}{d\varepsilon} J(u + \varepsilon v)\big|_{\varepsilon = 0} = \int_\Omega \nabla u\cdot \nabla v - f v dx$$

Die Gâteaux-Ableitung ist linear in $v$ und somit die erste Variation des Funktionals $J : H_0^1(\Omega) \to \mathbb{R} $.

Mit der Forderung für einen Extremalpunkt $J'(u) v = 0$ für alle $v\in H_0^1(\Omega)$ haben wir die schwache Gleichung (vgl. Einleitung Poisson Gleichung {eq}`eq:weakPoisson`)

$$\int_\Omega \nabla u\cdot\nabla v dx = \int f(x) v dx \quad \forall v \in H_0^1(\Omega).$$

Mit Hilfe des Gauss-Theorem oder partiellen Integration im mehrdimensionalen folgt die starke Gleichung (partielle Differentialgleichung)

$$-\Delta u = f(x)\quad x\in\Omega$$

mit der Dirichlet Randbedingung $u(x) = 0$ auf dem Rand $\partial\Omega$.

```{note}
Halten wir fest: Durch *Testen* der partiellen Differentialgleichung erhalten wir die schwache Gleichung, welche dem notwendigen Kriterium (Fréchet-Ableitung ist Null) eines zu (typischerweise) minimierenden Funktionals entspricht.
```

## Euler-Lagrange Gleichung

Für die Herleitung der formalen Euler-Lagrange Gleichung beschränkten wir uns auf eindimensionale Probleme. Wir betrachten das zu minimierende Funktional

$$J(u) = \int_a^b F(x,u,u') dx,$$ (eq:EulerLagrageFunktional)

welches auf $C^1(a,b)$ definiert ist.

```{admonition} Definition: lokaler Minimierer
Sei $X$ ein Banachraum. Eine Funktion $u\in X$ heisst lokaler Minimierere für das Funktional $J : X \to \mathbb{R}$, falls

$$J(u) \le J(v)\quad\text{für alle}\ v\in X\ \text{mit}\ \|u-v\|_X < \delta$$

mit einer Konstante $\delta > 0$ gilt.
```

Analog wird ein lokaler Maximierer definiert. Man beschränkt sich in der Regel auf Minimierer und betrachtet für Maximierer das Funktional $-J$.

Wir betrachten nun das Funktional {eq}`eq:EulerLagrageFunktional` und berechnen die Gâteaux-Ableitung für beliebiges $v\in C_0^1(a,b)$

$$\begin{split}
\frac{d}{d\varepsilon} J(u+\varepsilon v)\big|_{\varepsilon=0} & = \frac{d}{d\varepsilon} \int_a^b F(x,u+\varepsilon v,u'+\varepsilon v') dx \Big|_{\varepsilon=0}\\
& = \int_a^b \left(\partial_{u} F(x,u,u')\cdot v + \partial_{u'} F(x,u,u')\cdot v' \right) dx
\end{split}$$

Partielle Integration des zweiten Summanden liefert

$$J'(u) v = \int_a^b \left(\partial_{u} F(x,u,u')\cdot v - \frac{d}{dx} \partial_{u'} F(x,u,u')\cdot v\right) dx + \left[ \partial_{u'} F(x,u,u')\cdot v \right]_a^b.$$

Da $v\in C_0^1(a,b)$ folgt $v(a) = v(b) = 0$ und damit

$$J'(u) v = \int_a^b \left(\partial_{u} F(x,u,u') - \frac{d}{dx} \partial_{u'} F(x,u,u')\right)\cdot v dx.$$

Wir können damit folgenden Satz festhalten:

```{admonition} Satz
Die Funktion $u\in C^1(a,b)$ sei ein lokaler Minimierer für das Funktional {eq}`eq:EulerLagrageFunktional` und die Lagrange-Funktion $F: [a,b] \times \mathbb{R} \times \mathbb{R}$ sei stetig und bezüglich der letzten beiden Variablen stetig partiell differenzierbar. Dann gilt die **Euler-Lagrange-Gleichung**

$$\frac{d}{dx} \partial_{u'} F(x,u,u') =  \partial_{u} F(x,u,u').$$ (eq:StarkeEulerLagrangeGleichung)
```

**Bemerkungen**:
* Wir haben den Satz für stetige Funktionen auf ganz $[a,b]$ notiert. Das ganze lässt sich leicht auf stückweise stetige Funktionen verallgemeinern.

* Im hier notierten Fall von Funktionen in *einer* Variablen $x$ gilt, dass eine Lösung der schwachen Gleichung

  $$\int_a^b \left(\partial_{u} F(x,u,u')\cdot v + \partial_{u'} F(x,u,u')\cdot v' \right) dx = 0\quad \forall\ v\in C_0^1(a,b)$$

  auch Lösung der starken Gleichung {eq}`eq:StarkeEulerLagrangeGleichung` ist.

## Beispiele zur Lösung der Euler-Lagrange-Gleichung

### Kürzeste Verbindung zweier Punkte

Wir betrachten in der Ebene zwei Punkte $(x_1,y_1), (x_2,y_2)$ und suchen nach der kürzesten Kurve, welche diese verbindet. Beschreibt man die Kurve mit einer Funktion $y: [x_1,x_2] \to \mathbb{R}$, so können wir die Länge mit dem Funktional

$$J(y) = \int_{x_1}^{x_2} \sqrt{1+y'(x)^2} dx$$

berechnen. Gesucht ist daher ein Minimierer des Funktionals $J(y)$. Die zugehörige Euler-Lagrange Gleichung lautet

$$\frac{d}{dx} \frac{y'(x)}{\sqrt{1+y'(x)^2}} = 0,$$

wobei ein Faktor $2$ schon gekürzt wurde. Integration liefert

$$\frac{y'(x)}{\sqrt{1+y'(x)^2}} = a.$$

Durch auflösen nach $y'$ und integrieren erhalten wir den Kanditaten für einen Minimierer

$$y(x) = \frac{a}{\sqrt{1-a^2}} x + b = m\, x + b,$$

also eine Gerade zwischen den beiden Punkten. Die beiden Integrationskonstanten können aus den Randbedingungen $y(x_1) = y_1$ und $y(x_2) = y_2$ bestimmt werden. Es stellt sich natürlich die Frage, ist die Lösung ein Minimierer und falls ja ein globaler?

Analog zur Analysis reeller Funktionen kann auf für Funktionale unter geeigneten Voraussetzungen ein Kriterium, die zweite Variation, berechnet werden. Es gilt der Satz

```{admonition} Satz
Sei die Lagrange-Funktion {eq}`eq:EulerLagrageFunktional` zweimal bezüglich den letzten beiden Variablen stetig partiell differenzierbar. Dann ist die zweite Variateion gegeben durch die Bilinearform

$$\delta^2 J(u)(v,v) = \int_a^b \partial_{u,u}F(x,u,u') v^2 + 2 \partial_{u,u'} F(x,u,u') v v' + \partial_{u',u'} F(x,u,u') v'^2 dx.$$

Für einen lokalen Minimierer  gilt die **notwendige Bedingung von Legendre**

$$\partial_{u',u'} F(x,u,u') \ge 0\quad \forall\ x\in [a,b].$$

Ist

