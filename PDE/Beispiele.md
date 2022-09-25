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

# Was ist eine partielle Differentialgleichung

Es sei $D$ ein Gebiet in $\mathbb{R}^n$ und $x = (x_1,...,x_n) \in D$. Unter einer partiellen Differentialgleichung der Ordnung $k$ für eine Funktion $u(x)$ in $D$ versteht man eine Gleichung der Form

$$F\left(x, u(x), \frac{\partial u}{\partial x_1}, \ldots,  \frac{\partial u}{\partial x_n}, \ldots ,  \frac{\partial^k u}{\partial x_n^k}\right) = 0.$$ (eq:defpde)

In $u$ treten also mehrere unabhängige Veränderliche, nämlich $x_1,\ldots , x_n$ auf, und in {eq}`eq:defpde` neben $x$ und $u$ partielle Ableitungen von $u$ bis zur Ordnung $k$. Im Spezialfall $n = 1$ liegt mit {eq}`eq:defpde` eine gewöhnliche Differentialgleichung vor. Die partielle Differentialgleichung heisst *linear*, wenn die Funktion $u(x)$ und ihre Ableitungen in Form einer Linearkombination in der Gleichung {eq}`eq:defpde` auftreten.

Um zu eindeutig bestimmten Lösungen von {eq}`eq:defpde` zu gelangen, sind zusätzliche Bedingungen zu stellen, etwa *Rand- und/oder Anfangsbedingungen* oder Abklingbedingungen im Unendlichen.

Wir betrachten hier ein paar ausgewählte Beispiele von (in der Physik auftretenden) partiellen Differentialgleichungen.

* **Transport**: 1. Ordnung, linear

  $$u_x + u_y = 0$$

* **Transport**: 1. Ordnung, linear

  $$u_x + y\, u_y = 0$$

* **Stosswelle**: 1. Ordnung, nichtlinear

  $$u_x + u\, u_y = 0$$

* **Laplace-Gleichung**: oder Potentialgleichung, 2. Ordnung, linear

  $$u_{xx} + u_{yy} = 0\quad\text{oder kompakt}\ \Delta u(x) = 0$$

* **Poisson-Gleichung**: 2. Ordnung, linear

  $$-u_{xx} - u_{yy} = f(x,y)\quad\text{oder kompakt}\ -\Delta u(x) = f(x)$$

* **Wärmeleitungsgleichung**: 2. Ordnung, linear, eindimensional im Ort $D\subset \mathop{R}$

  $$u_t-u_{xx} = q(x)$$

  oder mehrdimensional im Ort $D\subset\mathop{R}^n$

  $$u_t-\Delta u(x) = q(x)$$

* **Wellengleichung**: 2. Ordnung, linear, eindimensional im Ort $D\subset \mathop{R}$

  $$u_{tt}-u_{xx} = 0$$

  oder mehrdimensional im Ort $D\subset\mathop{R}^n$

  $$u_{tt}-\Delta u(x) = 0$$

* **Helmholtzsche Schwingungsgleichung**: 2. Ordnung, linear

  $$u_{xx}+k^2 u(x) = 0, \quad k\in\mathbb{C}$$

  oder mehrdimensional im Ort $D\subset\mathop{R}^n$

  $$\Delta u(x)+k^2 u(x) = 0, \quad k\in\mathbb{C}$$

* **Schwingender Stab**: 4. Ordnung, linear, eindimensional im Ort $D\subset \mathop{R}$

  $$u_{tt}-u_{xxxx} = 0$$

  oder mehrdimensional im Ort $D\subset\mathop{R}^n$

  $$u_{tt}-\Delta^2 u(x) = 0(x)$$

  $\Delta^2u(x) = \Delta\Delta u(x)$ nennt man den Biharmonischen Operator.

* **Maxwell Gleichungen**: System 1. Ordnung $D \subset \mathop{R}^3$

  $$\begin{array}{rll}
	\frac{\partial B}{\partial t} + \mathop{rot} E & = 0 & \quad\text{Faraday's Gesetz}\\
	\frac{\partial D}{\partial t} + J & = \mathop{rot} H & \quad\text{Maxwell-Ampère Gesetz}\\
	\mathop{div} D & = \rho  & \quad\text{Gauss'sches Gesetz}\\
	\mathop{div} B & = 0 & \quad\text{magnetisches Gauss'sches Gesetz}\\
	B = \mu\, H,\ D & = \epsilon\, E & \quad \text{Material Gesetz}
  \end{array}$$

* **Navier-Stokes Gleichungen**: System 2. Ordnung $D \subset \mathop{R}^3$

  $$\begin{split}
	-\Delta u + Re\,(u\, \nabla)u - \nabla p & = f\\
	\mathop{div} u & = 0
  \end{split}$$

Die Liste ist natürlich bei weitem nicht vollständig! Es zeigt sich, dass die Theorie der partiellen Differentialgleichungen ein recht umfangreiches mathematisches Gebiet darstellt, bei dem sehr unterschiedliche Verfahren und Methoden Verwendung finden und das auch unter Gesichtspunkten der mathematischen Forschung Aktualität besitzt. Dies macht die Entwicklung in den letzten Jahrzehnten deutlich.