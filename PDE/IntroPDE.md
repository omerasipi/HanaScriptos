# Einführung partielle Differentialgleichungen

In der Technik und  Naturwissenschaft haben wir es oft mit partiellen Differentialgleichungen (PDE) zu tun: So kann das mechanische Schwingungsverhalten, die Ausbreitung von Licht, Strömungen um Objekte wie Autos, Flugzeuge, etc., elektromagnetische Felder von Transformatoren und anderen Bauteilen berechnet und für die technische Anwendung optimiert werden. Insbesondere die numerischen Methoden für PDE sind in der Ingenieurwissenschaft nicht mehr weg zu denken. In der Entwicklung werden oft ganze Systeme mathematisch modelliert, numerisch simuliert und optimiert, bevor der erste Prototyp erstellt wird. In den meisten Fällen baut die Modellierung auf einer Formulierung mit Hilfe von PDE auf. 

Die Theorie der PDE ist ein sehr umfangreiches mathematisches Gebiet, bei dem sehr unterschiedliche Verfahren und Methoden Anwendung finden. Wir werden uns in diesem Kurs insbesondere auf eine in der Anwendung sehr weit verbreitete numerische Methoden für PDE konzentrieren, die Methode der finiten Elemente.


## Was ist eine partielle Differentialgleichung

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

## Klassifikation 

Wie wir in den Beispielen gesehen haben, führen wichtige in den Anwendungen auftretende Probleme auf lineare partielle Differentialgleichungen 2-ter Ordnung. 

Die generelle Form einer linearen PDE 2. Ordnung definieren wir wie folgt: finde $u: \Omega \subset \mathbb{R}^d \to \mathbb{R}$ so, dass

$$- \sum_{i,j=1}^d \frac{\partial}{\partial x_i}\left(a_{i,j}(x) \frac{\partial u(x)}{\partial x_j}\right) + \sum_{i=1}^d b_i(x) \frac{\partial u(x)}{\partial x_i} + c(x) u(x) = f(x)$$ (eq:generellPDE)

Die Koeffizienten $a_{i,j}(x), b_i(x), c(x)$ und die rechte Seite $f(x)$ sind gegebene Funktionen. Zusätzlich sind verschiedene Typen von Randbedingungen notwendig. Das Verhalten der PDE hängt im wesentlichen vom Typ des Differential Operators

$$L := \sum_{i,j=1}^d \frac{\partial}{\partial x_i}\left(a_{i,j} \frac{\partial }{\partial x_j}\right) + \sum_{i=1}^d b_i \frac{\partial}{\partial x_i} + c $$

ab. Ersetzt man $\frac{\partial}{\partial x_i}$ mit $s_i$, dann beschreibt

$$\sum_{i,j=1}^d s_i a_{i,j} s_j     + \sum_{i=1}^d b_i s_i + c = 0$$

eine quadratische Form im $\mathbb{R}^d$. Wir unterscheiden folgende Fälle:

* Im Fall, dass $a = (a_{i,j})$ eine (positive oder negative) definite Matrix ist, ist die Form eine Ellipse. Die entsprechende PDE heisst in dem Fall **elliptisch**. Als einfaches Beispiel betrachte $a = \mathbb{1}, b=0$ und $c=0$. In dem Fall folgt die Poisson Gleichung

  $$- \sum_{i=1}^d \frac{\partial^2u(x)}{\partial x_i^2} = f.$$

  Elliptische PDE benötigen Randbedingungen.

  ```{seealso}
  Beispiel elliptische Form {ref}`chap:ellipticexmp`
  ```

* Ist die Matrix $a$ semi-definite und $b^T\cdot \xi \not= 0 $ für $\xi\in \mathop{kern}(a)$, dann beschreibt die Form eine Parabel. In dem Fall sprechen wir von einer **parabolischen** PDE. Ein Beispiel dazu ist

  $$-\sum_{i=1}^{d-1} \frac{\partial^2u}{\partial x_i^2} + \frac{\partial u}{\partial x_d} = f$$

  Oft entspricht die separate Richtung der Zeit und wird typischer weise wie folgt geschrieben
  
  $$\dot{u}(t,x) - \Delta u(t,x) = f,$$

  wobei hier die Reihenfolge der Summanden vertauscht wurde.
  Dieser Typ von Gleichungen benötigt Randbedingungen auf dem $d-1$ dimensionalen Rand und Anfangswerte für die andere Richtung.

  ```{seealso}
  Beispiel parabolische Form {ref}`chap:parabolexmp`
  ```

* Hat die Matrix $a$ $d-1$ positive und einen negativen (oder visa versa) Eigenwerte, dann beschreibt die Form ein Hyperboloid. Man nennt die PDE **hyperbolisch**. Das einfachste Beispiel ist gegeben durch

  $$-\sum_{i=1}^{d-1} \frac{\partial^2u}{\partial x_i^2} + \frac{\partial^2 u}{\partial x_d^2} = f$$

  Auch hier ist die andere Richtung typischerweise die physikalische Zeit:

  $$\ddot{u}(t,x) - \Delta u(t,x) = f.$$ (eq:exmpparabolclass)

  In dem Fall werden zwei Anfangswerte benötigt.

  ```{seealso}
  Beispiel hyperbolische Form {ref}`chap:hyperexmp`
  ```

* Falls die Matrix $a$ null ist, degeneriert die PDE zu einer PDE erster Ordnung

  $$\sum_{i=1}^d b_i \frac{\partial u}{\partial x_i} + c u = f$$ (eq:exmphyperbolclass)
  
  Auf einem Teil des Randes sind in dem Fall Randwerte erforderlich.

Diese Fälle verhalten sich mathematisch sehr unterschiedlich. Entsprechend unterschiedlich sind auch die physikalischen Anwendungen, welche oft zugrunde liegen. Mit dem Beispiel {eq}`eq:exmpparabolclass` kann ein zeitabhängiges Temperaturfeld und mit {eq}`eq:exmpparabolclass` die Wellenausbreitung modelliert werden.

