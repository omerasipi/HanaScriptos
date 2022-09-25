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

# Klassifikation 

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
