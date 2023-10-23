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

(chap:L2Raum)=
# Hilbertraum $L_2(\Omega)$

Da der Fokus in der Bearbeitung von partiellen Differentialgleichungen liegt, werden wir im Folgenden Funktionen auf beliebigen $n$ dimensionalen Gebieten $\Omega \subset \mathbb{R}^n$ und nicht nur zwingend auf Intervallen $[a,b] \subset \mathbb{R}$ betrachten.

Sei daher das Integrationsgebiet bzw. der Definitionsbereich der Funktion gegeben durch $\Omega$ eine beliebige (nichtleere) offene Menge in $\mathbb{R}^n$. Mit $C(\Omega)$ bezeichnen wir die Menge aller stetigen Funktionen auf $\Omega$.

```{prf:definition} Support
:label: my-def-support

Für $f\in C(\Omega)$ definieren wir den *Träger* oder *Support* von $f$ durch

$$\mathop{Tr} f := \overline{\{x\in\mathbb{R}^n | f(x) \not=0 \}}.$$

wobei mit $\overline{A}$ die Abschliessung einer Menge $A \subset \mathbb{R}^n$ bezeichnet.
```

Mit $C_0(\Omega)$ bezeichnen wir stetige Funktionen mit *beschränktem Support in* $\Omega$. Das Integral für $f\in C_0(\Omega)$ existiert: mit $a>0$ hinreichend gross gilt

$$\int_\Omega f(x)dx = \int_{-a}^a \left(\ldots\int_{-a}^a f(x_1,\ldots, x_n) dx_1 \ldots \right) dx_n.$$

```{figure} MehrdimensionalesIntegral.png
---
align: center
height: 250px
name: fig-mehrdimensionalesIntegral
---
Berechnung von $\int_\Omega f(x) dx$ in $\mathbb{R}^2$.
```

```{prf:definition} $C_0^\infty(\Omega)$ und $L_2(\Omega)$
:label: my-def-C0inftyL2

Mit $C_0^\infty(\Omega)$ bezeichnen wir die Menge aller Funktionen, welche in $\Omega$ beliebig oft stetig differenzierbar sind und einen beschränkten in $\Omega$ enthaltenen Support haben.

Mit der Quadratnorm

$$\|u\|_2 = \sqrt{\int_\Omega |u(x)|^2 dx}$$

definieren wir $L_2(\Omega)$ als den zu $(C_0^\infty(\Omega),\|\cdot\|_2)$ konjugierten (dualen) Raum:

$$L_2(\Omega) := (C_0^\infty(\Omega), \|\cdot\|_2)^*.$$
```

Man kann zeigen, dass sich die klassischen Funktionen aus $C_0^\infty(\Omega)$ in $L_2(\Omega)$ wiederfinden. Mathematisch bedeutet das, dass $C_0^\infty(\Omega)$ als Unterraum von $L_2(\Omega)$ aufgefasst werden kann.

$$C_0^\infty(\Omega) \subset L_2(\Omega)$$

Man sagt in dem Fall, $C_0^\infty(\Omega)$ ist in $L_2(\Omega)$ *eingebettet*. Es gilt

```{prf:theorem}
:label: my-thm-L2

Es gilt:
* $C_0^\infty(\Omega)$ liegt dicht in $L_2(\Omega)$
* $\overline{C_0^\infty(\Omega)}$ ist vollständig.
* $L_2(\Omega) = \overline{C_0^\infty(\Omega)}$
```
