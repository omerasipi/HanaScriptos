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

# Lineare Operatoren

Viele Aufgaben in der Mathematik und den Anwendungen führen auf Gleichungen der Form

$$T(x) = T x = y,$$

wobei $T: V \to W$ eine "lineare Abbildung", $V, W$ normierte Räume sind und $y\in W$ ein gegeben ist.

```{prf:definition} lineare Abbildung
:label: my-def-linabb

Die _Abbildung_ $T$ des normierten Raumes $V$ in den normierten Raum $W$ heisst _linear_, wenn für alle $x,y \in V$ und alle $\alpha \in \mathbb{K}$

$$T(x+y) = T x + T y$$
```
