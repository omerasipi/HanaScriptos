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

# Iterative Gleichungslöser

In Bezug auf Rechenzeit und **Speicherbedarf** oft viel effizientere Alternativen zur direkten Lösung von linearen Gleichungssystemen sind **iterative Verfahren** (insbesondere im 3D). Ausgehend von einem Startwert wird die näherungsweise Lösung iterativ verbessert. Die Verfahren benötigen nur Matrix-Vektorprodukte $A\cdot x$.

Werden von der Matrix $A$ nur die Nicht-Null Elemente gespeichert, so ist der Rechenaufwand für das Matrix-Vektorprodukt lediglich $O(N)$. Als mögliches Format für dünn besetzte Matrizen (sparse matrices) sei das [Compressed Sparse Row (CSR)](https://de.wikipedia.org/wiki/Compressed_Row_Storage) (vgl. auch [scipy sparse](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html)) erwähnt. Wir werden dieses in den Anwendungen benutzen.


In den folgenden Anwendungen unterschiedlicher iterativer Verfahren betrachten wir das Randwertproblem auf dem Einheitsquadrat

$$-\Delta u + 10\, u = 1\quad x\in\Omega = [0,1]^2$$

mit Dirichlet Randwerte $u=0$ auf $\partial\Omega$. Das schwache Problem lautet: gesucht ist $u\in H_0^1(\Omega)$:

$$\int_\Omega \big(\nabla u \nabla v\, + 10\,u\, v\big) dx = \int_\Omega 1\,v\,dx\quad \forall\ v\in H_0^1(\Omega).$$

```{tableofcontents}
```

```{glue:figure} FEM_RichardsonVerfahren_fig
---
figwidth: 400px
name: fig-FEM_RichardsonVerfahren_fig
---

Konvergenz Richardson Verfahrens
```

```{glue:figure} FEM_JacobiVerfahren_fig
---
figwidth: 400px
name: fig-FEM_JacobiVerfahren_fig
---

Konvergenz Jacobi Verfahrens
```

```{glue:figure} FEM_Gauss-SeidelVerfahren_fig
---
figwidth: 400px
name: fig-FEM_Gauss-SeidelVerfahren_fig
---

Konvergenz Gauss-Seidel Verfahrens
```

```{glue:figure} FEM_GradientenVerfahren_fig
---
figwidth: 400px
name: fig-FEM_GradientenVerfahren_fig
---

Konvergenz Gradienten Verfahrens
```

```{glue:figure} FEM_CGVerfahren_fig
---
figwidth: 400px
name: fig-FEM_CGVerfahren_fig
---

Konvergenz CG Verfahrens
```
