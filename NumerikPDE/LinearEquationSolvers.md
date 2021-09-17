# Lineare Gleichungslöser

Die Diskretisierung von linearen partiellen Differentialgleichungen führt zu **linearen Gleichungssystemen**

$$A\,x = b.$$

Der Rechenaufwand beim numerischen Lösen von PDE's steckt daher im
* Assembling (Berechnen der Matrizen und Vektoren)
* Lösen des linearen Gleichungssystems.

Im Fall, dass die Gleichung nichtlinear ist, muss z.B. durch anwenden des Newton-Gauss Verfahrens mehrfach ein lineares Gleichungssystem gelöst werden. Die Aufgabe bleibt daher auch in dem Fall die selbe.

Typische Eigenschaften der Gleichungssysteme sind
* Sehr **grosse Dimension $N$**. Tägliche Praxis für 3D Probleme ist ca. 10^6, bei Spezialanwendungen 10^9-10^12.
* Die Matrix hat **Bandstruktur**. Direkte Gleichungslöser wie z.B. LR-Verfahren für tridiagonal Matrizen haben mit Verwendung der Bandstruktur einen Rechenaufwand von 1D $O(N)$, 2D $O(N^2)$ und in 3D $O(N^{7/3})$.
* Die Matrix ist oft **symmetrisch und positiv definit**, das ist bei vielen FEM Diskretisierungen von elliptischen, parabolischen und hyperbolischen PDE's der Fall.

## Direkte Gleichungslöser

Die FEM Matrix der Dimension hat nach geeignetem Umordnen (Nummerierung der Knoten, Dreiecke, etc.) oft eine Block-Struktur. Im Gegensatz zu direkten Löser basierend auf Faktorisierungsmethoden, wie z.B. die oben erwähnte LR-Zerlegung, können Block-eliminations Methoden noch effizienter sein.

```{seealso}
[Beispiele zur Matrixstruktur](MatrixStrukturMesh.ipynb)
```

Wir teilen die Unbekannten in zwei Gruppen und schreiben die Gleichung $A\,u = f$ neu als Blocksystem

$$\begin{pmatrix}
A_{11} & A_{12}\\
A_{21} & A_{22}\end{pmatrix}\ \begin{pmatrix}u_1\\u_2\end{pmatrix} = \begin{pmatrix}f_1\\f_2\end{pmatrix}.$$

Als erstes berechnen wir $u_1$

$$u_1 = A_{11}^{-1}\,(f_1-A_{12}u_2)$$

und die Schur-Komplement Gleichung um $u_2$ zu bestimmen

$$\underbrace{(A_{22} - A_{21}A_{11}^{-1}A_{12})}_{=S}\,u_2 = f_2 - A_{21} A_{11}^{-1} f_1.$$

```{admonition} Definition: Schur-Komplement
Die Matrix

$$S = A_{22} - A_{21}A_{11}^{-1}A_{12}$$

heisst **Schur-Komplement** $S_{11}$ in $A$.
````

Diese Block-Faktorisierung wird in *sub-structuring* Algorithmen benutzt. Dabei zerlegt man das Gebiet in $m\times m$ Untergebiete und teilt die Unbekannten in innere und koppelnde Unbekannte auf. Nun kann die Faktorisierung der Blöcke mit den inneren Unbekannten und die Berechnung des Schur-Komplements  parallel berechnet werden. Die Eliminierung der internen Unbekannten kann mit sogenannten hierarchischen sub-structuring Algorithmen (z.B. nested dissection) durchgeführt werden.

```{figure} hierachischeElimination.png
---
figwidth: 400px
name: fig-hierachischeElimination
---

Hierarchische Elimination
```

Hierarchische Elimination funktioniert bei strukturierten Gitter sehr gut. Man kann sie jedoch auch unstrukturierten Gitter anwenden. In dem Fall basiert die Ordnung auf dem minimalen Grad der Kopplung. Nacheinander werden die Unbekannten mit den wenigsten Verbindungen im Matrixgraphen eliminiert.

In 2D ist eine direkte Methode mit optimaler Anordnung sehr effizient. In 3D ist die Situation für den direkten Löser schlechter.

```{note}
In der Praxis bedeutet das in aller Regel, dass 2D Probleme einfach mit direkten Löser berechnet werden können. Hingegen bei 3D Problemen benötigen wir zwingend einen anderen Ansatz.
```

    Help on method Inverse in module ngsolve.la:

    Inverse(...) method of ngsolve.la.SparseMatrixd instance
        Inverse(self: ngsolve.la.BaseMatrix, freedofs: pyngcore.BitArray = None, inverse: str = '') -> BaseMatrix
        
        Calculate inverse of sparse matrix
        Parameters:
        
        freedofs : BitArray
        If set, invert only the rows/columns the matrix defined by the bit array, otherwise invert the whole matrix
        
        inverse : string
        Solver to use, allowed values are:
            sparsecholesky - internal solver of NGSolve for symmetric matrices
            umfpack        - solver by Suitesparse/UMFPACK (if NGSolve was configured with USE_UMFPACK=ON)
            pardiso        - PARDISO, either provided by libpardiso (USE_PARDISO=ON) or Intel MKL (USE_MKL=ON).
                             If neither Pardiso nor Intel MKL was linked at compile-time, NGSolve will
                             look for libmkl_rt in LD_LIBRARY_PATH (Unix) or PATH (Windows) at run-time.


## Iterative Gleichungslöser

In Bezug auf Rechenzeit und **Speicherbedarf** oft viel effizientere Alternativen zur direkten Lösung von linearen Gleichungssystemen sind **iterative Verfahren** (insbesondere im 3D). Ausgehend von einem Startwert wird die näherungsweise Lösung iterativ verbessert. Die Verfahren benötigen nur Matrix-Vektorprodukte $A\cdot x$.

Werden von der Matrix $A$ nur die Nicht-Null Elemente gespeichert, so ist der Rechenaufwand für das Matrix-Vektorprodukt lediglich $O(N)$.

Im folgenden betrachten wir das Randwertproblem auf dem Einheitsquadrat

$$-\Delta u + 10\, u = 1\quad x\in\Omega = [0,1]^2$$

mit Dirichlet Randwerte $u=0$ auf $\partial\Omega$. Das schwache Problem lautet: gesucht ist $u\in H_0^1(\Omega)$:

$$\int_\Omega \big(\nabla u \nabla v\, + 10\,u\, v\big) dx = \int_\Omega 1\,v\,dx\quad \forall\ v\in H_0^1(\Omega).$$

```{seealso}
* [Richardson Verfahren](RichardsonVerfahren.ipynb)
* [Jacobi-Verfahren und Gauss-Seidel Verfahren](JacobiGauss_SeideVerfahren.ipynb)
* [Gradienten Verfahren](JacobiGauss_SeideVerfahren.ipynb)
* [Konjugiertes Gradienten Verfahren](GradientenVerfahren.ipynb)
* [CG - Verfahren](CGVerfahren.ipynb)
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
