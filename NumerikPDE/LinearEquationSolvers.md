# Lineare Gleichungslöser

Die Diskretisierung von linearen partiellen Differentialgleichungen führt zu **linearen Gleichungssystemen**

$$A\,x = b.$$

Der Rechenaufwand beim numerischen Lösen von PDE's steckt daher im
* Assembling (Berechnen der Matrizen und Vektoren)
* Lösen des linearen Gleichungssystems.

Im Fall, dass die Gleichung nichtlinear ist, muss z.B. durch anwenden des Newton-Gauss Verfahrens mehrfach ein lineares Gleichungssystem gelöst werden. Die Aufgabe bleibt daher auch in dem Fall die selbe.

Typische Eigenschaften der Gleichungssysteme sind
* Sehr **grosse Dimension $N$**. Tägliche Praxis für 3D Probleme ist ca. $10^6$, bei Spezialanwendungen $10^{9}-10^{12}$.
* Die Matrix hat **Bandstruktur**. Direkte Gleichungslöser wie z.B. LR-Verfahren für tridiagonal Matrizen haben mit Verwendung der Bandstruktur einen Rechenaufwand von 1D $O(N)$, 2D $O(N^2)$ und in 3D $O(N^{7/3})$.
* Die Matrix ist oft **symmetrisch und positiv definit**, das ist bei vielen FEM Diskretisierungen von elliptischen, parabolischen und hyperbolischen PDE's der Fall.
