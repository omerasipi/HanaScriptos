# Numerik für partielle Differentialgleichungen

Gängige Methoden für die numerische Berechnung von Lösungen partieller Differentialgleichungen sind

* Methode der finite Differenzen (FDM)

  Die Umsetzung ist in der Regel einfach, wobei vorausgesetzt wird, dass das Gebiet rechteckförmig ist.

* Methode der finiten Elemente (FEM)

  Die FEM ist sehr flexibel was komplexe Geometrien und Randbedingungen betrifft ebenso können Methoden höhrer Ordnung hergeleitet werden. Die Analysis der FEM basiert direkt auf der Funktionalanalysis. Die Implementierung effizienter Methoden erfordert in der Regel versierte Programmiertechniken.

* Methode der finiten Volumen (FVM)

  Die FVM sind flexibel im Umgang mit komplexen Geometrien und komplizierten Randbedingungen. Sie erhalten physikalische Gesetze auch in diskretisierter Form und sind daher sehr beliebt für das Lösen von Strömungsprobleme. Der Nachteil besteht darin, dass sich nur schwer Verfahren höherer Ordnung herleiten lassen.

* Methode der Randelemente (BEM)

  Die "Boundary Element Method" reduziert das Problem um eine Dimension. Sind daher jedoch typischerweise beschränkt auf lineare elliptische und parabolische Gleichungen. Die Implementierung ist wie bei FEM in der Regel nicht einfach und erfordert mehr mathematisches Wissen eine gute und äquivalente Integralform zu finden. Neuste Verfahren sind jedoch sehr effizient.

* Spektrale Methoden

  Die spektralen Methoden setzen üblicherweise auch eine einfache Geometrie voraus. Die Behandlung nichtlinearer Probleme sowie komplexer Geometrien ist nicht einfach. 

In diesem Kurs beschränken wir uns auf die **Methode der finiten Elemente**.
