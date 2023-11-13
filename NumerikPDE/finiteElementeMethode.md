# Effiziente Berechnung der FEM-Systemmatrix

Die Idee der Methode der finiten Elemente ist nun bekannt. Im einführenden Beispiel wurde direkt mit Hilfe der sogenannten globalen Basisfunktionen die Systemmatrix bzw. Rechteseite berechnet. Die Frage stellt sich jedoch, wie man diese möglichst universal anwendbar, daher für beliebige Gebietszerlegungen beschreiben kann, ohne komplizierte Konstruktion der (globalen) Basisfunktionen.

Wir stellen hier das Konzept des Assemblings mit Hilfe von lokalen Basisfunktionen auf Einheitsgebieten vor. Im Eindimensionalen benutzt man als Einheitsgebiet typischerweise das Einheitsintervall und im Zweidimensionalen das Einheitsdreieck oder -viereck.

Im Kapitel wird das Vorgehen im 1d und 2d detailliert mit Hilfe von numpy betrachtet. Parallel dazu wird die Implementierung ebenso mit NGSolve gezeigt.

```{tableofcontents}
```
