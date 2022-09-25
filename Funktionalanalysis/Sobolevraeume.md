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

# Sobolevräume

Im Kapitel stellen wir die für die Behandlung von partiellen Differentialgleichungen geeigneten *Funktionenräume* vor.
In den Beispielen für Hilberträume wurde erwähnt, dass der Funktionenraum der stetigen Funktionen versehen mit der Quadratnorm $\|x\|_2 = \sqrt{\int_a^b |x(t)|^2 dt}$ nicht vollständig und damit kein Banachraum ist. Das Problem muss für die Behandlung der partiellen Differentialgleichungen gelöst werden und führt uns zu den sogenannten Sobolevräumen.
