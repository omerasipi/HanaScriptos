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

# Einführung partielle Differentialgleichungen

In der Technik und  Naturwissenschaft haben wir es oft mit partiellen Differentialgleichungen (PDE) zu tun: So kann das mechanische Schwingungsverhalten, die Ausbreitung von Licht, Strömungen um Objekte wie Autos, Flugzeuge, etc., elektromagnetische Felder von Transformatoren und anderen Bauteilen berechnet und für die technische Anwendung optimiert werden. Insbesondere die numerischen Methoden für PDE sind in der Ingenieurwissenschaft nicht mehr weg zu denken. In der Entwicklung werden oft ganze Systeme mathematisch modelliert, numerisch simuliert und optimiert, bevor der erste Prototyp erstellt wird. In den meisten Fällen baut die Modellierung auf einer Formulierung mit Hilfe von PDE auf. 

Die Theorie der PDE ist ein sehr umfangreiches mathematisches Gebiet, bei dem sehr unterschiedliche Verfahren und Methoden Anwendung finden. Wir werden uns in diesem Kurs insbesondere auf eine in der Anwendung sehr weit verbreitete numerische Methoden für PDE konzentrieren, die Methode der finiten Elemente.
