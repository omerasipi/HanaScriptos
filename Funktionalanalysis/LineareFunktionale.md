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

(chap:lineareFunktionale)=
# Lineare Funktionale

Wir betrachten nun spezielle lineare Operatoren, welche in den zugrunde liegenden Zahlenkörper abbilden.

```{prf:definition} lineares Funktional
:label: my-def-linfunktional

Sei $V$ ein normierter Raum. Dann nennt man den linearen Operator $F: V \to \mathbb{K}$ ($\mathbb{R}$ oder $\mathbb{C}$) *lineares Funktional*.
```

Wir werden den Begriff **Linearform** für lineare Funktionale und die **Bilinearform** im Folgenden häufig antreffen. 

```{prf:definition} Linearform, Bilinearform
:label: my-def-linbilinform

* Linearform $f(\cdot)$ ist eine Abbildung

  $$\begin{split} f: V & \to \mathbb{K}\\
  v & \mapsto f(v)\quad \text{linear}\end{split}$$  (eq:Linearform)

* Bilinearform $A(\cdot, \cdot)$ ist eine Abbildung

  $$\begin{split} A: V \times V & \to \mathbb{K}\\
  (u,v) & \mapsto A(u,v)\quad \text{linear in $u$ und $v$.}\end{split}$$ (eq:Bilinearform)
```

```{prf:remark}
:label: my-rm-linearform

* Die Linearform ist ein äquivalenter Begriff für ein lineares Funktional.

* Der Stetigkeitsbegriff linearer Operatoren auf das lineare Funktional bzw. Linearform und auf die Bilinearform angewandt bedeutet:

  * die Linearform $f$ heisst *stetig*, falls

    $$\exists\ C>0:\quad |f(v)| \le C\ \|v\|_V\quad \forall\ v\in V$$

  * die Bilinearform $A(u,v)$ heisst *stetig*, falls

    $$\exists\ C>0:\quad |A(u,v)| \le C\ \|u\|_V\,\|v\|_V\quad \forall\ u,v\in V.$$
```

Normiert man $\mathbb{K}$ mit

$$\|\alpha\| := |\alpha|,\quad\alpha\in\mathbb{K}$$

so ist $\mathbb{K}$ ein Banachraum und damit die Menge aller beschränkten linearen Funktionale auf $V$ ein Banachraum. Dieser Raum ist insbesondere im Zusammenhang mit partiellen Differentialgleichungen sehr wichtig.

```{admonition} Dualraum von $V$
Der Banachraumm $L(V,\mathbb{K})$ aller beschränkten linearen Funktionale auf $V$ heisst der zu $V$ *konjugierte* oder *duale Raum* und wird mit $V^*$ oder $V'$ bezeichnet.
```

**Beispiel**: Sei $V$ ein Hilbertraum und $y_0$ ein beliebiges (festes) Element aus $V$. Für $x\in V$ wird durch

$$\begin{split}
F : V & \to \mathbb{R}\\
 x & \mapsto y = F x := (x,y_0)\end{split}$$

ein lineares Funktional $F$ definiert. (Die Linearität folgt direkt aus den Eigenschaften des Skalarprodukts.) Mit Hilfe der Schwarzschen Ungleichung folgt

$$\|F x\| = |F x| = |(x,y_0)| \le \|x\|\,\|y_0\|\quad\forall\ x\in V.$$

Damit folgt

$$\frac{\|F x\|}{\|x\|} \le \|y_0\|\quad \forall x\in V,\ \text{mit}\,x\not= 0$$

sprich $F$ ist ein beschränktes lineares Funktional, $F\in V^*$ mit $\|F\| \le \|y_0\|$. Da für $x=y_0$

$$\|F y_0\| = |(y_0,y_0)| = \|y_0\|^2 = \|y_0\|\, \|y_0\|$$

gilt, folgt

$$\|F\| = \sup_{\substack{x\in V\\x\not= 0}} \frac{\|F x\|}{\|x\|} = \|y_0\|.$$

Wir kommen nun zum Rieszschen Darstellungssatz: Beschränkte lineare Funktionale eines Hilbertraumes $V$ lassen sich besonders einfach darstellen. Die Darstellung aus obigem Beispiel erfasst **alle** beschränkten linearen Funktionale. Es gilt

```{prf:theorem} Darstellungssatz von Riesz
:label: my-thm-Riesz

Sei $V$ ein Hilbertraum und $F\in V^*$ beliebig. Dann gibt es ein *eindeutig* bestimmtes $y\in V$, so dass $F$ die Darstellung

$$F x = (x,y)\quad\forall x\in V$$

besitzt.
```

```{prf:remark}
:label: my-rm-ergebnishilbertraumtheorie

Dieser Satz ist das zentrale Ergebnis der Hilbertraum-Theorie. Neben seiner Bedeutung als Darstellungssatz kann er auch als Existenz- und Eindeutigkeitsprinzip aufgefasst werden. Diese Bedeutung des Rieszschen Satzes ist Grundlage für die moderne Theorie der elliptischen partiellen Differentialgleichungen (vgl. auch {numref}`chap:konvergenzanalyse`).
```
