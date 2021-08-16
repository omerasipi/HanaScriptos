# Funktionalanalysis

In diesem Kapitel wird das notwendige Rüstzeug für die Behandlung partieller Differentialgleichungen und numerischen Methoden bereit gestellt.

## Grundlegende Räume

In den Modulen linearen Algebra, Analysis und Numerik wurden Vektorräume und Normen mit dem $\mathbb{R}^n$ eingeführt und benutzt. Wir definieren hier der Vollständigkeit halber die Begriffe nochmals und erweitern die Anwendung auf allgemeinere Räume, insbesondere müssen diese nicht endlich dimensional sein.

Beginnen wir mit dem Begriff der Vektorräume. Um diesen definieren zu können, benötigen wir einen Zahlenkörper. In aller Regel benutzen wir die reellen Zahlen $\mathbb{R}$. Wir legen mit dem folgenden Axiom fest, was die reellen Zahlen sind:

```{admonition} Definition: reelle Zahlen

Es existiere eine Menge $\mathbb{R}$ mit folgenden Eigenschaften:
- Es existieren Operationen

  $$\begin{split}  + : \mathbb{R} \times \mathbb{R} \ & \to \ \mathbb{R}\quad\text{Addition}\\
  \cdot : \mathbb{R} \times \mathbb{R} \ & \to \ \mathbb{R}\quad\text{Multiplikation}\end{split}$$
  mit den Eigenschaften:
  - Assoziativgesetze

    $$\begin{split}
    (a+b) + c & = a + (b+c)\quad \forall\ a,b,c \in \mathbb{R}\\
    (a\cdot b) \cdot c & = a \cdot (b\cdot c)\quad \forall\ a,b,c \in \mathbb{R}
    \end{split}$$
 
  - Kommutativgesetze
    
    $$\begin{split}
    a+b & = b + a \quad \forall\ a,b \in \mathbb{R}\\
    a\cdot b & = b \cdot a\quad \forall\ a,b \in \mathbb{R}
    \end{split}$$

  - Neutrale Elemente
    
    $$\begin{split}
    \exists 0: a + 0 = a\quad \forall a \in \mathbb{R}\\
    \exists 1: a \cdot 1 = a\quad \forall a \in \mathbb{R}
    \end{split}$$
  
  - Inverse Elemente
    
    $$\begin{split}
    \forall a\in \mathbb{R}\ \exists a'\in\mathbb{R} : a + a' = 0\\
    \forall a\in \mathbb{R}\setminus \{ 0\}\ \exists \tilde{a}\in\mathbb{R} : a \cdot \tilde{a} = 1
    \end{split}$$
    Schreibweise: $-a := a', a - a := a + (-a), \frac{1}{a} := \tilde{a}, \frac{a}{a} := a\cdot \frac{1}{a}$
  
  - Distributivgesetz
    $$a\cdot (b+c) = a\cdot b + a\cdot c$$

- $\mathbb{R}$ ist geordnet, d.h. es existiert eine Relation $\le$ so, dass gilt
  - $\mathbb{R}$ ist totalgeordnet, d.h.
    1. $\mathbb{R}$ ist teilgeordnet, d.h.
       - Für alle $x\in\mathbb{R}$ gilt: $x \le x$.
       - Ist $x \le y$ und $y\le z$ für $x,y,z\in\mathbb{R}$, so ist $x \le z$.
       - Ist $x \le y$ und $y \le x$ für $x,y \in \mathbb{R}$, so ist $x = y$.
       
       (Mit $x<y$ bezeichnen wir den Fall $x \le y$ und $x \not= y$.)
    2. Für je zwei Elemente $x,y \in\mathbb{R}$ gilt
       
       $$x \le y\quad \text{oder}\quad y\le x.$$
  - Die Ordnung ist verträglich mit Addition und Multiplikation, d.h. für $a,b,c\in\mathbb{R}$ gilt
  
  $$\begin{array}{c}
  a \le b \Rightarrow a+c \le b+c\\
  a \le b, 0 \le c \Rightarrow a\cdot c \le b\cdot c
  \end{array}$$

- $\mathbb{R}$ ist vollständig. Das heisst, dass jede nicht leere nach oben beschränkte Menge reeller Zahlen eine kleinste obere Schranke besitzt.
```

### Metrische Räume

```{admonition} Definition

Hier definieren wir eine Metrik.

```

### Normierte Räume. Banachräume

### Skalarprodukträume. Hilberträume

## Lineare Operationen

### Beschränkte lineare Operatoren
