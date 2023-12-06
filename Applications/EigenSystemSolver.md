---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Eigenwertprobleme

```{code-cell} ipython3
:tags: [hide-cell, remove-output]

from ngsolve import *
from ngsolve.webgui import Draw
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
from numpy import random
```

## Rayleigh-Quotient

Wir betrachten das Eigenwertproblem

$$-\Delta u = \lambda u$$ (eq:eigenproblem)

auf dem Einheitsquadrat für unterschiedliche Randbedingungen
- Dirichlet Randbedingung

  $$u(x) = 0\qquad\text{für}\ x\in \partial\Omega$$

- und die vom Parameter $\alpha$ abhängige Robin Randbedingung

  $$\partial_n u(x) + \alpha\,u(x) = 0\qquad\text{für}\ x\in \partial\Omega.$$

Diese Probleme besitzen abzählbar viele positive Eigenwerte

$$0 < \lambda_1 < \lambda_2 \le \lambda_3 \le \ldots \lambda_n \to \infty.$$

Die Grundfrequenz der Dirichlet- und Robin-Probleme spielt eine zentrale Rolle in der Analysis und in physikalischen Anwendungen. Sie kann durch das Rayleigh-Prinzip charakterisiert werden

$$\lambda_1^D = \min_{u\in H_0^1(\Omega)} R(u) = \min_{u\in H_0^1(\Omega)} \frac{\int_\Omega |\nabla u|^2 dx}{\int_\Omega u^2 dx}$$ (eq:RayleighQuotientDirichlet)

für das Dirichlet Eigenwertproblem und

$$\lambda_1^R(\alpha) = \min_{u\in H^1(\Omega)} R(\alpha, u) = \min_{u\in H_0^1(\Omega)} \frac{\int_\Omega |\nabla u|^2 dx + \alpha \int_{\partial\Omega} u^2 ds}{\int_\Omega u^2 dx}$$ (eq:RayleighQuotientRobin)

für das Robin Eigenwertproblem. In beiden Fällen erhalten wir die erste Eigenfunktion.

## Numerische Berechnung von Eigenwerten

### Inverse Iteration

Die Diskretisierung des Randwertproblems {eq}`eq:eigenproblem` liefert das allgemeine Eigenwertproblem

$$A\cdot u + \lambda\, M\cdot u = 0,$$

im $\mathbb{R}^{n}$ wobei $A$ und $M$ gegeben sind durch die Matrizen

$$\begin{split}
A_{i,j} & = \int_\Omega \nabla\varphi_i(x)\cdot\nabla\varphi_j(x) dx\qquad\text{Steifigkeitsmatrix}\\
M_{i,j} & = \int_\Omega \varphi_i(x)\cdot\varphi_j(x) dx\qquad\text{Massenmatrix}\\
\end{split}$$

mit $i,j \in I_{\text{free dofs}}$. Die Inverse Iteration (**INVIT**)

$$u^{(k+1)} = A^{-1} M u^{(k)}$$

liefert den Rayleigh Quotient

$$\lambda_1^{(k)} = \frac{\langle A u^{(k)}, u^{(k)}\rangle}{\langle M u^{(k)}, u^{(k)}\rangle},$$

daher den kleinsten Eigenwert.

```{code-cell} ipython3
mesh = Mesh(unit_square.GenerateMesh(maxh=0.02))
```

Für das **Robin** Eigenwertproblem

```{code-cell} ipython3
:tags: [remove-output]

VR = H1(mesh, order=4)
uR,vR = VR.TnT()
gfuR = GridFunction(VR)

alpha = Parameter(1)

aR = BilinearForm(VR)
aR += grad(uR)*grad(vR)*dx
aR += alpha*uR*vR*ds

mR = BilinearForm(VR)
mR += uR*vR*dx

aR.Assemble()
mR.Assemble();
```

Wir starten die Inverse Iteration mit einem zufälligen Startwert:

```{code-cell} ipython3
:tags: [remove-output]

# Zufällige Werte für den Lösungsvektor
gfuR.vec.data = random.rand(VR.ndof)

tol = 1e-16
errR = []
lamR = []
r = gfuR.vec.CreateVector()
for k in range(20):
    with TaskManager():
        r.data = mR.mat*gfuR.vec
        gfuR.vec.data = aR.mat.Inverse(freedofs=VR.FreeDofs())*r
        uAu = InnerProduct(aR.mat*gfuR.vec,gfuR.vec)
        uMu = InnerProduct(mR.mat*gfuR.vec,gfuR.vec)
        lamR.append(uAu/uMu)
        r.data = aR.mat*gfuR.vec-lamR[-1]*mR.mat*gfuR.vec
        errR.append(Norm(r))
    print(k,lamR[-1],errR[-1])
    if errR[-1]/errR[0] < tol:
        break
errR = np.array(errR)
lamR = np.array(lamR)

# normieren der Lösung
gfuR.vec.data *= 1/np.sqrt(uMu)
```

```{code-cell} ipython3
Draw(gfuR)
```

Für das **Dirichlet** Eigenwertproblem

```{code-cell} ipython3
:tags: [remove-output]

VD = H1(mesh, order=4, dirichlet='bottom|right|top|left')
uD,vD = VD.TnT()
gfuD = GridFunction(VD)

aD = BilinearForm(VD)
aD += grad(uD)*grad(vD)*dx

mD = BilinearForm(VD)
mD += uD*vD*dx

aD.Assemble()
mD.Assemble();
```

Wir starten die Inverse Iteration mit einem zufälligen Startwert:

```{code-cell} ipython3
gfuD.vec.data = random.rand(VD.ndof)
```

```{code-cell} ipython3
:tags: [remove-output]

tol = 1e-16
errD = []
lamD = []
r = gfuD.vec.CreateVector()
for k in range(20):
    with TaskManager():
        r.data = mD.mat*gfuD.vec
        gfuD.vec.data = aD.mat.Inverse(freedofs=VD.FreeDofs())*r
        uAu = InnerProduct(aD.mat*gfuD.vec,gfuD.vec)
        uMu = InnerProduct(mD.mat*gfuD.vec,gfuD.vec)
        lamD.append(uAu/uMu)
        r.data = aD.mat*gfuD.vec-lamD[-1]*mD.mat*gfuD.vec
        errD.append(Norm(r))
    print(k,lamD[-1],errD[-1])
    if errD[-1]/errD[0] < tol:
        break
errD = np.array(errD)
lamD = np.array(lamD)

# normieren der Lösung
gfuD.vec.data *= 1/np.sqrt(uMu)
```

```{code-cell} ipython3
Draw(gfuD)
```

```{code-cell} ipython3
:tags: [hide-input]

plt.semilogy(errD, label='Residuum Dirichlet-Eigenvektor')
plt.semilogy(lamD-(2*np.pi**2), label='Fehler Dirichlet-Eigenwert')
plt.semilogy(errR, label='Residuum Robin-Eigenvektor')
plt.legend()
plt.grid()
plt.xlabel('Iteration')
plt.ylabel('Residuum / Fehler')
plt.xticks(range(0,errD.shape[0]))
plt.show()
```

### Preconditioned Inverse Iteration (PINVIT)

In der Inversen Iteration beinhaltet jede Iteration das Lösen des Linearen Gleichungssystems

$$A u^{(k+1)} = M u^{(k)},$$

was numerisch teuer ist. Die Preconditioned Inverse Iteration (**PINVIT**) vermeidet diesen Schritt in dem eine approximierte Inverse $C^{-1}$ benutzt wird:

$$\begin{split}
\lambda_1^{(k)} & = \frac{\langle A u^{(k)}, u^{(k)}\rangle}{\langle M u^{(k)}, u^{(k)}\rangle}\\
w^{(k)} & = C^{-1} (A\, u^{(k)} - \lambda_1^{(k)} M\, u^{(k)})\\
u^{(k+1)} & = u^{(k)} + \beta\, w^{(k)}
\end{split}$$

Die optimale Schrittweite $\beta$ kann durch Minimierung des Rayleigh-Quotienten in einem zweidimensionalen Raum gefunden werden:

$$u^{(k+1)} = \underset{v\in\{u^{(k)},w^{(k)}\}}{\text{arg min}} \frac{\langle A v, v\rangle}{\langle M v, v\rangle}$$

Dieses Minimierungsproblem kann durch ein kleines Eigenwertproblem

$$\tilde{A} y = \lambda \tilde{M} y$$

mit den Matrizen

$$\tilde{A} = \begin{pmatrix}
\langle A u^{(k)},u^{(k)}\rangle & \langle A u^{(k)},w^{(k)}\rangle\\
\langle A w^{(k)},u^{(k)}\rangle & \langle A w^{(k)},w^{(k)}\rangle\end{pmatrix},\qquad \tilde{M} = \begin{pmatrix}
\langle M u^{(k)},u^{(k)}\rangle & \langle M u^{(k)},w^{(k)}\rangle\\
\langle M w^{(k)},u^{(k)}\rangle & \langle M w^{(k)},w^{(k)}\rangle\end{pmatrix}$$

gelöst werden. Der neue Vektor ist gegeben durch

$$u^{(k+1)} = y_1\,u^{(k)} + y_2\,w^{(k)},$$

wobei $(y_1,y_2)^T$ der zum kleinsten Eigenwert zugehörige Eigenvektor sei. Wir wenden das Verfahren wiederum auf beide Eigenwertprobleme an.

```{code-cell} ipython3
def PINVIT(gfu, a, m, pre, maxIter=40, tol=1e-16):
    r = gfu.vec.CreateVector()
    w = gfu.vec.CreateVector()
    Mu = gfu.vec.CreateVector()
    Au = gfu.vec.CreateVector()
    Mw = gfu.vec.CreateVector()
    Aw = gfu.vec.CreateVector()
    
    proj = Projector(gfu.space.FreeDofs(), True)
    r.data = random.rand(gfu.space.ndof)
    gfu.vec.data = proj * r
    
    err = []
    with TaskManager():
        for k in range(maxIter):
            Au.data = a.mat * gfu.vec
            Mu.data = m.mat * gfu.vec
            auu = InnerProduct(Au, gfu.vec)
            muu = InnerProduct(Mu, gfu.vec)

            # Rayleigh quotient
            lam = auu/muu

            # residual
            r.data = Au - lam * Mu
            w.data = pre.mat * r.data
            w.data = 1/Norm(w) * w
            Aw.data = a.mat * w
            Mw.data = m.mat * w

            # setup and solve 2x2 small eigenvalue problem
            asmall = np.zeros((2,2))
            asmall[0,0] = auu
            asmall[0,1] = asmall[1,0] = InnerProduct(Au, w)
            asmall[1,1] = InnerProduct(Aw, w)

            msmall = np.zeros((2,2))
            msmall[0,0] = muu
            msmall[0,1] = msmall[1,0] = InnerProduct(Mu, w)
            msmall[1,1] = InnerProduct(Mw, w)

            eval,evec = scipy.linalg.eigh(a=asmall, b=msmall)

            gfu.vec.data = float(evec[0,0]) * gfu.vec + float(evec[1,0]) * w
            
            err.append(Norm(proj*r))
            if err[-1] < tol:
                break
    return lam, err
```

Für das **Robin** Eigenwertproblem:

Als Vorkonditionierer benutzen wir einen Multigrid Vorkonditionierer. Um diesen korrekt zur Verfügung zu haben, müssen wir die Bilinearform $a$ neu Assemblieren.

```{code-cell} ipython3
:tags: [remove-output]

preR = Preconditioner(aR, "multigrid")
aR.Assemble();
```

```{code-cell} ipython3
lamR,errR2 = PINVIT(gfuR,aR,mR,preR)
```

Für das **Dirichlet** Eigenwertproblem:

Ebenso für den Dirichlet-Fall müssen wir die Bilinearform neu Assemblieren, um den Vorkonditionierer korrekt zur Verfügung zu haben.

```{code-cell} ipython3
:tags: [remove-output]

preD = Preconditioner(aD, "multigrid")
aD.Assemble();
```

```{code-cell} ipython3
lamD,errD2 = PINVIT(gfuD,aD,mD,preD)
lamD
```

```{code-cell} ipython3
lamD-2*(np.pi**2)
```

```{code-cell} ipython3
:tags: [hide-input]

plt.semilogy(errD, label='Residuum Dirichlet - INVIT')
plt.semilogy(errD2, label='Residuum Dirichlet - PINVIT')
plt.semilogy(errR, '--', label='Residuum Robin - INVIT')
plt.semilogy(errR2, '--', label='Residuum Robin - PINVIT')
plt.legend()
plt.grid()
plt.xlabel('Iteration')
plt.ylabel('Residuum')
plt.show()
```

### Numerische Berechnung von mehreren Eigenwerte und -vektoren

Das Verfahren kann erweitert werden, um mehrere Eigenwerte zu berechnen. Dazu initialisieren wir eine GridFunction mit mehrfachen Komponenten für das Speichern der Eigenvektoren. Entsprechend werden alle Vektoren mit zufälligen Werte initialisiert. Der wesentliche Unterschied besteht nun darin, dass wir ein **kleines** Eigenwertproblem der Dimension $2 \times$Anzahl Eigenwerte lösen.

```{code-cell} ipython3
def PINVIT(gfu, a, m, pre, maxIter=20, tol=1e-16):
    num = len(gfu.vecs)
    r = gfu.vec.CreateVector()
    Av = gfu.vec.CreateVector()
    Mv = gfu.vec.CreateVector()
    proj = Projector(gfu.space.FreeDofs(), True)
    
    vecs = []
    for i in range(2*num):
        vecs.append (gfu.vec.CreateVector())

    for v in gfu.vecs:
        r.data = random.rand(gfu.space.ndof)
        v.data = proj * r
    
    asmall = np.zeros((2*num, 2*num))
    msmall = np.zeros((2*num, 2*num))
    lams = num * [1]

    with TaskManager():
        for k in range(maxIter):
            maxErr = 1e4
            for j in range(num):
                vecs[j].data = gfu.vecs[j]
                r.data = a.mat * vecs[j] - lams[j] * m.mat * vecs[j]
                vecs[num+j].data = pre.mat * r
                err = Norm(proj*r)
                if err < maxErr:
                    maxErr = err

            for j in range(2*num):
                Av.data = a.mat * vecs[j]
                Mv.data = m.mat * vecs[j]
                for i in range(2*num):
                    asmall[j,i] = InnerProduct(Av, vecs[i])
                    msmall[j,i] = InnerProduct(Mv, vecs[i])

            ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)
            lams[:] = ev[0:num]

            for j in range(num):
                gfu.vecs[j][:] = 0.0
                for i in range(2*num):
                    gfu.vecs[j].data += float(evec[i,j]) * vecs[i]
            if maxErr < tol:
                break
    
    return np.array(lams)
```

Damit folgt für die 5 ersten Eigenwerte:

```{code-cell} ipython3
num = 5
gfuR = GridFunction(VR, multidim=num)
gfuD = GridFunction(VD, multidim=num)
```

Für das **Robin**-Eigenwertproblem:

```{code-cell} ipython3
lamR = PINVIT(gfuR,aR,mR,preR)
lamR
```

Wir stellen noch sicher, dass die erste Eigenfunktion positiv ist:

```{code-cell} ipython3
gfuR.vecs[0].data *= np.sign(gfuR(mesh(0.5,0.5)))
Draw(gfuR)
```

Für das **Dirichlet**-Eigenwertproblem:

```{code-cell} ipython3
lamD = PINVIT(gfuD,aD,mD,preD)
lamD
```

```{code-cell} ipython3
gfuD.vecs[0].data *= np.sign(gfuD(mesh(0.5,0.5)))

Draw(gfuD)
```

## Zurück zum Rayleigh-Quotient

Der Rayleigh-Quotient für das Robin-Eigenwertproblem {eq}`eq:RayleighQuotientRobin` ist vom Parameter $\alpha$ abhängig. Für $\alpha \to \infty$ konvergiert $\lambda_1^R(\alpha) \to \lambda_1^D$, dem Rayleigh-Quotienten für das Dirichlet-Eigenwertproblem {eq}`eq:RayleighQuotientDirichlet`. Wir berechnen diese Abhängigkeit für das Einheitsquadrat:

```{code-cell} ipython3
:tags: [hide-input]

alphas = np.linspace(1e-4,500,200)
lamRalphas = []
gfuR = GridFunction(VR)

for ai in alphas:
  alpha.Set(ai)
  with TaskManager():
    aR.Assemble()
  lamRalphas.append(PINVIT(gfuR,aR,mR,preR)[0])

plt.plot(alphas,lamRalphas)
plt.axhline(lamD[0],c='gray')
plt.ylim(0,22)
plt.xlim(-4,500)
plt.grid()
plt.text(503,lamD[0],r'$\lambda_1^D(\Omega_R) = 2 \pi^2$')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\lambda_1^R(\alpha)$')
plt.show()
```

```{prf:remark}
Die Rayleigh-Faber-Krahn-Ungleichung besagt, dass unter allen Gebieten $\Omega$ mit gegebenem Volumen $\lambda_1(\Omega)$ minimal für die Kugel ist. Es kann sogar gezeigt werden, dass die Kugel das einzige minimale Gebiet ist.
```

```{prf:remark}
Mit einem $\alpha$ gross genug können Dirichlet-Randwerte mit Robin-Randwerte gut approximiert werden. Dies gilt nicht nur für das Eigenwertproblem sondern im Allgemeinen. 
```
