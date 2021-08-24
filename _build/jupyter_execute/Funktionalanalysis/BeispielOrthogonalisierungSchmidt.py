#!/usr/bin/env python
# coding: utf-8

# # Beispiel Orthogonalisierungsverfahren
# 
# nach Erhard Schmidt

# In[1]:


import numpy as np
from sympy import symbols, integrate, lambdify
import matplotlib.pyplot as plt


# Der Raum der quadratisch integrierbaren (nach Lebesgues) Funktionen $L_2[-1,1]$ ist mit dem Skalarprodukt
# 
# $$(x,y) = \int_{-1}^{1} x(t) y(t) dt$$
# und der induzierten Norm
# 
# $$\|x\|_2 = \sqrt{(x,x)}$$
# ein Hilbertraum. Wir definieren daher das Skalarprodukt (dot-product) und die norm wie folgt:

# In[2]:


t = symbols('t')
def dot(x,y):
    return integrate(x*y,(t,-1,1))
def norm(x):
    return dot(x,x)**(1/2)


# Wir betrachten die Folge $\{t^n\}_{n\in\mathbb{N}}\subset L^2[-1,1]$ von Monomen:

# In[3]:


yi = []
for i in range(5):
    yi.append(t**i)
print(yi)


# Die Monome sind nicht orthogonal, jedoch linear unabhängig. Falls dem so wäre müsste die Einheitsmatrix entstehen:

# In[4]:


m = []
for i in range(5):
    mi = []
    for j in range(5):
        mi.append(dot(yi[i],yi[j]))
    m.append(mi)
m


# Wir berechnen nun ein orthonormales System, basierend auf den Monomen nach dem Verfahren von Schmidt:

# In[5]:


xi = [yi[0]/norm(yi[0])]
for i in range(1,5):
    s = 0
    for j in range(i):
        s += dot(yi[i],xi[j])*xi[j]
    zi=yi[i]-s
    xi.append(zi/norm(zi))


# Das Orthonormalsystem ist damit gegeben durch:

# In[6]:


xi


# Test des Orthonormalsystems ergibt die Einheitsmatrix:

# In[7]:


m = []
for i in range(5):
    mi = []
    for j in range(5):
        mi.append(dot(xi[i],xi[j]))
    m.append(mi)
np.round(np.array(m,dtype=float),8)


# In[8]:


tp = np.linspace(-1,1,400)
n = 0
plt.plot(tp,np.ones_like(tp),label='$P_'+str(n)+'(t)$')
for xii in xi[1:]:
    n += 1
    f = lambdify(t,xii,'numpy')
    plt.plot(tp, f(tp),label='$P_'+str(n)+'(t)$')
plt.grid()
plt.legend()
plt.show()


# Das Resultat sind bezüglich der $L_2$-Norm normierte Legendre'sche Polynome. Üblicher weise werden die Polynome mit $P_n(1) = 1$ normiert:

# In[9]:


tp = np.linspace(-1,1,400)
n = 0
plt.plot(tp,np.ones_like(tp),label='$P_'+str(n)+'(t)$')
for xii in xi[1:]:
    n += 1
    f = lambdify(t,xii,'numpy')
    plt.plot(tp, f(tp)/f(1),label='$P_'+str(n)+'(t)$')
plt.grid()
plt.legend()
plt.show()


# Die ersten fünf Legendre'sche Polynome sind gegeben durch:

# In[10]:


lp = [xii/xii.subs(t,1) for xii in xi]
print(lp)


# In[ ]:




