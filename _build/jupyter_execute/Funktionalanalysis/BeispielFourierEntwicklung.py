#!/usr/bin/env python
# coding: utf-8

# # Beispiel zur Fourierentwicklung

# In[1]:


import numpy as np
from sympy import symbols, integrate, lambdify, exp, sin, cos, pi
import scipy.integrate as scint
import matplotlib.pyplot as plt


# Für das Beispiel benutzen wir zwei verschiedene orthonormale Basen für den $L_2[-1,1]$ 

# In[2]:


t = symbols('t')

def dot(x,y):
    return integrate(x*y,(t,-1,1))
def norm(x):
    return dot(x,x)**(1/2)
def dotN(x,y):
    xy = lambdify(t, x*y,'numpy')
    return scint.quad(xy, -1, 1)[0]
def normN(x):
    return dotN(x,x)**(1/2)


# ## Legendre'sche Polynome

# In[3]:


# Orthonormalbasis nach Schmidt
N = 9
yi = [t**i for i in range(N)]
xi = [yi[0]/norm(yi[0])]
for i in range(1,N):
    s = 0
    for j in range(i):
        s += dot(yi[i],xi[j])*xi[j]
    zi=yi[i]-s
    xi.append(zi/norm(zi))


# In[4]:


xi


# In[5]:


f = exp(-5*t**2)
fn = lambdify(t, f, 'numpy')


# In[6]:


tp = np.linspace(-1,1,400)
plt.plot(tp, fn(tp))
plt.grid()
plt.show()


# In[7]:


alpha = [dotN(f, xii) for xii in xi]
alpha


# In[8]:


s = 0
for alphai, xii in zip(alpha,xi):
    s += alphai*xii
sn = lambdify(t, s, 'numpy')
s


# In[9]:


tp = np.linspace(-1,1,400)
plt.plot(tp, fn(tp),label='f(t)') 
plt.plot(tp, sn(tp),label='Legendre Basis')
plt.legend()
plt.grid()
plt.show()


# Parsevallsche Gleichung (Vollständigkeitsrelation):
# 
# $$\sum_{k=1}^N |(x,x_k)|^2 \le \sum_{k=1}^\infty |(x,x_k)|^2 = \|x\|^2$$

# In[10]:


np.sum(np.array(alpha)**2)-normN(f)**2


# ## Trigonometrische Funktionen

# In[11]:


# Orthonormalbasis nach Schmidt
yi = []
N = 5
yi.append(1/2)
for i in range(1,N):
    yi.append(cos(pi*i*t))
    yi.append(sin(pi*i*t))
print('Anzahl Basisfunktionen: '+str(len(yi)))
xi = [yi[0]/normN(yi[0])]
for i in range(1,2*N-1):
    s = 0
    for j in range(i):
        s += dotN(yi[i],xi[j])*xi[j]
    zi=yi[i]-s
    xi.append(zi/normN(zi))


# In[12]:


xi


# In[13]:


alpha2 = [dotN(f, xii) for xii in xi]
alpha2


# In[14]:


s2 = 0
for alphai, xii in zip(alpha2,xi):
    s2 += alphai*xii
s2n = lambdify(t, s2, 'numpy')
s2


# In[15]:


tp = np.linspace(-1,1,400)
plt.plot(tp, fn(tp),label='f(t)') 
plt.plot(tp, sn(tp),label='Legendre Basis')
plt.plot(tp, s2n(tp),label='trigonometrische Basis')
plt.legend()
plt.grid()
plt.show()


# Parsevallsche Gleichung (Vollständigkeitsrelation):
# 
# $$\sum_{k=1}^N |(x,x_k)|^2 \le \sum_{k=1}^\infty |(x,x_k)|^2 = \|x\|^2$$

# In[16]:


np.sum(np.array(alpha2)**2)-normN(f)**2


# ## Vergleich

# In[17]:


plt.plot(tp, sn(tp)-fn(tp),label='Legendre Basis')
plt.plot(tp, s2n(tp)-fn(tp),label='trigonometrische Basis')
plt.legend()
plt.grid()
plt.show()

