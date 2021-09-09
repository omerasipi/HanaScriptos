#!/usr/bin/env python
# coding: utf-8

# # Beispiele zu quadratischen Formen

# In[1]:


import matplotlib.pyplot as plt
import numpy as np


# In[2]:


from scipy.linalg import solve_triangular, eig, eigvals, norm, qr, svd


# In[3]:


from sympy import symbols, lambdify


# In[4]:


xp = np.linspace(-2,2,30)
Xp,Yp = np.meshgrid(xp,xp)


# In[5]:


def getQuadForm(a,b,c, lam = True):
    x,y = symbols('x,y')
    u = np.array([x,y])
    expr=u@a@u+b@u+c
    if lam:
        return lambdify((x,y),expr)
    else:
        return expr


# (chap:ellipticexmp)=
# ## Elliptisch

# Die Matrix $a$ muss positiv definit sein:

# In[6]:


a = np.array([[2,1],[0,1]])
b = np.array([0,0])
c = 0


# In[7]:


getQuadForm(a,b,c, lam = False)


# In[8]:


plt.contour(Xp,Yp,getQuadForm(a,b,c)(Xp,Yp))
plt.gca().set_aspect(1)
plt.grid()
plt.show()


# In[9]:


a = np.array([[-2,1],[0,-1]])
b = np.array([0,0])
c = 0


# In[10]:


getQuadForm(a,b,c, lam = False)


# In[11]:


eigvals(a)


# In[12]:


plt.contour(Xp,Yp,getQuadForm(a,b,c)(Xp,Yp))
plt.gca().set_aspect(1)
plt.grid()
plt.show()


# (chap:parabolexmp)=
# ## Parabolisch

# ### Ein kleiner Ausrutscher in die lineare Algebra: Berechnung des Kern einer Matrix

# Der Kern einer linearen Abbildung $A: \mathbb{R}^n \to \mathbb{R}^m$ ist definiert durch die Menge der Vektoren, welche die lineare Abbildung auf den Nullvektor abbildet. Daher
# 
# $$\mathop{kern} A := \{x \in \mathbb{R}^d | A\cdot x = 0\}$$
# 

# Eine numerische Approximation des Kern einer (quadratischen) Matrix können wir mit Hilfe der zum Eigenwert $\lambda = 0$ zugehörigen Eigenvektoren bestimmen. Das funktioniert solange die Vielfachheit eins ist. Im Fall einer Vielfachheit, ist der Eigenvektor komplexkonjugiert, womit wir mit dem Ansatz nicht alle Basisvektoren für den Nullraum erhalten. In dem Fall müssten man die Hauptvektoren bestimmen.
# 
# Als erstes erstellen wir eine 5x5 Matrix mit Rang 3:

# In[13]:


a = np.random.randint(-9,9,size=(5,3))
# build quadratic matrix with lower rank
a = np.vstack((a.T,a[:,0]-a[:,1],a[:,2]-a[:,1])).T
a


# In[14]:


# Beispiel einer Matrix mit komplex konjugierten Eigenvektoren für den Eigenwert 0.
#a = np.array([[ 3, -7,  6, 10, 13],
#       [ 8,  1,  8,  7,  7],
#       [ 7, -8, -5, 15,  3],
#       [ 2, -7, -8,  9, -1],
#       [-5, -6, -8,  1, -2]])


# Für die Eigenwerte bzw. Eigenvektoren zum Eigenwert 0 folgt

# In[15]:


ew,ev = eig(a)
ns1=ev[:,(np.abs(ew)<1e-13)] # wir sind für das weitere nur am Realanteil interessiert
ns1


# In[16]:


a@ns1


# In[17]:


[norm(a@ni,np.inf) for ni in ns1.T]


# Ein analoger Zugang, jedoch mit der Singulärwertzerlegung liefert immer eine reellwertige Lösung und den gesammten Nullraum:

# In[18]:


atol=1e-13 # absolute Toleranz
rtol=0 # relative Toleranz
u, s, vh = svd(a) # Singulärwerte sind 
tol = max(atol, rtol * s[0])
nnz = (s >= tol).sum()
ns2 = vh[nnz:].conj().T
ns2


# In[19]:


a@ns2


# In[20]:


[norm(a@ni,np.inf) for ni in ns2.T]


# Wir zeigen, dass die Vektoren jeweilen durch Linearkombinationen der anderen darstellbar sind. Womit gezeigt ist, dass die Beschreibungen des Nullraums identisch sind.

# In[21]:


q,r = qr(ns1,mode='economic')
r


# In[22]:


sol=solve_triangular(r,q.T@ns2)
sol


# In[23]:


ns1@sol-ns2


# Der Vorteil der Singulärwertzerlegung sieht man im Resultat der QR-Zerlegung: die beiden Vektoren sind orthogonal. Die $R$ Matrix ist diagonal.

# In[24]:


q,r = qr(ns2,mode='economic')
r


# In[25]:


sol=solve_triangular(r,q.T@ns1)
sol


# In[26]:


ns2@sol-ns1


# Die beiden Beschreibungen sind bis auf numerische Rundung identisch.

# ### Zurück zur quadratischen Form

# Betrachten wir das erwähnte Beispiel: Der Nullraum ist in dem Fall offensichtlich gegeben durch $(0,1)^T$ und entsprechend kompatibel mit dem $b$ Vektor.

# In[27]:


a = np.array([[1,0],[0,0]])
b = np.array([0,1])
c = 0


# In[28]:


getQuadForm(a,b,c, lam = False)


# In[29]:


plt.contour(Xp,Yp,getQuadForm(a,b,c)(Xp,Yp))
plt.gca().set_aspect(1)
plt.grid()
plt.show()


# (chap:hyperexmp)=
# ## Hyperbolisch

# Wir betrachten das Beispiel:

# In[30]:


a = np.array([[1,0],[0,-1]])
b = np.array([0,0])
c = 0


# In[31]:


getQuadForm(a,b,c, lam = False)


# In[32]:


plt.contour(Xp,Yp,getQuadForm(a,b,c)(Xp,Yp))
plt.gca().set_aspect(1)
plt.grid()
plt.show()


# In[ ]:




