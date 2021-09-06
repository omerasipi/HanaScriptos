#!/usr/bin/env python
# coding: utf-8

# # Variable Koeffizienten

# Das folgende Beispiel folgt dem analogen aus dem interaktiven Kurs von Joachim Schöberl {cite}`schoeberliFEM`.$\DeclareMathOperator{\opdiv}{div}$ $\DeclareMathOperator{\setR}{R}$

# In[1]:


from ngsolve import *
from ngsolve.webgui import Draw


# ## Problem
# 
# Ein Setup mit unterschiedlicher Wärmeleitungskoeffizienten wir mit Hilfe der Gleichung
# 
# $$
# -\opdiv \lambda(x) \nabla u(x) = f(x)
# $$
# 
# modelliert, wobei $\lambda=\lambda(x)$ der ortsabhängige Wärmeleitungskoeffizient ist. Der zugehörige Wärmefluss  
# 
# $$
# q = -\lambda \nabla u
# $$
# 
# ist gegeben durch den Gradient der Temperatur $\nabla u$. Im Fall, dass $\lambda$ diskontinuierlich ist, wird die Gleichung im Sinne von Distributionen verstanden. Dies beinhaltet Interface Bedingungen: Die Temperatur auf der linken und rechten Seite sind gleich und der Wärmefluss der linken in die rechte Seite müssen gleich gross sein:  
# 
# $$\begin{split}
# u_l & = & u_r \\
# \lambda_l \frac{\partial u_r}{\partial n} & = & \lambda_r \frac{\partial u_r}{\partial n}
# \end{split}$$
# 
# Die variationelle Form des Problems ist gegeben durch: finde $u \in H^1(\Omega)$ mit
# 
# $$
# \int_\Omega \lambda(x) \nabla u \nabla v dx = \int_\Omega f v dx
# $$
# 
# Diskontinuierliche Koeffizienten bilden kein Problem. Beide Interface Bedingungen werden erfüllt:
# * Stetigkeit der Temperatur $u$ durch die Stetigkeit der Ansatzfunktionenraumes
# * Stetigkeit des Wärmeflusses in einem schwachen Sinne, durch Neumann Randbedingungen.

# ## Geometrie und Mesh

# In[2]:


from netgen.geom2d import *
geo = SplineGeometry()
geo.AddRectangle( (0,0), (1,1), leftdomain=1, rightdomain=0, 
                 bcs=['b','r','t','l'])
geo.AddCircle( (0.3,0.7), 0.1, leftdomain=2, rightdomain=1)
geo.AddRectangle ( (0.2,0.2), (0.9,0.3), leftdomain=3, rightdomain=1)
geo.SetMaterial(1, "air")
geo.SetMaterial(2, "source")
geo.SetMaterial(3, "bar")

mesh = Mesh(geo.GenerateMesh(maxh=0.03))
mesh.Curve(3)
Draw (mesh,)


# Die Subdomains (Materials) sind gegeben durch:

# In[3]:


print (mesh.GetMaterials())


# Die Geometrie hat folgende Ränder (innere und äussere):

# In[4]:


print (mesh.GetBoundaries())


# In[5]:


V = H1(mesh, order=3, dirichlet="b|r")
u = V.TrialFunction()
v = V.TestFunction()


# Wir definieren nun das Material in dem wir den Wärmekoeffizient stückweise konstant ansetzen:

# In[6]:


lamvalues = { "air" : 1, "bar" : 1e2, "source" : 2 }
lam = CoefficientFunction( 
    [lamvalues[mat] for mat in mesh.GetMaterials()])
Draw (log(lam), mesh, "log lambda")


# Wir nützen die Wärmeleitfähigkeit für die Definition der Bilinearform:

# In[7]:


a = BilinearForm(V)
a += lam*grad(u)*grad(v)*dx

f = LinearForm(V)
f += 1*v*dx("source")
a.Assemble()
f.Assemble()


# und lösen das System:

# In[8]:


gfu = GridFunction(V)
gfu.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec


# Für die Wärmeverteilung erhalten wir:

# In[9]:


Draw (gfu, mesh, "temperature")


# Für den Gradient der Wärmeverteilung $\nabla u$

# In[10]:


Draw (grad(gfu), mesh, "gradient",vectors=grad(gfu))


# und für den Wärmefluss $-\lambda \nabla u$:

# In[11]:


Draw (-lam*grad(gfu), mesh, "heatflux")


# In[ ]:




