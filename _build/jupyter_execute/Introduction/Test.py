#!/usr/bin/env python
# coding: utf-8

# # Test

# In[1]:


from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.webgui import Draw


# In[2]:


mesh = Mesh(unit_square.GenerateMesh(maxh=1/3))


# In[3]:


Draw(mesh);


# In[ ]:




