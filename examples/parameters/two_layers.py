# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     notebook_metadata_filter: all,-toc,-latex_envs
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.5.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.7.6
# ---

# %% [markdown]
# <img src="media/image1.jpeg" width="30%">

# %% [markdown]
# <img src="media/image2.png" width="30%">

# %% [markdown]
# <img src="media/image6.png" width="30%">

# %% [markdown]
# <img src="media/image8.png" width="30%">

# %% [markdown]
# <img src="media/image13.jpeg" width="30%">

# %%
"""
   layers.py provides two functions:
      do_two(Sol=341.,epsilon=0.55,albedo=0.3)
         -- this solves the two layer model analytically

      do_two_matrix(Sol=341.,albedo=0.3,epsilon1=0.55,epsilon2=0.55)         
          -- this solves the two layer model numerically
"""

# %%
from __future__ import division
from __future__ import print_function
import numpy as np
from scipy import linalg

# %%
#show the difference between importing and running this
#module by printing the module namespace
#
print("my namespace is: ", __name__)

# %%
sigma=5.67e-8  #W/m^2/K^4


# %%
def do_two(Sol=341.,epsilon=0.55,albedo=0.3):
    """
       do_two(Sol=341.,epsilon=0.55,albedo=0.3)
       returns [Fg,F1,F2]   -- layer fluxes in W/m^2
    """
    Sol=341.
    epsilon=0.55
    Fg=(1-albedo)*Sol/(1-epsilon/2*(1+epsilon*(1-epsilon)/2)/ \
                         (1-epsilon*epsilon/4)-(1-epsilon)*epsilon/2*((1-epsilon)+ \
                        epsilon/2)/(1-epsilon*epsilon/4))
    F2=Fg*epsilon/2*((1-epsilon)+epsilon/2)/(1-epsilon*epsilon/4)
    F1=Fg*epsilon/2*(1+epsilon*(1-epsilon)/2)/(1-epsilon*epsilon/4)
    #check balances
    TOA= Sol*(1 - albedo) - F2  - (1-epsilon1)*F1 - (1-epsilon1)*(1-epsilon2)*Fg
    Lay1=Sol*(1 - albedo) + F2 - F1 - (1 - epsilon1)*Fg
    Ground=Sol*(1 - albedo) + F1 + (1-epsilon1)*F2 - Fg
    return (Fg,F1,F2)


# %%
def do_two_matrix(Sol=341.,albedo=0.3,epsilon1=0.55,epsilon2=0.55):
    """
       do_two_matrix(Sol=341.,albedo=0.3,epsilon1=0.55,epsilon2=0.55)
       returns [Fg,F1,F2]   -- layer fluxes in W/m^2
    """   
    Sol=Sol*(1-albedo)
    abs1=epsilon1
    abs2=epsilon2
    Tr1=1. - abs1
    Tr2=1. - abs2
    #layer 2 budget
    #dF2/dt = abs2*Tr1*Fg + abs2*F1 - 2*F2
    #layer 1 budget
    #dF1/dt = abs1*Fg - 2*F1 + abs1*F2
    #Ground budget
    #dFg/dt = Sol - Fg + F1 + Tr1*F2
    the_array=[[abs2*Tr1, abs2, -2.], \
               [abs1, -2., abs1],\
               [-1., 1., Tr1]]
    the_array=np.array(the_array)
    rhs=[0,0,-Sol]
    the_inv=linalg.inv(the_array)
    fluxes=np.dot(the_inv,rhs)
    return fluxes


# %%
def find_temps(fluxes,epsilon1=0.55,epsilon2=0.55):
    """
       find_temps(fluxes,epsilon1=0.55,epsilon2=0.55)
       fluxes=(Fg,F1,F2)
       returns (Tg,T1,T2)
    """
    Tg=(fluxes[0]/sigma)**0.25
    T1=(fluxes[1]/(sigma*epsilon1))**0.25
    T2=(fluxes[2]/(sigma*epsilon2))**0.25
    return (Tg,T1,T2)   


# %% [markdown]
# the following boiler plate runs test code when
# layers.py is executed from the command line and it's namespace
# is "__main__", but not when it is imported
# by other code and it's namespace is "layers"

# %%
if __name__=="__main__":
    #analytic_fluxes=do_two()
    numeric_fluxes=do_two_matrix(epsilon1=0.6,epsilon2=0.4,Sol=240.,albedo=0.)
    #print("analytic temperatures: ",find_temps(analytic_fluxes))
    print("numeric temperatues: ",find_temps(numeric_fluxes))
    print("numeric fluxes: ",numeric_fluxes)

# %%
