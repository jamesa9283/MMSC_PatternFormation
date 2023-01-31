#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 17:45:24 2023

@author: arthur
"""

from matplotlib import pyplot as plt
import itertools
import numpy as np
from numpy import arange
from numpy import meshgrid

a = 0.1
b = 0.9

delta = 0.025
urange = arange(0, 5.0, delta)
vrange = arange(0, 7.0, delta)
U, V = meshgrid(urange,vrange)

# F is for the u_t equation. G for the v_t
F = a - U + U**2*V
G = b - U**2*V

# plt1 = plt.contour(U, V, F, [0], colors="r", label="F")
# plt2 = plt.contour(U, V, G, [0], colors="b", label="G")
# plt.plot([1], [0.9], marker="o", markersize=7, markeredgecolor="red", markerfacecolor="green")
# plt.text(1.1,1.1,'(1, 0.9)')
# plt.title("Nullclines of Schnakenberg Model")
# plt.legend([plt1, plt2], ["F", "G"], loc="upper left")
# plt.show()

fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)
cs1 = axes.contour(U, V, F, [0], colors="r", label="F")
cs2 = axes.contour(U, V, G, [0], colors="b", label="G")
axes.plot([1], [0.9], marker="o", markersize=7, markeredgecolor="red", markerfacecolor="green")
proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_edgecolor()[0]) for pc in itertools.chain(cs1.collections, cs2.collections)]
plt.legend(proxy, ["F", "G"])
plt.show()
