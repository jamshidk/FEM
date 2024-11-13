#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:02:13 2024

@author: jamkavosi
"""

import numpy as np
import math
from functions import *

x=np.array([[0,0],
           [1,0],
           [1,2],
           [0,2]])

GPE=1 #4

K=stiffness(x, GPE)



