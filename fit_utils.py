from scipy.optimize import curve_fit

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class fitClass:

    def __init__(self):
        pass

    def sigmoid(self, t, s):
        return self.p0*np.exp(s*t)/(1 + self.p0*(np.exp(s*t) - 1))
    
    def delay_sigmoid(self,t,s,t_s):
        
        return self.p0*(t<=t_s) + self.p0*np.exp(s*(t - t_s))/(1 + self.p0*(np.exp(s*(t - t_s)) - 1))*(t>t_s)