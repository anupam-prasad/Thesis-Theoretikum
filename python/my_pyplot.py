# proper setup of matplotlib ("import guard")
import matplotlib
try:
    defined_TkAgg
except:
    defined_TkAgg=True
    matplotlib.use('TkAgg') # BEFORE importing pylab (suppresse deprecated use warning)
import matplotlib.pyplot as plt
