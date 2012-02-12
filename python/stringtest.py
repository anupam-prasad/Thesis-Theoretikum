#!/usr/bin/env python
import numpy as np
from math import *

str1='(hello world) 1'
print str1.split('(')[0]
print str1.split('(')[0].strip(')').split(',')

