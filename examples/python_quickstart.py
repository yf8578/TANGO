import numpy as np

from tangoqtl import tango_test

z = np.array([2.1, 1.7, -0.4, 0.8, 1.2])
corr = np.eye(len(z))

result = tango_test(z, corr=corr)
print(result.as_dict())
