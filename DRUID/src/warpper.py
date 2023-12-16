import ctypes

area_lib = ctypes.CDLL('./DRUID/src/area.o')

# specify the function return type it returns std::vector<std::vector<bool>>
area_lib.get_mask.restype = ctypes.vector(ctypes.c_bool)


import numpy as np

test_image = np.random.rand(0, 255, (100, 100), dtype=np.float64)

mask = area_lib.get_mask(test_image, 10, 5)
print(mask)