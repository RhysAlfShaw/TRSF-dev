import numpy as np
cimport numpy as np

cpdef np.ndarray[bint, ndim=2] compute(np.ndarray[bint, ndim=2] non_overlapping,
                                       np.ndarray[float, ndim=2] Image,
                                       np.ndarray[bint, ndim=2] componet,
                                       float min_val):
    cdef int i, j
    cdef double prop_pixel_val, mean_without
    cdef list current_neightbour_vals

    for i in range(non_overlapping.shape[0]):
        for j in range(non_overlapping.shape[1]):
            if not non_overlapping[i, j]:
                continue

            temp_img = Image * componet
            prop_pixel_val = Image[i, j]
            current_neightbour_vals = []

            # left pixel
            try:
                current_neightbour_vals.append(temp_img[i, j - 1])
            except:
                continue
            # right pixel
            try:
                current_neightbour_vals.append(temp_img[i, j + 1])
            except:
                continue
            # top pixel
            try:
                current_neightbour_vals.append(temp_img[i - 1, j])
            except:
                continue
            # bottom pixel
            try:
                current_neightbour_vals.append(temp_img[i + 1, j])
            except:
                continue
            # top left pixel
            try:
                current_neightbour_vals.append(temp_img[i - 1, j - 1])
            except:
                continue
            # top right pixel
            try:
                current_neightbour_vals.append(temp_img[i - 1, j + 1])
            except:
                continue
            # bottom left pixel
            try:
                current_neightbour_vals.append(temp_img[i + 1, j - 1])
            except:
                continue
            # bottom right pixel
            try:
                current_neightbour_vals.append(temp_img[i + 1, j + 1])
            except:
                continue

            # calculate the mean of the current_neightbour_vals
            current_neightbour_vals = [x for x in current_neightbour_vals if x != 0]
            mean_without = np.mean(current_neightbour_vals)
            current_neightbour_vals.append(prop_pixel_val)

            if (mean_without >= prop_pixel_val) & (prop_pixel_val >= min_val):
                non_overlapping[i, j] = True
            else:
                non_overlapping[i, j] = False

    return non_overlapping
