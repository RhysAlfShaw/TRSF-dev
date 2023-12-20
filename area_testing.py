import numpy as np
import pandas 
import matplotlib.pyplot as plt
import time
from multiprocess import Pool

def calculate_area(row, img):

    mask = get_mask(row,img)
    
   #t0 = time.time()
    area = np.sum(mask)
    #t1_area = time.time()
    #print('Time to get area: ',t1_area-t0)
    return area


def get_enclosing_mask(x, y, mask):
    '''
    Returns the mask of the enclosed area of the point (x,y) in the mask.
    '''
    
    # Ensure the point is inside the mask
    #print(mask[y, x])
    if not mask[y, x]:
        return None
    # Create a copy of the mask
    enclosed_mask = np.copy(mask)
    # Perform a flood fill starting from the point
    h, w = mask.shape
    stack = [(x, y)]
    while stack:
        x, y = stack.pop()
        if enclosed_mask[y, x]:
            enclosed_mask[y, x] = False
            if x > 0:
                stack.append((x - 1, y))
            if x < w - 1:
                stack.append((x + 1, y))
            if y > 0:
                stack.append((x, y - 1))
            if y < h - 1:
                stack.append((x, y + 1))
    
    return mask & (~enclosed_mask)

def process(i):
    return calculate_area(pd.iloc[i], img)

def get_mask(row, img):
    mask = np.zeros(img.shape)
    birth = row.Birth
    death = row.Death
    
    #t0 = time.time()
    # bottleneck
    # use sparse matrix?
    
    np.logical_and(img <= birth, img > death, out=mask)
    mask = np.where(mask, 1, 0)
    #mask = np.logical_or(mask,np.logical_and(img <= birth,img > death))
    #t1_mask_expr = time.time()
    #print('Time to get mask: ',t1_mask_expr-t0)

    #t0 = time.time()
    mask_enclosed = get_enclosing_mask(int(row.y1),int(row.x1),mask)
    #print(mask_enclosed)
    #t1_getenclosed = time.time()
    #print('Time to get enclosed: ',t1_getenclosed-t0)
    return mask_enclosed


pd = pandas.read_csv('test_pd_area.csv')
print(pd.head())
img= np.load('are_test_image.npy')


i = 8

t0 = time.time()
row = pd.iloc[i]
#t1_row = time.time()
#print('Time to get row: ',t1_row-t0)
t0 = time.time()

with Pool() as p:
    area_list = list(p.imap(process, range(len(pd))))
t1_area_mp = time.time()
print('Time to get area Parrellel: ',t1_area_mp-t0)

# tranditionall
t0 = time.time()
area_list = []
for i in range(len(pd)):
    area_list.append(calculate_area(pd.iloc[i], img))
t1_area = time.time()
print('Time to get area: ',t1_area-t0)

print("AREA :" ,area_list)



