import numpy as np
import pandas as pd
from TRSF.homology.cripser_homol import (compute_ph_cripser, 
        apply_confusion_limit, classify_point, get_enclosing_mask, death_correc)
import pandas as pd

def test_compute_ph_cripser():
    # float image with one peak in the middle
    img = np.array([
        [0.01, 0.01, 0.01, 0.01, 0.01],
        [0.01, 1.5, 1.5, 1.5, 0.01],
        [0.01, 1.5, 2.0, 1.5, 0.01],
        [0.01, 1.5, 1.5, 1.5, 0.01],
        [0.01, 0.01, 0.01, 0.01, 0.01]
    ])

    local_bg = 0.01
    sigma = 5
    maxdim = 0

    result = compute_ph_cripser(img, local_bg, sigma, maxdim)
    print(result)

    assert isinstance(result, pd.DataFrame)

    expected_columns = ['Birth','Death','x1','y1','x2','y2','lifetime']
    assert list(result.columns) == expected_columns

    expected_values = [
        [2.0, 0.05, 2.0, 2.0, 0.0, 0.0, -1.95]]
    
    assert result.values.tolist() == expected_values



def test_apply_confusion_limit():
    # Create a persistence diagram DataFrame
    data = {
        'Birth': [-1.0, -1.0, -1.0, -1.0, -0.5, -0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0],
        'Death': [-0.5, -0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5],
        'x1': [0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2],
        'y1': [0, 0, 1, 3, 1, 0, 1, 3, 1, 0, 1, 3],
        'x2': [0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2],
        'y2': [1, 1, 2, 3, 2, 1, 2, 3, 2, 1, 2, 3],
        'lifetime': [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
    }
    p = pd.DataFrame(data)

    # Set the confusion limit
    confusion_limit = 2.0

    # Call the function to apply the confusion limit
    result = apply_confusion_limit(p, confusion_limit)

    # Check if the result is a pandas DataFrame
    assert isinstance(result, pd.DataFrame)

    # Check the column names of the DataFrame
    expected_columns = ['Birth', 'Death', 'x1', 'y1', 'x2', 'y2', 'lifetime']
    assert list(result.columns) == expected_columns

    # Check the values in the DataFrame
    expected_values = [[0.0,0.5,2.0,3.0,2.0,3.0,0.0]]
    print(result.values.tolist())
    assert result.values.tolist() == expected_values

import pandas as pd

def test_classify_point():
    # Create a row DataFrame
    result = pd.DataFrame({
        'parent_tag': [1, 2, 3],
        'name': [1, 2, 3],
        'len_enclosed': [2, 1, 2],
        'new_row': [0, 0, 1]
    })

    # Call the function to classify the point
    result['class'] = result.apply(lambda row: classify_point(row),axis=1)
    # Check the result
    expected_result = [4, 0, 5]
    assert result['class'].tolist() == expected_result



def test_get_enclosing_mask():
    # Create a mask
    mask = np.array([[False, True, False, False, True],
                     [True, True, True, False, True],
                     [False, True, False,False, True]], dtype=bool)

    # Call the function to get the enclosing mask for point (1, 1)
    result = get_enclosing_mask(1, 1, mask)

    # Check if the result is a numpy array
    assert isinstance(result, np.ndarray)

    # Check the shape of the result
    assert result.shape == mask.shape

    # Check the values in the result
    expected_result = np.array([[False, True, False, False, False],
                                [True, True, True, False, False],
                                [False, True, False,False,False]], dtype=bool)

    assert np.array_equal(result, expected_result)


