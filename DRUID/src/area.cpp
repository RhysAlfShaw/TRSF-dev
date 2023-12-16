//     for i in tqdm(range(0,len(pd)),total=len(pd),desc='Calculating area',disable=not output):
//             area = calculate_area(pd.iloc[i],img)
//             #print(area)
//             area_list.append(area)

#include <vector>

int sumMatrixElements(const std::vector<std::vector<int>>& matrix) {
    int sum = 0;
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            sum += element;
        }
    }
    return sum;
}

// function that takes value x,y, birth and death and a image 

std::vector<std::vector<bool>> get_mask(const std::vector<std::vector<float>>& img, int birth, int death) {
    
    std::vector<std::vector<bool>> mask(img.size(), std::vector<bool>(img[0].size(), false));

    for (size_t i = 0; i < img.size(); ++i) {
        for (size_t j = 0; j < img[i].size(); ++j) {
            if (img[i][j] <= birth && img[i][j] > death) {
                mask[i][j] = true;
            }
        }
    }

    return mask;
}
