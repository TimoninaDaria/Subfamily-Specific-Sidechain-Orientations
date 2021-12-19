import numpy as np

def find_elbow(data, theta):
    # make rotation matrix
    co = np.cos(theta)
    si = np.sin(theta)
    rotation_matrix = np.array(((co, -si), (si, co)))

    # rotate data vector
    rotated_vector = data.dot(rotation_matrix)
    # return index of elbow
    return np.where(rotated_vector == rotated_vector.min())[0][0]

def get_data_radiant(data):
    return np.arctan2(data[:, 1].max() - data[:, 1].min(), 
                    data[:, 0].max() - data[:, 0].min())



