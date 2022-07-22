import numpy as np

class MakeDataset():
    
    def __init__(self, tsp_path, with_speed_limits = False):
        self.with_speed_limits = with_speed_limits
        self.tsp_path = tsp_path
    
    def get_dist_matrix(self):
        m = []
        with open(self.tsp_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                curLine = []
                for cur_dist in line.split(' '):
                    curLine.append(cur_dist)
                m.append(curLine)
        return np.array(m)
    
    def get_speed_limits(self, distance_martix):
        cnt = distance_martix.shape[0]
        v = np.random.rand(cnt)
        return v

    def get_data(self):
        if self.with_speed_limits:
            distance_matrix = self.get_dist_matrix()
            speed_limits = self.get_speed_limits(distance_matrix)
            return distance_matrix, speed_limits
        else:
            return self.get_dist_matrix()