import os
import glob
import trimesh
import numpy as np
from matplotlib import pyplot as plt
# import plotly.express as px
import pandas as pd
import math
 

def get_dataset(DATA_DIR, temperature = 24, surface = ['FOCTS', 'PEG', 'ODS', 'GLAS'], load_NIPAM = True, Flag_zreverse= True):
    available_temperature = os.listdir(DATA_DIR)
    if (str(temperature) not in available_temperature):
        raise Exception("input temperature has not any proper data")
    list_surface = []
    available_surface = os.listdir(DATA_DIR+'/'+str(temperature))
    for surf in available_surface:
        if surf in surface:
            list_surface.append(surf)
    selected_polymer = "NiPAM" if load_NIPAM == True else "NiPMAM"
    dataset_dict = {'Temperature':temperature , 'selected_polymer':selected_polymer}
    for surf in list_surface:
        PC_files   = glob.glob(os.path.join(DATA_DIR, str(temperature), surf, selected_polymer, "*.ply"))
        point_list = []
        for file in PC_files:            
            temp      = trimesh.load(file)[:]
            if Flag_zreverse:
                temp[:,2] = 1-temp[:,2]
            point_list.append(temp)
        dataset_dict[surf] = point_list
    return dataset_dict
