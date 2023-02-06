import numpy as np
import math
import cv2
import pandas as pd
import mahotas
import struct
from utils import dataset_loader, my_utils, Operations, Moments
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage.measure import regionprops
from skimage.transform import rescale, resize
from skimage.morphology import thin
from skimage.morphology import skeletonize
from skimage.morphology import medial_axis
from matplotlib.ticker import MultipleLocator

def generate_gaussian(sigma_x, sigma_y):
    """
    generate gaussian distribution and sampling from this distribution
    
    Parameters
    ----------
    sigma_x : float
        sigma distribution in x dimension
    sigma_y : float
        sigma distribution in x dimension

    Returns
    -------
    xyz_points: 2d numpy array
        a numpy array with 3 columns x, y, z
    """
    
    x    = np.linspace(-3, 3, num = 100)
    y    = np.linspace(-3, 3, num = 100)
    x, y = np.meshgrid(x, y)
    z    = np.exp(-sigma_x*x**2-sigma_y*y**2)
    width, height = x.shape
    x     = 100*x
    x_new = np.asarray(x)
    x_new = x_new.reshape(width*height)
    y     = 100*y
    y_new = np.asarray(y)
    y_new = y_new.reshape(width*height)
    bias  = 100 * np.ones((width, height))
    z     = 300*z
    z_new = np.asarray(z)
    z_new = z_new.reshape(width*height)
    x     = x_new.reshape(width,height)
    y     = y_new.reshape(width,height)
    z     = z_new.reshape(width,height)
    fig   = plt.figure()
    ax    = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x,y,z)#, cmap=cm.jet
    plt.show()
    xyz_points = np.vstack([x_new, y_new, z_new]).T
    xyz_points.shape
    return xyz_points
	
def write_pointcloud(filename, point_cloud, rgb_points = None):
    """ 
    Write header of .ply file
    assert point_cloud.shape[1] == 3,''
    assert point_cloud.shape == rgb_points.shape,'Input RGB colors should be Nx3 float array and have same size as input XYZ points'

    Parameters
    ----------
    filename : string
        path to write ply point_cloud
    point_cloud : 2d numpy array
        Input point cloud should be Nx3 float array
    rgb_points : 
        to add color to point cloud data

    Returns
    -------
    """
    
    if rgb_points is None:
        rgb_points = np.ones(point_cloud.shape).astype(np.uint8)*255
    fid = open(filename,'wb')
    fid.write(bytes('ply\n', 'utf-8'))
    fid.write(bytes('format binary_little_endian 1.0\n', 'utf-8'))
    fid.write(bytes('element vertex %d\n'%point_cloud.shape[0], 'utf-8'))
    fid.write(bytes('property float x\n', 'utf-8'))
    fid.write(bytes('property float y\n', 'utf-8'))
    fid.write(bytes('property float z\n', 'utf-8'))
    fid.write(bytes('property uchar red\n', 'utf-8'))
    fid.write(bytes('property uchar green\n', 'utf-8'))
    fid.write(bytes('property uchar blue\n', 'utf-8'))
    fid.write(bytes('end_header\n', 'utf-8'))
    for i in range(point_cloud.shape[0]):
        fid.write(bytearray(struct.pack("fffccc",point_cloud[i,0],point_cloud[i,1],point_cloud[i,2],
                                        rgb_points[i,0].tostring(),rgb_points[i,1].tostring(),
                                        rgb_points[i,2].tostring())))
    fid.close()
    return
