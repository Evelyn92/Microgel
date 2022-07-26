import numpy as np
import copy
import math
from utils import dataset_loader,Operations
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator

global surface
surface = ['FOCTS', 'PEG', 'ODS']

def visualize_list_point_clouds(dataset_dict, sample_per_class = 4, colors = ['g', 'b', 'k', 'y', 'm', 'c']):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         PC_list: a list of point clouds
         class_label: class label
         monomer: a string that display monomer type in diagram(NiPAM or NiPMAM)
         colors: list of 6 colors for depicting each row    
    """
    i = -1
    for surf in surface:
        fig = plt.figure(figsize=(25,5))
        i = i+1
        pcs = dataset_dict[surf]
        sample_per_class = sample_per_class if (len(pcs)>sample_per_class) else len(pcs)
        for idx in range(sample_per_class):
            points = pcs[idx]
            ax1 = fig.add_subplot(1, sample_per_class, idx+1, projection="3d")
            mue = np.mean(points, axis=0)
            ax1.scatter(points[:, 0], points[:, 1], points[:, 2], s=1, c=colors[i], marker="s", facecolor="red", lw=0, alpha=1)
            ax1.scatter(mue[0], mue[1], mue[2], s=50, c='r', marker="*", )
            ax1.set_xlabel('X Label')
            ax1.set_ylabel('Y Label')
            ax1.set_zlabel('Z Label')
            ax1.set_title(dataset_dict['selected_polymer']+' Class:{}, #:{}'.format(surf,idx+1))
        plt.show()
    return

def calc_mue(M):
    """
     Args:
         M: point cloud
    returns mean of a point cloud
    """
    return np.mean(M, axis=0)

def calc_dist(p, mue): 
    """    
    At first, minus mue from point clouds, then calculate square 2 of them
     Args:
         p:   point cloud
         mue: mean of point cloud
    
    returns: calculated distance
    """
    return np.linalg.norm(p-mue, 2)

def depict_cmf(dataset_dict, sample_per_class = 4, colors = ['g', 'b', 'k', 'y', 'm', 'c']):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws 
    cumulative distribution function (CDF) all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         PC_list: a list of point clouds
         colors: list of 6 colors for depicting each row    
         class_label: class label
    """
    i = -1
    all_cmf = []
    for surf in surface:
        fig = plt.figure(figsize=(25,5))
        i = i+1
        pcs = dataset_dict[surf]
        sample_per_class = sample_per_class if (len(pcs)>sample_per_class) else len(pcs)
        cmf = []
        for idx in range(sample_per_class):
            points = pcs[idx]
            ax1 = fig.add_subplot(1, 4, idx+1)
            mue = calc_mue(points)
            sorted_dists = sorted([calc_dist(p, mue) for p in points])
            XY = [(0, 0)]
            t = len(sorted_dists)
            for j, d in enumerate(sorted_dists):
                XY.append((d,1.*j/t))
            x, y = [xy[0] for xy in XY], [xy[1] for xy in XY]
            cmf.append(y)
            ax1.plot(x, y, colors[i] ,label='CMF')
            ax1.set_xlabel('Distance')
            ax1.set_ylabel('Probablity')
            ax1.set_title('CMF '+dataset_dict['selected_polymer']+' class:{}, #:{} '.format(surf, idx))
            ax1.legend()		
        plt.show()
        all_cmf.append(cmf)
    return all_cmf

def depict_cmf_group(dataset_dict, surfaces, sample_per_class = 4, colors = ['g', 'b', 'k', 'y', 'm', 'c']):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws 
    cumulative distribution function (CDF) all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         PC_list: a list of point clouds
         colors: list of 6 colors for depicting each row    
         class_label: class label
    """
    i = 0
    group_cmf = []
    fig = plt.figure(figsize=(7,7))
    for surf in surfaces:
        points = dataset_dict[i]
        ax1 = fig.add_subplot(1, 1, 1)
        mue = calc_mue(points)
        sorted_dists = sorted([calc_dist(p, mue) for p in points])
        XY = [(0, 0)]
        t = len(sorted_dists)
        for j, d in enumerate(sorted_dists):
            XY.append((d,1.*j/t))
        x, y = [xy[0] for xy in XY], [xy[1] for xy in XY]
        x = x/np.max(x)
        group_cmf.append(y)
        ax1.plot(x, y, colors[i], label="CMF "+surf)
        i += 1

    plt.legend(loc="upper right")
    ax1.set_xlabel('Distance')
    ax1.set_ylabel('Probablity')
    ax1.set_title('CMF ')
    plt.show()
    return group_cmf

def depict_normalzed_cmf(dataset_dict, sample_per_class = 4, colors = ['g', 'b', 'k', 'y', 'm', 'c']):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws 
    normalized cumulative distribution function (CDF) all point clouds(between 0 and 100).
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         PC_list: a list of point clouds
         colors: list of 6 colors for depicting each row    
         class_label: class label
    """
    i = -1
    all_cmf_normalize = []
    for surf in surface:
        fig = plt.figure(figsize=(25,5))
        i = i+1
        pcs = dataset_dict[surf]
        cmf = []
        sample_per_class = sample_per_class if (len(pcs)>sample_per_class) else len(pcs)
        for idx in range(sample_per_class):
            points = pcs[idx]
            ax1 = fig.add_subplot(1, 4, idx+1)            
            mue = calc_mue(points)
            sorted_dists = sorted([calc_dist(p, mue) for p in points])
            XY = [(0, 0)]
            t = len(sorted_dists)
            for j, d in enumerate(sorted_dists):
                XY.append((d,1.*j/t))
            x, y = [xy[0] for xy in XY], [xy[1] for xy in XY]
            x = x/np.max(x)*100
            cmf.append(y)
            ax1.plot(x, y, colors[i] ,label='CMF')
            ax1.set_xlabel('Distance')
            ax1.set_ylabel('Probablity')
            ax1.set_title('Normalized CMF '+dataset_dict['selected_polymer']+' class:{}, #:{} '.format(surf, idx))
            ax1.legend()
        plt.show()
        all_cmf_normalize.append(cmf)
    return all_cmf_normalize

def depict_cmf_group_aggregation(dataset_dict, colors = ['g', 'b', 'k', 'y', 'm', 'c']):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws 
    cumulative distribution function (CDF) all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         PC_list: a list of point clouds
         colors: list of 6 colors for depicting each row    
         class_label: class label
    """
    i = 0
    fig = plt.figure(figsize=(7,7))
    for surf in surface:
        pcs = dataset_dict[surf]
        points = np.concatenate(pcs, axis= 0)
        Temp = []
        Temp.append(points)
        dataset_dict[surf]  = Temp 
        ax1 = fig.add_subplot(1, 1, 1)
        mue = calc_mue(points)
        sorted_dists = sorted([calc_dist(p, mue) for p in points])
        XY = [(0, 0)]
        t = len(sorted_dists)
        for j, d in enumerate(sorted_dists):
            XY.append((d,1.*j/t))
        x, y = [xy[0] for xy in XY], [xy[1] for xy in XY]
        x = x/np.max(x)
        ax1.plot(x, y, colors[i], label="CMF "+surf)
        i += 1

    plt.legend(loc="upper right")
    ax1.set_xlabel('Distance')
    ax1.set_ylabel('Probablity')
    ax1.set_title('CMF '+dataset_dict['selected_polymer'])
    plt.show()
    return dataset_dict

def depict_pdf(dataset_dict, sample_per_class = 4, colors = ['g', 'b', 'k', 'y', 'm', 'c']):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws 
    probability density function (PDF). Also it depicts mean of each PDF in red.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         PC_list: a list of point clouds
         colors: list of 6 colors for depicting each row    
         class_label: class label
    """
    i = -1
    for surf in surface:
        fig = plt.figure(figsize=(25,5))
        i = i+1
        pcs = dataset_dict[surf]
        sample_per_class = sample_per_class if (len(pcs)>sample_per_class) else len(pcs)
        for idx in range(sample_per_class):
            points = pcs[idx]
            ax1 = fig.add_subplot(1, sample_per_class, idx+1)            
            mue = calc_mue(points)
            sorted_dists = sorted([calc_dist(p, mue) for p in points])
            x_mean = np.mean(sorted_dists)
            XY = [(0, 0)]
            t = len(sorted_dists)
            for j, d in enumerate(sorted_dists):
                XY.append((d,1.*j/t))
            x, y = [xy[0] for xy in XY], [xy[1] for xy in XY]
            x = x/np.max(x)*100
            ax1.hist(sorted_dists, density=True, bins=100)
            ax1.set_title('PDF, class:{}, #:{}'.format(surf, idx))
            ax1.set_xlabel('Distance')
            ax1.set_ylabel('Probablity')
            ax1.legend()
            ax1.axvline(x_mean, color='red', linestyle='dashed', linewidth=1)
        plt.show()
    return
	
def Refine_PC(dataset_dict, sample_per_class = 4, confidence = 0.99, colors = ['g', 'b', 'k', 'y', 'm', 'c'], do_normalize = False):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) , centralized them, removes ouliers, then 
    draws all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors.     
    
     Args:
         PC_list: a list of point clouds
         class_label: class label
         monomer: a string that display monomer type in diagram(NiPAM or NiPMAM)
         confidence: level of confidence for removing points out ot it.
         colors: list of 6 colors for depicting each row  
     Returns:
         processed_points: processed point clouds
    """
    normalized_dataset_dict = copy.deepcopy(dataset_dict)
    i = -1
    for surf in surface:
        fig = plt.figure(figsize=(25,5))
        i   = i+1
        pcs = dataset_dict[surf]
        processed_points = []
        sample_per_class = sample_per_class if (len(pcs)>sample_per_class) else len(pcs)
        for idx in range(len(pcs)):
            points         = pcs[idx]            
            if do_normalize:
                Inlaier_points = []
                mue            = calc_mue(points)
                sorted_dists   = sorted([calc_dist(p, mue) for p in points])
                XY             = [(0, 0)]
                t              = len(sorted_dists)
                for j, d in enumerate(sorted_dists):
                    XY.append(( d, 1.*j/t ))
                    if 1.*j/t >= confidence:
                        break
                sorted_idx     = sorted([(i, calc_dist(p, mue)) for i, p in enumerate(points)], key= lambda d:(d[1], d[0]))
                inlaier        = [d[0] for d in sorted_idx[:j]]
                inlaier_points = [points[j] for j in inlaier]
                points         = np.array(inlaier_points)
                points        -= mue
            processed_points.append(points)
            mue            = np.mean(points, axis=0)
            if(idx<4):
                ax1        = fig.add_subplot(1, sample_per_class, idx+1, projection="3d")
                ax1.scatter(points[:, 0], points[:, 1], points[:, 2], s=1, c=colors[i], marker="*", facecolor="red", lw=0, alpha=1)
                ax1.scatter(mue[0], mue[1], mue[2], s=100, c='r', marker="*")
                ax1.set_xlabel('X Label')
                ax1.set_ylabel('Y Label')
                ax1.set_zlabel('Z Label')
                ax1.set_xlim((-600,600))
                ax1.set_ylim((-600,600))
                ax1.set_zlim((-400,400))
                ax1.set_title('Normalized PC '+dataset_dict['selected_polymer']+' Class:{}, #:{}'.format(surf,idx+1))
        plt.show()
        normalized_dataset_dict[surf]=processed_points
    return normalized_dataset_dict
	
def depict_best_3D(PC_list, sample_per_class = 3, colors = ['g','b','k','y','m','c'], set_axis_limit = False):
    N = len(PC_list)
    idc=0
    scale = 8
    i = -1
    fig = plt.figure(figsize=(25,5))
    i = i+1
    sample_per_class = sample_per_class if ( len(PC_list) > sample_per_class) else len(PC_list)
    for idx in range(sample_per_class):
        points = PC_list[idx]
        X_list=[]
        for i in range(len(points)):
            X_list.append(i)
        Inlaier_points = []
        mue = calc_mue(points)
        ax1 = fig.add_subplot(1, 3, idx+1, projection="3d")
        ax1.scatter(points[:, 0], points[:, 1], points[:, 2], s=1, c=colors[1], marker="s", facecolor="red", lw=0, alpha=1)
        ax1.scatter(mue[0], mue[1], mue[2], s=100, c='r', marker="*")#, markersize=15
        ax1.set_xlabel('X Label')
        ax1.set_ylabel('Y Label')
        ax1.set_zlabel('Z Label')
        ax1.set_title('Class:{}'.format(surface[idx]))
        if(set_axis_limit == True):
            ax1.set_xlim((-600,600))
            ax1.set_ylim((-600,600))
            ax1.set_zlim((-400,400))
            ax1.set_title('Centralized Class:{}'.format(surface[idx]))
        ax1.grid() 
    plt.show()
    return
	
def depict_map_2D(Sel_Points, surface, sample_per_class = 3, colors = ['g','b','k','y','m','c'], set_axis_limit = False, Centralized = False):    
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N           = sample_per_class
    dim         = 0
    Descriptors = []
    while(dim<3):
        fig = plt.figure(figsize=(25,5))
        for idx in range(sample_per_class):
            points = Sel_Points[idx]
            ax1    = fig.add_subplot(1, 3, idx+1)
            if(dim==2):
                x=points[:,0]
                ax1.set_xlabel('X Label')
                y=points[:,1]
                ax1.set_ylabel('Y Label')
            elif(dim==1):
                x=points[:,0]
                ax1.set_xlabel('X Label')
                y=points[:,2]
                ax1.set_ylabel('Z Label')
            elif(dim==0):
                x=points[:,1]
                ax1.set_ylabel('Y Label')
                y=points[:,2]            
                ax1.set_ylabel('Z Label')
            ax1 = fig.add_subplot(1, 3, idx+1)
            
            spacing = 100 # This can be your user specified spacing. 
            minorLocator = MultipleLocator(spacing)
            # Set minor tick locations.
            ax1.yaxis.set_minor_locator(minorLocator)
            ax1.xaxis.set_minor_locator(minorLocator)
            # Set grid to use minor tick locations. 
            ax1.grid(which = 'minor')
            
            ax1.scatter(x, y, s= 2, c=colors[idx], marker="s", facecolor="red", lw=0, alpha=1.)
            if(Centralized):
                ax1.set_xlim((-500,500))
                ax1.set_ylim((-500,500))
                circle2 = plt.Circle((0, 0), 200, color='r', fill=False,linewidth=3)
                ax1.add_patch(circle2)
            ax1.grid(which='minor', alpha=0.2)
            ax1.grid(which='major', alpha=0.8)
            ax1.set_title('Class:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))
        dim=dim+1
    plt.show()
    return 
	
	
def draw_hist_angle_map_2d(hist_angle_out, CLASS_MAP, number_rad, colors = ['g', 'b', 'k','y','m','c']):    
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N = len(hist_angle_out)
    C = len(CLASS_MAP)
    idc=0
    scale = 8
    X_list=[]
    for i in range(number_rad):
        X_list.append(i)
    dim=0
    while(idc<3):
        fig = plt.figure(figsize=(25,5))
        for idx in range(3):
            points = hist_angle_out[idc*3+idx]
            ax1    = fig.add_subplot(1, 3, idx+1)
            ax1.plot(X_list, points/max(points), colors[idc] ,label='PDF')
            ax1.set_title('Density for :{} sector , class:{}, {}'.format(number_rad,CLASS_MAP[idx],CLASS_Plane[idc]))
            ax1.set_xlabel('Distance')
            ax1.set_ylabel('Probablity')
            ax1.set_xlim((0,45))
            #ax1.set_ylim((0,0.15))
            ax1.grid() 
            ax1.legend()
        idc=idc+1
    plt.show()
    return

def depict_pointcloud_2D(Sel_Points, CLASS_MAP, Max_Point_Sector, Centralized, number_bin, radius=200):
    CLASS_Plane=['YZ Plane','XZ Plane','XY Plane']  
    colors = ['g', 'b', 'k','y','m','c']
    N=len(Sel_Points)
    dim=0
    Descriptors = []
    while(dim<3):
        fig = plt.figure(figsize=(25,5))
        for idx in range(len(Sel_Points)):
            points=Sel_Points[idx]
            ax1 = fig.add_subplot(1, len(Sel_Points), idx+1)
            if(dim==2):
                x=points[:,0]
                ax1.set_xlabel('X Label')
                y=points[:,1]
                ax1.set_ylabel('Y Label')
            elif(dim==1):
                x=points[:,0]
                ax1.set_xlabel('X Label')
                y=points[:,2]
                ax1.set_ylabel('Z Label')
            elif(dim==0):
                x=points[:,1]
                ax1.set_ylabel('Y Label')
                y=points[:,2]            
                ax1.set_ylabel('Z Label')          
            spacing = 100 # This can be your user specified spacing. 
            minorLocator = MultipleLocator(spacing)
            ax1.xaxis.set_minor_locator(minorLocator)
            ax1.yaxis.set_minor_locator(minorLocator)
            ax1.scatter(x, y, s= 2, c=colors[idx], marker="s", facecolor="red", lw=0, alpha=1.)
            if(Centralized):
                ax1.set_xlim((-400,400))
                ax1.set_ylim((-400,400))
                circle2 = plt.Circle((0, 0), radius, color='r', fill=False,linewidth=3)
                ax1.add_patch(circle2)
                v_origin = [0,0]
                v = [350,0]
                for angle in range(number_bin):
                    if(angle!=0):
                        theta=360/number_bin
                        xp=v[0]*math.cos(math.radians(theta))-v[1]*math.sin(math.radians(theta))
                        yp=v[1]*math.cos(math.radians(theta))+v[0]*math.sin(math.radians(theta))
                        v=[xp,yp]
                    angular_uncert = 20.
                    v_angle = np.arctan2(v[1],v[0])
                    ax1.arrow(v_origin[0], v_origin[1] ,v[0], v[1], head_width=10,
                             head_length=10, lw=2, fc='#777777', ec='#777777', length_includes_head=True)
            #ax1.grid(which='minor', alpha=0.2)
            #ax1.grid(which='major', alpha=0.8)
            ax1.set_title('Class:{}, #Dimension:{}, Plane:{}'.format(CLASS_MAP[idx], dim,CLASS_Plane[dim]))
        dim=dim+1
    plt.show()
    return
	
	
def draw_histogram_pointcloud_plane(train_points, CLASS_MAP, radius = 200, colors = ['g', 'b', 'k','y','m','c']):
    processed_points = []
    N=len(train_points)
    C = len(CLASS_MAP)
    idc=0
    Descriptors = []
    while(idc<3):
        fig = plt.figure(figsize=(25,5))
        for idx in range(len(train_points)):
            points = train_points[idx]
            Inlaier_points = []
            mue = calc_mue(points)
            sorted_dists = sorted([calc_dist(p, mue) for p in points])
            sorted_idx   = sorted([(i, calc_dist(p, mue)) for i, p in enumerate(points)], key= lambda d:(d[1], d[0]))
            XY = [(0, 0)]
            t = len(sorted_dists)
            for i, d in enumerate(sorted_dists):
                XY.append(  ( d, 1.*i/t )  )
                if 1.*i/t >= 0.99:
                    break
            inlaier = [d[0] for d in sorted_idx[:i]]
            inlaier_points = [points[i] for i in inlaier]
            points = np.array(inlaier_points)
            points -= mue
            mue = np.mean(points, axis=0)
            points = Operations.map_3D(points, dim=idc)
            Descriptors.append(points)             
            ax1 = fig.add_subplot(1, 4, idx+1)
            ax1.scatter(range(radius), points, s=5, c=colors[1], marker="s", facecolor="red", lw=0, alpha=1.)
            ax1.set_xlabel('X Label')
            ax1.set_ylabel('Y Label')
            #ax1.set_ylim((0,0.05))
            ax1.grid()
            ax1.set_title('Class:{}, #Dimension:{}'.format(CLASS_MAP[idx], idc))
        idc=idc+1
    plt.show()
    Descriptors = np.array(Descriptors)
    return