import  numpy as np
import math
from matplotlib import pyplot as plt
from utils import my_utils
from matplotlib.ticker import MultipleLocator

def find_thresholds(arr, zoom):
    arr_min = np.min(arr)
    arr_max = np.max(arr)

    data_range = arr_max - arr_min
    zoom_bin = data_range / zoom

    if zoom%2 == 1:
        min_threshold = arr_min + (zoom_bin * int(zoom/2))
        max_threshold = arr_max - (zoom_bin * int(zoom/2))
    else:
        min_threshold = arr_min + (zoom_bin * (int(zoom/2)-0.5))
        max_threshold = arr_max - (zoom_bin * (int(zoom/2)-0.5))

    arr = np.ma.array(arr)
    arr_masked = \
        np.ma.masked_where(\
                (np.logical_xor((min_threshold < arr), (arr < max_threshold))), \
                arr, np.nan)

    return arr_masked

def zoom_function(xyz_values, zoom = 10):
    x_values_masked = find_thresholds(xyz_values[:, 0], zoom)
    y_values_masked = find_thresholds(xyz_values[:, 1], zoom)
    
    return x_values_masked, y_values_masked

def quatize_distribution(x, y, bin = 100):
    Features = []
    Features.append(0)
    for idx in range(bin-1):
        j = 0
        while x[j]<=idx+1:
            j += 1    
        Features.append(y[j])
    return np.arange(bin), np.array(Features)

def draw_histogram_angle(hist_angle, number_bin = 30, sample_per_class = 4, class_name = 'FOCTS', colors = ['g','b','k','y','m','c']):
    """
    This function receives a list of histograms point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         hist_angle: a list of point clouds histogram
         class_label: a list of strings
         number_bin: number of sectors
         colors: list of 6 colors for depicting each row    
    """
    X_list=[]
    for i in range(number_bin):
        X_list.append(i) 
    i = -1
    fig = plt.figure(figsize=(25,5))
    i = i+1
    pcs = hist_angle
    sample_per_class = sample_per_class if (len(pcs)>sample_per_class) else len(pcs)
    for idx in range(sample_per_class):
        points = pcs[idx]
        ax1 = fig.add_subplot(1, sample_per_class, idx+1)
        mue = np.mean(points, axis=0)
        ax1.plot(X_list, points/max(points), colors[0] ,label='PDF')
        ax1.set_title('Density Histogram Angle, class:{}, #:{}'.format(class_name, idx+1))
        ax1.set_xlabel('Number Histogram Angle')
        ax1.set_ylabel('Probablity')
        #ax1.set_ylim((0,0.05))
        ax1.grid() 
        ax1.legend()
    plt.show()
    return

def histogram_angle(PC_list, number_bin):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and number of bins.
    divide space to several sectors and calculate number of points in each sector.
	
    for calculating an angle it's necessary to have three points: "world_center,fixed_point" are reference points and "candidate_point" is candidate point.
    world_center is center of coordinate system.
    fixed_point is an arbitrary point in x axis, for example with unit length
    now having three points it's easy to calculate angle of each point respect to "world_center,fixed_point"
    np.arccos(Theta): returns angle between 0 and 180, but I should calculate angle between 0 and 360. So considering "y"
    element of "candidate_point" could be useful.
	
	           a
	         -|
	       -|
	     -|
	    b-------------c
    
    
     Args:
         PC_list(list): a list of point clouds
         number_bin(Integer): number of sectors
     Returns:
         hist_angle(list): list of calculated histogram         
    """
    world_center = np.array([0,0,0])
    fixed_point = np.array([1,0,0])
    bc = fixed_point - world_center
    hist_angle=[]
    for i in range(len(PC_list)):
        cur_hist=[0]*number_bin
        points_normalize=PC_list[i]
        for j in range(len(points_normalize)):
            candidate_point=points_normalize[j]
            ba = candidate_point - world_center
            cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            angle = np.arccos(cosine_angle)
            if(ba[1]<0):
                angle=360-np.degrees(angle)
            else:
                angle=np.degrees(angle)
            which_bin=math.floor(angle/(360/number_bin))
            cur_hist[which_bin]=cur_hist[which_bin]+1
        cur_hist=np.array(cur_hist)/len(points_normalize)
        hist_angle.append(cur_hist)
    return hist_angle

def compare_histogram_angle(PC_list, hist_angle, class_name = 'ODS', number_bin = 30, visualize = False):
    """
    This function receives a list of histograms point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws mod point cloud 
    and it's histogram.
    
    
     Args:
         hist_angle: a list of point clouds histogram
         PC_list: list all point cloud
         class_label: class label
         number_bin: number of sectors
     returns:
         index_max: index of point cloud mod in list of point clouds
    """  
    N = len(PC_list)
    idc=0
    max_val=0
    indexi=0
    indexj=0
    hist_mod=[0]*N
    sum_hist=0
    index_max=0
    for i in range(N):
        sum_cur=0
        for j in range( N):
            r = np.corrcoef(hist_angle[i], hist_angle[j])
            sum_cur=sum_cur+r[0, 1]
        if(sum_cur>sum_hist):
            sum_hist=sum_cur
            index_max=i
    points = hist_angle[index_max]    
    X_list=[]
    for i in range(number_bin):
        X_list.append(i) 
    if(visualize):
        fig = plt.figure(figsize=(25,5))
        ax1               = fig.add_subplot(1, 4,1)
        ax1.plot(X_list, points/max(points), 'b' ,label='PDF')
        ax1.set_title('Density Histogram Angle, class:{}, #:{}'.format(class_name, index_max))
        ax1.set_xlabel('Number Histogram Angle')
        ax1.set_ylabel('Probablity')
        #ax1.set_ylim((0,0.05))
        ax1.grid() 
        ax1.legend()
        plt.show()
        
    fig = plt.figure(figsize=(25,5))
    points = PC_list[index_max]
    mue = np.mean(points, axis=0)
    
    if(visualize):
        ax1 = fig.add_subplot(1, 4, 1, projection="3d")
        mue = np.mean(points, axis=0)
        ax1.scatter(points[:, 0], points[:, 1], points[:, 2], s=1, c='b', marker="s", facecolor="red", lw=0, alpha=1)
        ax1.scatter(mue[0], mue[1], mue[2], s=100, c='r', marker="*")
        ax1.set_xlabel('X Label')
        ax1.set_ylabel('Y Label')
        ax1.set_zlabel('Z Label')
        ax1.set_xlim((-600,600))
        ax1.set_ylim((-600,600))
        ax1.set_zlim((-400,400))
        print(index_max)
        ax1.set_title('Class:{}, #:{}'.format(class_name,index_max))
        plt.show()
    return index_max,hist_angle[index_max]

def histogram_radius(PC_list, number_rad = 60, distance_bet_radius = 10):
    """
    This function receives a list of point clouds for a class(PEG,ODS,FOCTS,GLAS) and number of radius.
    divide space to several of spherical and calculate number of points in each spherical.
    
     Args:
         PC_list: a list of point clouds
         number_rad: number of radii
     Returns:
         hist_radius: calculated histogram         
    """  
    hist_radius=[]
    for i in range(len(PC_list)):
        cur_hist=[0]*number_rad
        points_normalize=PC_list[i]
        for j in range(len(points_normalize)):
            a=points_normalize[j]
            r=np.sqrt(np.linalg.norm(a, ord=2) ** 2)
            which_r=math.floor(r/distance_bet_radius)
            cur_hist[which_r]=cur_hist[which_r]+1
        cur_hist=np.array(cur_hist)/len(points_normalize)
        hist_radius.append(cur_hist)
    return hist_radius

def draw_histogram_radius(hist_radius, number_rad = 60,sample_per_class = 4, class_name = 'FOCTS', colors = ['g','b','k','y','m','c']):
    """
    This function receives a list of histograms point clouds for a class(PEG,ODS,FOCTS,GLAS) and draws all point clouds.
    In each row draws 4 point clouds. 
    colors contain 6  diffrent colors. 
    
    
     Args:
         hist_radius: a list of point clouds histogram
         class_label: class label
         number_rad: number of sectors
         colors: list of 6 colors for depicting each row    
    """
    N = len(hist_radius)
    idc=0
    scale = 8
    X_list=[]
    for i in range(number_rad):
        X_list.append(i)
        
    i = -1
    fig = plt.figure(figsize=(25,5))
    i = i+1
    sample_per_class = sample_per_class if (len(hist_radius)>sample_per_class) else len(hist_radius)
    for idx in range(sample_per_class):
        points = hist_radius[idx]
        ax1    = fig.add_subplot(1, 4, idx+1)
        ax1.plot(X_list, points/max(points), colors[0] ,label='PDF')
        ax1.set_title('Density Histogram Radius, class:{}, #:{}'.format(class_name, idx+1))
        ax1.set_xlabel('Number Histogram Radius')
        ax1.set_ylabel('Probablity')
        #ax1.set_ylim((0,0.1))
        ax1.legend()
        ax1.grid() 
    plt.show()
    return

def compare_histogram_radius(PC_list, hist_radius, number_rad, class_name, confidence = 0.95, colors = ['g','b','k','y','m','c'], visualize = False):
    
    """
    This function receives a list of histograms point clouds for a class(PEG,ODS,FOCTS,GLAS) and mod point cloud and it's
    histogram.
    
    
     Args:
         hist_radius: a list of point clouds histograms
         PC_list: list all point cloud
         class_label: class label
         number_rad(integer): number of radius
         confidence: level of confidency for ignoring outliers
         visualize(Boolean): visualizing output diagram or not
     returns:
         index_max: index of point cloud mod in list of point clouds
    """  
    N = len(PC_list)
    idc=0
    max_val=0
    indexi=0
    indexj=0
    hist_mod=[0]*N
    sum_hist=0
    index_max=0
    for i in range(N):
        sum_cur=0
        for j in range( N):
            r = np.corrcoef(hist_radius[i],hist_radius[j])
            sum_cur=sum_cur+r[0, 1]
        if(sum_cur>sum_hist):
            sum_hist=sum_cur
            index_max=i
    print("mod is : ",index_max," Value is ",sum_hist/N)
    X_list=[]
    for i in range(number_rad):
        X_list.append(i)
    points = hist_radius[index_max]
    if(visualize):
        fig = plt.figure(figsize=(25,5))
        ax1 = fig.add_subplot(1, 4, 1)
        ax1.plot(X_list, points/max(points), 'b' ,label='PDF')        
        ax1.set_title('Density Histogram Radius, class:{}, #:{}'.format(class_name, index_max))
        ax1.set_xlabel('Number Histogram Radius')
        ax1.set_ylabel('Probablity')
        #ax1.set_ylim((0,0.07))
        ax1.grid() 
        ax1.legend()        
        plt.show()

    points = PC_list[index_max]
    Inlaier_points = []
    mue = my_utils.calc_mue(points)
    sorted_dists = sorted([my_utils.calc_dist(p, mue) for p in points])
    sorted_idx   = sorted([(i, my_utils.calc_dist(p, mue)) for i, p in enumerate(points)], key= lambda d:(d[1], d[0]))

    XY = [(0, 0)]
    t = len(sorted_dists)
    for i, d in enumerate(sorted_dists):
        XY.append(  ( d, 1.*i/t )  )
        if 1.*i/t >= confidence:
            break
    inlaier = [d[0] for d in sorted_idx[:i]]
    inlaier_points = [points[i] for i in inlaier]
    points = np.array(inlaier_points)
    points -= mue
    mue = np.mean(points, axis=0)
    if(visualize):
        fig = plt.figure(figsize=(25,5))
        ax1 = fig.add_subplot(1, 4, 1, projection="3d")
        mue = np.mean(points, axis=0)
        ax1.scatter(points[:, 0], points[:, 1], points[:, 2], s=1, c=colors[0], marker="s", facecolor="red", lw=0, alpha=1)
        ax1.scatter(mue[0], mue[1], mue[2], s=100, c='r', marker="*")
        ax1.set_xlabel('X Label')
        ax1.set_ylabel('Y Label')
        ax1.set_zlabel('Z Label')
        ax1.set_xlim((-600,600))
        ax1.set_ylim((-600,600))
        ax1.set_zlim((-400,400))
        ax1.set_title('Class:{}, #:{}'.format(class_name,index_max))
        plt.show()
    return index_max,hist_radius[index_max]

def map_3D(PC, dim = 2):
    """
    This function receives a list of point clouds and draws point cloud.   
    
     Args:
         PC: point cloud
         dim: map 3D point cloud by removing specific column.
     returns:
         X: mapped point cloud
    """  
    if dim==0:
       temp1      = PC[:, 0]
       PC[:, 0:2] = PC[:, 1:3]
       PC[:, 2]   = temp1
    elif dim==1:
       temp1      = PC[:, 1]
       PC[:, 1] = PC[:, 2]
       PC[:, 2]   = temp1  
    Total    = len(PC)*1.0
    scatter = np.linspace(-300, 300, num=200, endpoint=True, retstep=False, dtype=None, axis=0)
    XY      = []
    X       = []
    total = 0
    for idx, sc in enumerate(scatter):   
        mask = abs(PC[:, 2]-sc)<=5
        XY.append(PC[mask, :2])


    for idx, sc in enumerate(scatter):   
        X.append(len(XY[idx])/Total)
    return X

def hist_angle_map(Sel_Points, number_bin, radius = 200):
    b = np.array([0,0])
    c = np.array([1,0])
    bc = c - b
    hist_angle_map_2d = []
    dim=0
    while(dim<3):
        points   = Sel_Points
        cur_hist = [0]*number_bin
        if(dim == 2):
            x = points[:,0]
            y = points[:,1]
        elif(dim == 1):
            x = points[:,0]
            y = points[:,2]
        elif(dim == 0):
            x = points[:,1]
            y = points[:,2]
        map_points = np.column_stack((x, y))
        count_point_out = 0
        for j in range(len(map_points)):
            a = map_points[j]
            r = np.sqrt(np.linalg.norm(a, ord = 2) ** 2)
            if(r > radius):
                count_point_out = count_point_out+1
                ba = a - b
                cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                angle = np.arccos(cosine_angle)
                if(ba[1] < 0):
                    angle = 360-np.degrees(angle)
                else:
                    angle = np.degrees(angle)
                which_bin = math.floor(angle/(360/number_bin))
                cur_hist[which_bin] = cur_hist[which_bin]+1
        cur_hist = np.array(cur_hist)/count_point_out
        hist_angle_map_2d.append(cur_hist)
        dim = dim+1
    return hist_angle_map_2d


def hist_angle_map_new(Sel_Points, number_bin, radius = 200, n_farest=10):
    b  = np.array([0,0])
    c  = np.array([1,0])
    bc = c - b
    hist_angle_map_2d = []
    dim = 0
    max_dist_2d  = []
    mean_dist_2d = []
    which_point  = []
    while(dim < 3):
        max_point = [[0]*2]*number_bin
        max_dist  = [0]*number_bin
        mean_dist = [0]*number_bin
        cur_hist  = [0]*number_bin
        distance  = [[0]*n_farest]*number_bin
        distance  = np.asarray(distance)
        points    = Sel_Points
        count_point_out = 0
        if(dim == 2):
            x = points[:,0]
            y = points[:,1]
        elif(dim == 1):
            x = points[:,0]
            y = points[:,2]
        elif(dim == 0):
            x = points[:,1]
            y = points[:,2]
        map_points = np.column_stack((x, y))
        for j in range(len(map_points)):
            a = map_points[j]
            r = np.sqrt(np.linalg.norm(a, ord = 2) ** 2)
            if(r >= radius):
                count_point_out = count_point_out + 1
                ba = a - b
                cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                angle = np.arccos(cosine_angle)
                if(ba[1]<0):
                    angle = 360-np.degrees(angle)
                else:
                    angle = np.degrees(angle)
                which_bin = math.floor(angle/(360/number_bin))
                cur_hist[which_bin]  = cur_hist[which_bin] + 1
                mean_dist[which_bin]+= r-radius                
                dist = distance[which_bin]
                dist.sort()
                if(dist[0]<r-radius):
                    dist[0] = r-radius
                    distance[which_bin] = dist
                if( r - radius > max_dist[which_bin]):
                    max_dist[which_bin]  = r-radius
                    max_point[which_bin] = a
        for index in range(number_bin):
            if(cur_hist[index] > 0):
                if (np.count_nonzero(distance[index]) > 0):
                    mean_dist[index] = sum(distance[index])/np.count_nonzero(distance[index])
                else:
                    mean_dist[index] = 0
        max_dist_2d.append(max_dist)
        mean_dist_2d.append(mean_dist)
        which_point.append(max_point)
        cur_hist = np.array(cur_hist)/count_point_out
        hist_angle_map_2d.append(cur_hist)
        dim = dim + 1
    return hist_angle_map_2d, max_dist_2d, mean_dist_2d, which_point


def convex_hull_2D(Sel_Points, surface, ratio_list, hist_map_2d_distance, Max_Point_Sector, Centralized, number_bin = 60, colors = ['g', 'b', 'k','y','m','c'], radius = 200, mean_dist = False):
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N           = 4
    C           = len(surface)
    dim         = 0
    while(dim < 3):
        fig   = plt.figure(figsize=(25,5))
        for idx in range(len(surface)):
            maximum_point = Max_Point_Sector[dim*3+idx]
            ratio         = ratio_list[dim*3 + idx]
            points = Sel_Points[idx]
            ax1    = fig.add_subplot(1, len(surface), idx+1)
            if(dim == 2 ):
                x = points[:,0]
                ax1.set_xlabel('X Label')
                y = points[:,1]
                ax1.set_ylabel('Y Label')
            elif(dim == 1):
                x = points[:,0]
                ax1.set_xlabel('X Label')
                y = points[:,2]
                ax1.set_ylabel('Z Label')
            elif(dim == 0):
                x = points[:,1]
                ax1.set_ylabel('Y Label')
                y = points[:,2]            
                ax1.set_ylabel('Z Label')
            ax1          = fig.add_subplot(1, len(surface), idx+1)            
            spacing      = 100
            minorLocator = MultipleLocator(spacing)
            ax1.yaxis.set_minor_locator(minorLocator)
            ax1.xaxis.set_minor_locator(minorLocator)
            ax1.grid(which = 'minor')            
            ax1.scatter(x, y, s = 2, c = colors[idx], marker = "s", facecolor = "red", lw = 0, alpha = 1.)
            if(Centralized):
                ax1.set_xlim((-800,800))
                ax1.set_ylim((-800,800))
                circle2 = plt.Circle((0, 0), radius, color = 'r', fill = False, linewidth = 3)
                ax1.add_patch(circle2)
                v_origin = [radius,0]
                v        = [350,0]
                x_pose   = []
                y_pose   = []
                for angle in range(number_bin):
                    k    = ratio[angle]
                    v[0] = (maximum_point[angle])[0]
                    v[1] = (maximum_point[angle])[1]                    
                    if( np.sqrt(v[0]**2+v[1]**2) >= radius):
                        xp       = (radius*v[0])/(np.sqrt(v[0]**2+v[1]**2))
                        yp       = (radius*v[1])/(np.sqrt(v[0]**2+v[1]**2))
                        v_origin = [xp, yp]
                    else:
                        theta    = (angle+0.5)*360/number_bin
                        xp       = radius*math.cos(math.radians(theta))
                        yp       = radius*math.sin(math.radians(theta))
                        v_origin = [xp,yp]
                        v        = v_origin
                    if mean_dist:
                        v[0] = xp + k*(v[0] - xp)                  
                        v[1] = yp + k*(v[1] - yp)
                    x_pose.append(v[0])
                    y_pose.append(v[1])
                    ax1.arrow(v_origin[0], v_origin[1] ,v[0]-v_origin[0], v[1]-v_origin[1],
                              head_width = 10, head_length = 10,
                              lw = 2,  fc = '#777777',  ec = '#777777')#,  length_includes_head = True)
                ax1.plot(x_pose, y_pose)
            ax1.grid(which = 'minor', alpha = 0.2)
            ax1.grid(which = 'major', alpha = 0.8)
            if mean_dist:
                ax1.set_title('convex hull mean distance Class:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))
            else:
                ax1.set_title('convex hull max distanceClass:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))                
        dim=dim+1
    plt.show()
    return

def get_convex_hull_pose(Sel_Points, surface, ratio_list, Max_Point_Sector, Centralized, number_bin = 60, colors = ['g', 'b', 'k','y','m','c'], radius = 200, mean_dist = False):
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N           = 4
    C           = len(surface)
    dim         = 0
    all_pose_convex_hull = []
    while(dim<3):
        for idx in range(len(surface)):
            maximum_point = Max_Point_Sector[dim*3+idx]
            ratio         = ratio_list[dim*3 + idx]
            points = Sel_Points[idx]
            if(dim==2):
                x = points[:,0]
                y = points[:,1]
            elif(dim==1):
                x = points[:,0]
                y = points[:,2]
            elif(dim==0):
                x = points[:,1]
                y = points[:,2]
            spacing      = 100
            minorLocator = MultipleLocator(spacing)
            if(Centralized):
                v_origin = [radius,0]
                v        = [350,0]
                x_pose   = []
                y_pose   = []
                for angle in range(number_bin):
                    k    = ratio[angle]
                    v[0] = (maximum_point[angle])[0]
                    v[1] = (maximum_point[angle])[1]
                    if(np.sqrt(v[0]**2+v[1]**2)>=radius):
                        xp       = (radius*v[0])/(np.sqrt(v[0]**2+v[1]**2))
                        yp       = (radius*v[1])/(np.sqrt(v[0]**2+v[1]**2))
                        v_origin = [xp,yp]
                    else:
                        theta    = (angle+0.5)*360/number_bin
                        xp       = radius*math.cos(math.radians(theta))
                        yp       = radius*math.sin(math.radians(theta))
                        v_origin = [xp,yp]
                        v        = v_origin
                    if mean_dist:
                        v[0] = xp + k*(v[0] - xp)                  
                        v[1] = yp + k*(v[1] - yp)
                    x_pose.append(v[0])
                    y_pose.append(v[1])
                x_pose = np.array(x_pose)
                y_pose = np.array(y_pose)
                temp = np.vstack((x_pose,y_pose))
                all_pose_convex_hull.append(temp)
        dim=dim+1
    return all_pose_convex_hull

def histogram_convex_hull_2D(Sel_Points, surface, ratio_list, Max_Point_Sector, Centralized, number_bin = 60, colors = ['g', 'b', 'k','y','m','c'], radius = 200, mean_dist = False):
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N           = 4
    C           = len(surface)
    dim         = 0
    while(dim<3):
        fig = plt.figure(figsize=(25,5))
        for idx in range(len(surface)):
            maximum_point = Max_Point_Sector[dim*3+idx]
            ratio         = ratio_list[dim*3 + idx]
            hist_convex_hull=[]
            points = Sel_Points[idx]
            ax1    = fig.add_subplot(1, len(surface), idx+1)
            if(dim==2):
                x = points[:,0]
                y = points[:,1]
            elif(dim==1):
                x = points[:,0]
                y = points[:,2]
            elif(dim==0):
                x = points[:,1]
                y = points[:,2]            
            ax1          = fig.add_subplot(1, len(surface), idx+1)            
            spacing      = 100
            minorLocator = MultipleLocator(spacing)
            if(Centralized):
                v_origin = [radius,0]
                v        = [350,0]
                x_pose   = []
                y_pose   = []
                maximum_point = Max_Point_Sector[dim*3+idx]
                for angle in range(number_bin):
                    k    = ratio[angle]
                    v[0] = (maximum_point[angle])[0]
                    v[1] = (maximum_point[angle])[1]
                    if(np.sqrt(v[0]**2+v[1]**2)>=radius):
                        xp       = (radius*v[0])/(np.sqrt(v[0]**2+v[1]**2))
                        yp       = (radius*v[1])/(np.sqrt(v[0]**2+v[1]**2))
                        v_origin = [xp,yp]
                    else:
                        theta    = (angle+0.5)*360/number_bin
                        xp       = radius*math.cos(math.radians(theta))
                        yp       = radius*math.sin(math.radians(theta))
                        v_origin = [xp,yp]
                        v        = v_origin
                    if mean_dist:
                        v[0] = xp + k*(v[0] - xp)                  
                        v[1] = yp + k*(v[1] - yp)
                    x_pose.append(v[0])
                    y_pose.append(v[1])                    
                    hist_convex_hull.append(np.sqrt( (v[0]-v_origin[0])**2 + (v[1]-v_origin[1])**2 ) )
                center = np.linspace(0, number_bin-1, num = number_bin)
                ax1.bar(center, hist_convex_hull/np.max(hist_convex_hull))        
            ax1.grid(which = 'minor', alpha = 0.2)
            ax1.grid(which = 'major', alpha = 0.8)
            ax1.set_ylabel('Frequeny')
            ax1.set_xlabel('bins')
            if mean_dist:
                ax1.set_title('convex hull mean distance Class:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))
            else:
                ax1.set_title('convex hull max distance Class:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))
        dim=dim+1
    plt.show()
    return

def get_histogram_convex_hull_2D(Sel_Points, surface, Max_Point_Sector, Centralized, number_bin = 60, colors = ['g', 'b', 'k','y','m','c'], radius = 200):
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N           = 4
    C           = len(surface)
    dim         = 0
    HIST = []
    while(dim<3):
        fig = plt.figure(figsize=(25,5))
        for idx in range(len(surface)):
            hist_convex_hull=[]
            points = Sel_Points[idx]
            if(dim==2):
                x = points[:,0]
                y = points[:,1]
            elif(dim==1):
                x = points[:,0]
                y = points[:,2]
            elif(dim==0):
                x = points[:,1]
                y = points[:,2]            
            ax1          = fig.add_subplot(1, len(surface), idx+1)            
            spacing      = 100
            minorLocator = MultipleLocator(spacing)
            if(Centralized):
                v_origin = [radius,0]
                v        = [350,0]
                x_pose   = []
                y_pose   = []
                maximum_point = Max_Point_Sector[dim*3+idx]
                for angle in range(number_bin):
                    v[0]          = (maximum_point[angle])[0]
                    v[1]          = (maximum_point[angle])[1]
                    if(np.sqrt(v[0]**2+v[1]**2)>=radius):
                        xp       = (radius*v[0])/(np.sqrt(v[0]**2+v[1]**2))
                        yp       = (radius*v[1])/(np.sqrt(v[0]**2+v[1]**2))
                        v_origin = [xp,yp]
                    else:
                        theta    = (angle+0.5)*360/number_bin
                        xp       = radius*math.cos(math.radians(theta))
                        yp       = radius*math.sin(math.radians(theta))
                        v_origin = [xp,yp]
                        v        = v_origin
                    x_pose.append(v[0])
                    y_pose.append(v[1])
                    hist_convex_hull.append(np.sqrt( (v[0]-v_origin[0])**2 + (v[1]-v_origin[1])**2 ) )
                center = np.linspace(0, number_bin-1, num = number_bin)
                HIST.append(hist_convex_hull/np.max(hist_convex_hull))
        dim=dim+1
    return HIST

def display_convex_image(all_pose_convex_hull, surface, planes):
    k=0
    while (k < len(all_pose_convex_hull)):
        fig = plt.figure(figsize=(25,5))
        for l in range(len(surface)):
            ax1    = fig.add_subplot(1, len(surface), l+1)
            minimum=np.amin(all_pose_convex_hull[k],axis=1)
            new_calue = all_pose_convex_hull[k]-2*(np.min(minimum))
            new_calue = np.asarray(new_calue)
            new_calue = np.floor(new_calue)
            width, height = 2000 , 2000
            image = np.zeros((height, width))
            line_thickness = 2
            X = new_calue[0,:]
            Y = new_calue[1,:]
            X1 = X[0]
            Y1 = Y[0]
            i  = 1
            while(i < len(X) ):
                X2 = X[i]
                Y2 = Y[i]
                image = cv2.line(image, (int(X1), int(Y1)), (int(X2), int(Y2)), (255, 255, 255), thickness=line_thickness)
                X1 = X2
                Y1 = Y2
                i  = i+1
            X2 = X[0]
            Y2 = Y[0]
            image = cv2.line(image, (int(X1), int(Y1)), (int(X2), int(Y2)), (255, 255, 255), thickness=line_thickness)
            image = np.flipud(image)
            kernel = np.ones((5,5), np.uint8)
            img_dilation = cv2.dilate(image, kernel, iterations=1)
            img_dilation = img_dilation.astype(np.uint8)
            contour,hier = cv2.findContours(img_dilation, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
            for cnt in contour:
                cv2.drawContours(img_dilation,[cnt],0,255,-1)
            gray = cv2.bitwise_not(img_dilation)
            plt.imshow(img_dilation)
            k = k+1
            if(k>=len(all_pose_convex_hull)):
                break
        if(k>=len(all_pose_convex_hull)):
            break
    return