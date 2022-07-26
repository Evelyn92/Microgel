import numpy as np
import math
import cv2
import pandas as pd
import mahotas
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage.measure import regionprops
from skimage.transform import rescale, resize
from skimage.morphology import thin
from skimage.morphology import skeletonize
from skimage.morphology import medial_axis
from matplotlib.ticker import MultipleLocator

def zernike_moment(img, radius = 10, degree = 8):
    """
    degree = 8 generates 25 samples moment while degree = 81 generates
    """    
#     regions = regionprops(img)
#     for props in regions:
#         print(f'props.centroid is {props.centroid} and props.bbox is {props.bbox}')        
#         print(props.moments_hu)
#         print(props.moments_central)
#         print(props.moments_normalized)
    value  = mahotas.features.zernike_moments(img, radius, degree = degree)
#     mahotas.features.zernike_moments(im, radius, degree=8, cm={center_of_mass(im)})
    return value
	
def centralize_convex(img, resizeing = True):
    """
    degree = 8 generates 25 moment while degree = 16 generates  81 moments
    
    """
    radius = 10
    regions = regionprops(img)
    for props in regions:
        center       = props.centroid
        bbox         = props.bbox
        a            = [center[0]-bbox[0], center[1]-bbox[1]]
        lower_radius = np.sqrt(np.linalg.norm(a, ord=2) ** 2)
        b            = [bbox[2]-center[0], bbox[3]-center[1]]
        upper_radius = np.sqrt(np.linalg.norm(b, ord=2) ** 2)
        radius       = int(max(lower_radius,upper_radius))
    crop_image       = img[int(center[0])- radius:int(center[0])+ radius ,
                           int(center[1])- radius:int(center[1])+ radius]
    if resizeing :
        crop_image   = resize(crop_image, (1300, 1300), anti_aliasing=True)
    crop_image_size  = crop_image.shape
    new_image        = np.zeros(img.shape)
    new_image_size   = new_image.shape
    new_image [new_image_size[0]//2-crop_image_size[0]//2:new_image_size[0]//2+crop_image_size[0]//2,
               new_image_size[1]//2-crop_image_size[1]//2:new_image_size[1]//2+crop_image_size[1]//2] = crop_image
    return new_image, radius

def calculate_moment(im):
    # Calculate Moments
    moments = cv2.moments(im)
    # Calculate Hu Moments
    huMoments = cv2.HuMoments(moments)
    # Log scale hu moments
    for i in range(0,7):
        huMoments[i] = -1* math.copysign(1.0, huMoments[i]) * math.log10(abs(huMoments[i]))
    return moments, huMoments
	
def display_convex_image(all_pose_convex_hull, surface, planes, radius = 10, degree = 8, resizeing = True):
    k=0
    Hu_moment_dataframe      = []
    zernike_moment_dataframe = []
    list_image_convex_hull   = []
    list_contours            = []
    while (k < len(all_pose_convex_hull)):
        fig = plt.figure(figsize=(25,5))
        hu_moment = []
        zernike_mom = []
        for l in range(len(surface)):
            minimum   = np.amin(all_pose_convex_hull[k],axis=1)
            new_calue = all_pose_convex_hull[k]-2*(np.min(minimum))
            new_calue = np.asarray(new_calue)
            new_calue = np.floor(new_calue)
            width, height  = 2000 , 2000
            image          = np.zeros((height, width))
            line_thickness = 2
            X  = new_calue[0,:]
            Y  = new_calue[1,:]
            X1 = X[0]
            Y1 = Y[0]
            i  = 1
            while(i < len(X) ):
                X2    = X[i]
                Y2    = Y[i]
                image = cv2.line(image, (int(X1), int(Y1)), (int(X2), int(Y2)), (255, 255, 255), thickness=line_thickness)
                X1    = X2
                Y1    = Y2
                i     = i+1
            X2 = X[0]
            Y2 = Y[0]
            image  = cv2.line(image, (int(X1), int(Y1)), (int(X2), int(Y2)), (255, 255, 255), thickness=line_thickness)
            image  = np.flipud(image)
            kernel = np.ones((5,5), np.uint8)
            img_dilation  = cv2.dilate(image, kernel, iterations=1)
            img_dilation  = img_dilation.astype(np.uint8)
            contour, hier = cv2.findContours(img_dilation, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
            list_contours.append(contour)
            for cnt in contour:
                cv2.drawContours(img_dilation,[cnt],0,255,-1)
            gray = cv2.bitwise_not(img_dilation)
            new_image, radius = centralize_convex(img_dilation, resizeing)
            
            zernike_mom.append(zernike_moment(new_image, radius, degree))
            moments, huMoments = calculate_moment(new_image)
#             img_dilation = np.concatenate( img_dilation, axis=0 )
#             huMoments = huMoments.tolist()
#             print(f'moments are {moments}')
#             print(f'huMoments are {huMoments}, shape is {type(huMoments)}')
            hu_moment.append(huMoments)
#             df = pd.DataFrame(huMoments[0], index =['a', 'b', 'c', 'd', 'e', 'f', 'g'], columns =[str(surface[l])])
            
            list_image_convex_hull.append(new_image)
        
            ax1       = fig.add_subplot(1, len(surface), l+1)
            ax1.set_title('Convex Hull Class:{}, Plane:{}'.format(surface[l], planes[k//len(surface)]))
            plt.imshow(new_image)
            new_image = 255 * new_image
            new_image = new_image.astype(np.uint8)
#             cv2.imwrite(surface[l]+planes[k//len(surface)]+'.jpg',  new_image)
            k = k+1
            if(k>=len(all_pose_convex_hull)):
                break
        df = pd.DataFrame(list(zip(hu_moment[0], hu_moment[1], hu_moment[2])), 
                          columns =[surface[0], surface[1], surface[2]])
        Hu_moment_dataframe.append(df)
        
        df = pd.DataFrame(list(zip(zernike_mom[0], zernike_mom[1], zernike_mom[2])), 
                          columns =[surface[0], surface[1], surface[2]])
        zernike_moment_dataframe.append(df)
        if(k>=len(all_pose_convex_hull)):
            break
    return Hu_moment_dataframe, zernike_moment_dataframe, list_image_convex_hull, list_contours
	
def convert_arrayofarray_to_array(arr):
    output = []
    for item in arr:
        output.append(item[0])
    return output

def histogram_moment(moment_dataframes, surface, colors = ['g', 'b', 'k','y','m','c'], hu_moment = True ):
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    N           = 4
    C           = len(surface)
    dim         = 0
    bar_width   = 0.5
    while(dim<3):
        dataframe = moment_dataframes[dim]
        fig = plt.figure(figsize=(25,5))
        for idx in range(len(surface)):
            hist_moment  = dataframe[surface[idx]]
            ax1          = fig.add_subplot(1, len(surface), idx+1)
            if hu_moment:
                hist_moment  = convert_arrayofarray_to_array(hist_moment)
                ax1.set_title('Hu Moment Class:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))
            else :
                ax1.set_title('Zernike Moment Class:{}, Plane:{}'.format(surface[idx], CLASS_Plane[dim]))
            spacing      = 100
            minorLocator = MultipleLocator(spacing)
            center = np.linspace(0, len(hist_moment), num = len(hist_moment))
            index = np.arange(len(hist_moment))
            ax1.bar(index, hist_moment/np.max(hist_moment), bar_width, fc='b') #width=0.5)
            ax1.grid(which = 'minor', alpha = 0.2)
            ax1.grid(which = 'major', alpha = 0.8)
            ax1.set_ylabel('Value')
            ax1.set_xlabel('moment')                
        dim=dim+1
    plt.show()
    return
	
def contour_distance(list_image_convex_hull, surface):
    CLASS_Plane = ['YZ Plane','XZ Plane','XY Plane']
    dim = 0    
    while (dim < 3):
        FOCTS    = list_image_convex_hull[dim*3 + 0]
        PEG      = list_image_convex_hull[dim*3 + 1]
        ODS      = list_image_convex_hull[dim*3 + 2]
        dis_FP   = cv2.matchShapes(FOCTS, PEG, cv2.CONTOURS_MATCH_I3, 0.0)
        dis_FO   = cv2.matchShapes(FOCTS, ODS, cv2.CONTOURS_MATCH_I3, 0.0)
        dis_PO   = cv2.matchShapes(PEG,   ODS, cv2.CONTOURS_MATCH_I3, 0.0)
        print(f'plane {CLASS_Plane[dim]} distance FOCTS and PEG is {dis_FP}')
        print(f'plane {CLASS_Plane[dim]} distance FOCTS and ODS is {dis_FO}')        
        print(f'plane {CLASS_Plane[dim]} distance PEG and ODS is {dis_PO}')
        
#         d2 = cv2.matchShapes(im1,im2,cv2.CONTOURS_MATCH_I2,0)
#         d3 = cv2.matchShapes(im1,im2,cv2.CONTOURS_MATCH_I3,0)
        dim += 1
    return
	
