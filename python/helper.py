import numpy as np
from cv2 import blur, distanceTransform, erode
from cv2 import DIST_L2, DIST_MASK_PRECISE
from skimage.feature import peak_local_max
from scipy.stats import scoreatpercentile
from collections import Counter
import cvxpy as cp
import matplotlib.pyplot as plt

def yy(imgSize, centered = 0):
    ny, nx = imgSize
    if centered:
        y = np.arange(-(ny//2),(ny+1)//2).reshape(-1,1)
    else:
        y = np.arange(ny).reshape(-1,1)
    return np.tile(y, [1, nx])  


def xx(imgSize, centered = 0):
    ny, nx = imgSize
    if centered:
        x = np.arange(-(nx//2),(nx+1)//2)
    else:
        x = np.arange(nx)
    return np.tile(x, [ny, 1])  


def elipsBound(cy, cx, a, b, d, imgSize):
    ny, nx = imgSize
    dy = np.abs((b**2*np.cos(d) + a**2*np.tan(d)*np.sin(d)) /
                np.sqrt(b**2 + a**2*np.tan(d)**2))
    dx = np.abs((a**2*np.cos(d) + b**2*np.tan(d)*np.sin(d)) / 
                np.sqrt(a**2 + b**2*np.tan(d)**2))
    dy[np.isnan(dy)] = a[np.isnan(dy)]
    dx[np.isnan(dx)] = b[np.isnan(dx)]
    y_min = np.maximum(np.floor(cy-dy),0).astype(np.int64)
    x_min = np.maximum(np.floor(cx-dx),0).astype(np.int64)
    y_max = np.minimum(np.ceil(cy+dy),ny-1).astype(np.int64)
    x_max = np.minimum(np.ceil(cx+dx),nx-1).astype(np.int64)
    return y_min, x_min, y_max, x_max


def concaveDetect(poly):
    diff1 = poly - np.roll(poly, 1, axis = 0)
    diff2 = np.roll(diff1, -1, axis = 0)
    idx = diff1[:,:,0] * diff2[:,:,1] - diff1[:,:,1] * diff2[:,:,0] > 0
    cocv = poly[np.squeeze(idx),0,:]
    return {(x[0], x[1]) for x in cocv}


def findEllipse(img, rto = np.arange(1, 4, 0.2), ang = np.arange(0, 180, 5),
                minArea = None):
    if minArea == None:
        minArea = img.sum() / 100
    
    npair = len(ang) * len(rto)
    angle = np.tile(ang, [1, len(rto)]).flatten()
    angle = angle / 180 * np.pi
    ratio = np.tile(rto.reshape(-1,1), [1, len(ang)]).flatten()

    cos, sin = np.cos(angle), np.sin(angle)
    cos2, sin2 = cos ** 2, sin ** 2
    cossin = cos * sin
    Ta, Tinva = (1 / ratio - 1), ratio - 1
    T, Tinv = np.zeros([npair,2,2]), np.zeros([npair,2,2])
    T[:,0,0],T[:,1,1],T[:,0,1] = Ta*sin2+1,Ta*cos2+1,Ta*cossin
    T[:,1,0] = T[:,0,1]
    Tinv[:,0,0],Tinv[:,1,1],Tinv[:,0,1] = Tinva*sin2+1, Tinva*cos2+1, Tinva*cossin
    Tinv[:,1,0] = Tinv[:,0,1]
    
    ycord, xcord = yy(img.shape, centered = 1), xx(img.shape, centered = 1)
    ref = np.reshape([ycord[0,0], xcord[0,0]], [1,2,1])
    imgBool = img.astype(bool)
    Tcord = np.round(T.dot(np.stack((ycord[imgBool], xcord[imgBool])))) - ref
    upLeft = np.squeeze(np.min(Tcord, axis=2));
    lowRight = np.squeeze(np.max(Tcord, axis=2));
    Tcord = (Tcord - upLeft[:,:,np.newaxis]).astype(np.int64)
    imageSize = (lowRight - upLeft + 1).astype(np.int64)
    parID = np.arange(npair)[imageSize.min(axis=1) > 1]
    
    flagRatio = 0
    elist = []
    for i in parID:
        if ratio[i] == 1:
            if flagRatio:
                continue
            else:
                flagRatio = 1
        edt = np.zeros(imageSize[i], dtype = np.uint8)
        edt.ravel()[np.ravel_multi_index(Tcord[i], imageSize[i])] = 1
        edt = distanceTransform(edt, DIST_L2, 3)
        peak = peak_local_max(blur(edt, (3,3)))
        cord = peak + upLeft[i,:] + ref.squeeze()
        cord = Tinv[i,:,:].dot(cord.T) - ref.squeeze(axis=0)                    # round ?
        minor = edt.ravel()[np.ravel_multi_index(peak.T, imageSize[i])]
        major = minor * ratio[i]
        idxArea = (np.pi * major * minor) > minArea
        elips = np.zeros([idxArea.sum(), 5])
        elips[:,:2] = cord.T[idxArea,:]
        elips[:,2] = major[idxArea]
        elips[:,3] = minor[idxArea]
        elips[:,4] = angle[i]
        elist.append(elips)
    return np.concatenate(elist)


def overlay(cntrImg, elist, prcntl = [30,50,70], coverRate = 0.95):
    ny, nx = cntrImg.shape
    cy, cx, a, b, d = elist.T
    y_min, x_min, y_max, x_max = elipsBound(cy, cx, a, b, d, cntrImg.shape)
    edt = distanceTransform(1-cntrImg, DIST_L2, DIST_MASK_PRECISE)
    ycord, xcord = yy(cntrImg.shape), xx(cntrImg.shape)
    scoreAll = np.full([len(prcntl), ny, nx], np.inf)
    areaAll = np.zeros([len(prcntl), ny, nx])
    overlap = -1 * np.ones([len(prcntl), ny, nx], dtype = np.int32)
    idx3D = ny * nx * np.arange(len(prcntl)).reshape(-1,1)
    for i in range(len(elist)):
        edti = edt[y_min[i]:(y_max[i]+1), x_min[i]:(x_max[i]+1)]
        y = ycord[y_min[i]:(y_max[i]+1), x_min[i]:(x_max[i]+1)]
        x = xcord[y_min[i]:(y_max[i]+1), x_min[i]:(x_max[i]+1)] 
        cy, cx, a, b, d = elist[i]
        elipsIn = (((x-cx)*np.cos(d)+(y-cy)*np.sin(d))/a)**2 + \
            (((y-cy)*np.cos(d)-(x-cx)*np.sin(d))/b)**2 < 1
        elipsOut = (((x-cx)*np.cos(d)+(y-cy)*np.sin(d))/(a-1))**2 + \
            (((y-cy)*np.cos(d)-(x-cx)*np.sin(d))/(b-1))**2 > 1
        score = scoreatpercentile(edti[np.logical_and(elipsIn, elipsOut)],
                                  prcntl)
        
        cordElips = np.array(elipsIn.nonzero()) + [[y_min[i]],[x_min[i]]]
        idx = np.ravel_multi_index(cordElips, cntrImg.shape) + idx3D
        score = np.tile(score, [idx.shape[1],1]).T
        idxCover = np.logical_and(scoreAll.ravel()[idx] == score,
                                  areaAll.ravel()[idx] < a*b)
        idxCover = np.logical_or(scoreAll.ravel()[idx] > score, idxCover)
        idx = idx[idxCover]
        scoreAll.ravel()[idx] = score[idxCover];
        areaAll.ravel()[idx] = a*b;
        overlap.ravel()[idx] = i;
    remainID = set()
    minArea = (1 - coverRate) * np.pi * elist[:,2] * elist[:,3]
    for i in range(len(prcntl)):
        remainID.update({j for j, A in Counter(overlap[i,:,:].ravel()).items()
                         if A > minArea[j]})
    remainID.remove(-1)
    return elist[list(remainID),:]


def computeDist(elist, cntrs, cocvs, imgSize):
    distMat = []
    y, x = yy(imgSize), xx(imgSize)
    for i, cntr in enumerate(cntrs):
        dist = np.zeros([len(cntr), len(elist)])
        cntr = [cntr[:,0,1], cntr[:,0,0]]
        for j in range(len(elist)):
            cy, cx, a, b, d = elist[j]
            elips = (((x-cx)*np.cos(d)+(y-cy)*np.sin(d))/a)**2 + \
                (((y-cy)*np.cos(d)-(x-cx)*np.sin(d))/b)**2 < 1;
            elipsCntr = 1 - elips * (1-erode(elips.astype(np.uint8), None))
            elipsEdt = distanceTransform(elipsCntr, DIST_L2, DIST_MASK_PRECISE)
            dist[:,j] = elipsEdt.ravel()[np.ravel_multi_index(cntr, imgSize)]
        if len(cocvs[i]) < 2:
            distMat.append(dist.sum(axis=0, keepdims=True))
        else:
            distSeg = np.zeros([len(cocvs[i]), len(elist)])
            for j in range(len(cocvs[i])-1):
                distSeg[j,:] = dist[cocvs[i][j]:cocvs[i][j+1],:].sum(axis=0)
            distSeg[-1,:] = dist[cocvs[i][-1]:,:].sum(axis=0) + \
                dist[:cocvs[i][0],:].sum(axis=0)
            distMat.append(distSeg)
    return np.concatenate(distMat)


def intergerProgamming(distMat, lam):
    M, N = distMat.shape
    x = cp.Variable(shape=(M,N), boolean=True)
    z = cp.Variable(shape=(1,N), boolean=True)
    obj = cp.sum(cp.multiply(distMat, x)) + lam * cp.sum(z)
    cnstrs = []
    cnstrs.append(cp.sum(x, axis=0, keepdims=True) >= z)
    cnstrs.append(cp.sum(x, axis=0, keepdims=True) <= z*M)
    cnstrs.append(cp.sum(x, axis=1) == 1)
    opt = cp.Problem(cp.Minimize(obj), cnstrs)
    _ = opt.solve('GLPK_MI')
    return np.nonzero(z.value.ravel())[0]


def plotEllipse(img, elist):
    img = (img != 0).astype(np.uint8)
    ny, nx = img.shape
    cy, cx, a, b, d = elist.T
    y_min, x_min, y_max, x_max = elipsBound(cy, cx, a, b, d, img.shape)
    ycord, xcord = yy(img.shape), xx(img.shape)
    for i in range(len(elist)):
        y = ycord[y_min[i]:(y_max[i]+1), x_min[i]:(x_max[i]+1)]
        x = xcord[y_min[i]:(y_max[i]+1), x_min[i]:(x_max[i]+1)] 
        cy, cx, a, b, d = elist[i]
        elips = (((x-cx)*np.cos(d)+(y-cy)*np.sin(d))/a)**2 + \
            (((y-cy)*np.cos(d)-(x-cx)*np.sin(d))/b)**2 < 1
        elipsCntr = np.logical_and(elips, 1-erode(elips.astype(np.uint8), None))
        cordElips = np.array(elipsCntr.nonzero()) + [[y_min[i]],[x_min[i]]]
        img.ravel()[np.ravel_multi_index(cordElips, img.shape)] = 2
    return img
    