import numpy as np
import math

def makeTranslation(x, y, z):
    """
    return the 4x4 numpy array transpose of the 
    3x4 
    """
    m = np.eye(4)
    m[:3,3] = x,y,z
    return m

def makeRotationX(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    m = np.eye(4)
    m[1,1] = c
    m[1,2] = -s
    m[2,1] = s
    m[2,2] = c
    return m

def makeRotationY(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    m = np.eye(4)
    m[0,0] = c
    m[0,2] = s
    m[2,0] = -s
    m[2,2] = c
    return m

def makeRotationZ(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    m = np.eye(4)
    m[0,0] = c
    m[0,1] = -s
    m[1,0] = s
    m[1,1] = c
    return m

def makeQuaternionMatrix(x, y, z, w):
    n = math.sqrt(x**2 + y**2 + z**2 + w**2)
    qx = x/n
    qy = y/n
    qz = z/n
    qw = w/n
    m = np.eye(3)
    m[0, 0] =  1 - 2*qy**2 - 2*qz**2
    m[0, 1] = 2*qx*qy - 2*qz*qw
    m[0, 2] = 2*qx*qz + 2*qy*qw
    m[1, 0] = 2*qx*qy + 2*qz*qw
    m[1, 1] = 1 - 2*qx**2 - 2*qz**2
    m[1, 2] = 2*qy*qz - 2*qx*qw
    m[2, 0] = 2*qx*qz - 2*qy*qw
    m[2, 1] = 2*qy*qz + 2*qx*qw
    m[2, 2] = 1 - 2*qx**2 - 2*qy**2
    return m
# end def

def applyQuaternion(coords, m3):
    """ assume we are applying rows of vectors
    [               
    [x1, y1, z1],
        ...             * m3.T 
    [xn, yn, zn]
    ]
    """
    return np.dot(coords, m3.T)
# end def

def applyGeomStack(coords, stack):
    pass
# end def

def applyTransform(coords, m4):
    """ assume we are applying rows of vectors
    [               
    [x1, y1, z1, 1.],
        ...             * m4[:3, :].T 
    [xn, yn, zn, 1.]
    ]

    where m4 is the total 4x4 transfomation matrix
    ['n11','n12', 'n13', 'n14',
     'n21', 'n22', 'n23', 'n24',
     'n31', 'n32', 'n33', 'n34',
     'n41', 'n42', 'n43', 'n44'])
    """

    rows, cols = coords.shape
    # print(rows, cols, coords.dtype)
    stacked_coords = np.hstack((coords, np.ones((rows, 1))))
    m4x3 = m4[:3,:].T   # cut off last row so we end up with a (rows , 3) product 
    return np.dot(stacked_coords, m4x3)


