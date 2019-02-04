from numpy import mat
from math import cos, sin, radians

def rot_a(angle):
    """Input: Rotation angle in degrees
    Returns: Rotation matrix required to rotate with the given angle around the a-axis
    """
    angle = radians(angle)
    return mat([[1, 0         ,           0],
                [0, cos(angle), -sin(angle)],
                [0, sin(angle),  cos(angle)]])


def rot_b(angle):
    """Input: Rotation angle in degrees
    Returns: Rotation matrix required to rotate with the given angle around the b-axis
    """
    angle = radians(angle)
    return mat([[cos(angle) ,0 ,sin(angle)],
                [0          ,1 ,         0],
                [-sin(angle),0 ,cos(angle)]])

def rot_c(angle):
    """Input: Rotation angle in degrees
    Returns: Rotation matrix required to rotate with the given angle around the b-axis
    """
    angle = radians(angle)
    return mat([[cos(angle), -sin(angle), 0],
                [sin(angle),  cos(angle), 0],
                [0         , 0          , 1]])
                  
def inversion():

    return mat([[-1,  0,  0],
                [ 0, -1,  0],
                [ 0,  0, -1]])
                  
def mirror_b():

    return mat([[1 ,  0,  0],
                [0 , -1,  0],
                [0 ,  0,  1]])