
import math as math
from numpy.lib.scimath import sqrt as csqrt

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import scipy.signal as sg
from scipy.special import sph_harm


class Target(object):
    """ This is a public, abstract class. All instantiations of targets will inherit and or modify
        a subset of the methods inside of the Target class
    """

    def __init__(self, name, desc, n_pix, anti_aliasing ):
        """ Target constructor
            The target is constructed using the characteristic dimension....

        :param name: Name of target
        :param desc: Description of target
        :param n_pix: The number pixels at which to render target
        :param anti_aliasing: oversampling to reduce aliasing

        :return:
        """
        self.name = name
        self.desc = desc
        self.n_pix = n_pix
        self.anti_aliasing = anti_aliasing

    def get_name(self):
        return self.name

    def get_desc(self):
        return self.desc
    
    def get_n_pix(self):
        return self.n_pix
    
    def get_AntiAliasing(self):
        return self.anti_aliasing

class GreyscaleLambertianBeachballTarget(Target):
    """ This class inherits the generic(abstract) Target class and adds to
        it specialised methods specific to this type of target. An instantiation
        of this class will render a greyscale beachball. This target is a spherical,lambertian reflector
        that is neutral (grey) across all wavelengths.
    """
    def __init__(self, diam=1.0, n_sector_pairs=12, min_refl=0.0, max_refl=1.0, n_pix=128, anti_aliasing=1, zero_padding):
        """ Create a greyscale lambertian beachball target.

        :param diam: Diameter of the target in m
        :param n_sector_pairs: Number of dark/light sector pairs
        :param min_refl: Reflectance of the dark sectors
        :param max_refl: Reflectance of the light sectors
        :param n_pix: Number of pixels (horizontal and vertical) to render target
        :param anti_aliasing: is set to 1 by default, this means there is no oversampling. Can be set to an integer
            value to oversample the target and then downsample to n_pix in order to reduce aliasing effects (jagged
            edges) along the sector boundaries.
        :return:
        """
        Target.__init__(self, 'BeachBall', 'Grey Spherical Lambertian', n_pix= n_pix, anti_aliasing= anti_aliasing)
        
        self.diam = diam
        self.n_sector_pairs = n_sector_pairs
        self.min_refl = min_refl
        self.max_refl = max_refl

    def get_diam(self):
        return self.diam

    def get_n_sector_pairs(self):
        return self.n_sector_pairs

    def __str__(self):
        printString = """ The target is a %s %s with %i sectors, rendered with %i pixels and oversampled by a factor %i."""         % (self.desc, self.name, self.n_sector_pairs, self.n_pix, self.anti_aliasing)
        return str(printString)
    
    def get_target(self):
        """ Build a reflectance image of the target and calculate surface normals at every pixel

        """
        # Set up pixel relative x-coodinates
        x = np.linspace(-1.0, 1.0, self.n_pix * self.anti_aliasing)
        # Meshgrid the x and y coordinates as the same
        x, y = np.meshgrid(x, x)
        # Calculate the polar coordinate theta of each pixel
        theta = np.arctan2(y, x)
        # Produce a siemens star target with the given number of spokes/sectors
        siemans_star = np.cos(self.Nsector_pairs * theta)
        # Create an array that will become a target mask
        mask = np.ones(siemans_star.shape)
        # Make the mask binary by taking the sign of the siemens star
        siemans_star =  np.sign(siemans_star)/2.0 + 0.5
        # Set all value outside the unit circle to the new value
        circle = np.sqrt(x**2 + y**2) > 1.0
        mask[circle] = 0.0
        # downsample if aliasing requested
        if self.anti_aliasing > 1.0:
            the_filter = np.ones((self.anti_aliasing, self.anti_aliasing)) / self.anti_aliasing**2
            # Apply 2D filter convolution
            siemans_star = sg.convolve2d(siemans_star, the_filter, mode = 'full')
            # Convole the mask as well
            mask = sg.convolve2d(mask, the_filter, mode = 'full')
            siemans_star = siemans_star[0.0::self.anti_aliasing, ::self.anti_aliasing]
            mask = mask[0.0::self.anti_aliasing, ::self.anti_aliasing]
        x = np.linspace(-1.0, 1.0, siemans_star.shape[0])
        y = np.linspace(-1.0, 1.0, siemans_star.shape[1]) 
        x, y = np.meshgrid(x, y)            
        # Create blank list of normal vector components
        normal_vector_list = [[],[],[]]
        normal_vector_list[0] = x
        normal_vector_list[1] = y
        normal_vector_list[2] = (csqrt(1.0 - x**2 - y**2)).real
        # Save reflectance map, normal vectors and mask in self
        self.normal_vector_list = normal_vector_list
        self.siemens_star = siemans_star
        self.mask = mask
        return normal_vector_list, mask, siemans_star

test = GreyscaleLambertianBeachballTarget(Npix=128*4, anti_aliasing=1)

vector, grey_scale_mask, star = test.get_target()

xn = vector[0] * grey_scale_mask
yn = vector[1] * grey_scale_mask
zn = vector[2] * grey_scale_mask

def render_target(self, vza, svaz, normal_vector, rad_env, rad_fields):
    """
    This function will render a target in a Radiant Environment Map (REM)
    
    """
    
    if vza >= 0.0 and vza <= 90.0:
        
        targ_REM = rad_env.compilation[1].rad_fields[0]
        
        polar_ang = np.linspace(0.0, np.pi, targ_REM.shape[0])  
        azi_ang = np.linspace(0.0, np.pi, targ_REM.shape[1])
        sph_fit = sph_harm(2, 2, azi_ang, polar_ang)
        
        targ_SHREM = rad_env.compilation[1].rad_fields[0] *sph_fit
        targ_edir = rad_env.compilation[1].edir
        
    elif vza >90.0 and vza <=180.0:
        
        targ_REM = rad_env.compilation[0].rad_fields[0]
        
        polar_ang = np.linspace(0.0, np.pi, targ_REM.shape[0])  
        azi_ang = np.linspace(0.0, np.pi, targ_REM.shape[1])
        sph_fit = sph_harm(2, 2, azi_ang, polar_ang)
        
        targ_SHREM = rad_env.compilation[0].rad_fields[0] *sph_fit
        targ_edir = rad_env.compilation[0].edir
    else:
    
        pass # log an error / warning here
        

    target_mask = [[]]
    at_target_radiance = [[]]
    
#   Compute the sun-pointing vector 

    sza = rad_env.sza
    sun_vectors = [[-1.0 * np.sin(np.deg2rad(sza))], [0.0], [np.cos(np.deg2rad(sza))]]
        REM = np.ones((5,5))
#   Compute the polar and azimuth angles for the REM
#     polar_ang = np.linspace(0.0, np.pi, REM.shape[0])  #
#     azi_ang = np.linspace(0.0, np.pi, REM.shape[1])
    
#   Do vector geometry
              
#   Rotate beachball surface normals into the VZA
#   Rotating into the VZA is a rotation about the x-axis
    
    x_rot = np.array(np.pi-vza)
    
#     Then rotate surface normals into the solar-relative view
#     azimuth angle. This is a rotation about the z-axis
    
    z_rot = svaz
              
    new_normal_vectors = ([[normal_vector[0]], [normal_vector[1]], [normal_vector[2]]] *  x_rot) * z_rot 
   
    cos_to_sun_vectors = new_normal_vectors * sun_vectors
    cos_to_sun_vectors[cos_to_sun_vectors < 0.0] = 0.0
    
    # PSEUDOCODE

    # if target is spectrally neutral then replicate the reflectance matrix in the spectral dimension use np.tile
    #
    # extract band number from band defintion from REM this is then saved in REMBands
    # extract band number from band model ,this is then saved in Targbands
    #
    # check bands function -> if REM bands in Targ bands return true else retutrn false, and issue warning
    # 
new_render = render_target(vector)