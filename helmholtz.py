import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal 

def masked_axicon_zplot(zlim=[-0.2,1], npoints=60, maskrad=1):
        ## Make the forward beam 
        light_forward = Light(type='axicon',N=512)
        ## Move to focus
        light_forward.propagate_after_lens(z=0)
        ## Apply the mask and propagate
        light_forward.apply_mask(maskrad=maskrad)
        light_forward.show_intensity()
        zlimf = zlim[1]
        npointsf = np.floor(npoints*zlimf/np.sum(np.abs(zlim)))
        zprofile_forward = light_forward.zplot_after_mask(zlim=zlimf, npoints = npointsf)

        ## Make the reverse beam 
        light_reverse= Light(type='axiconrev',N=512)
        ## Move to focus
        light_reverse.propagate_after_lens(z=0)
        ## Apply the mask and propagate
        light_reverse.apply_mask(maskrad=maskrad)
        zlimr = np.abs(zlim[0])
        npointsr = np.abs(np.floor(npoints*zlimr/np.sum(np.abs(zlim))))
        zprofile_reverse = light_reverse.zplot_after_mask(zlim=zlimr, npoints = npointsr)

        ## Put everything together
        full_light = np.hstack(([np.fliplr(zprofile_forward.transpose()),zprofile_reverse.transpose()]))

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        ax.imshow(full_light, aspect=0.2)

        plt.xlabel('z')
        plt.ylabel('r')
        plt.title('Propagation from axicon with a mask')

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)

        plt.show()

        return


class Light():
    
    def __init__(self, N=1024, type='axicon'):
        ## Initial properties
        self.wavelength = 0.532 # all units in um
        self.resolution = N
        
        ## Optical parameters:
        axicon_number = 10
        vortex_number = 6
        spiralaxicon_number = 10
        
        # Frequency space
        freqx = np.linspace(-10,10,N)-1/N
        self.freqx, self.freqy = np.meshgrid(freqx, freqx)
        
        # Real space
        x = np.linspace(-10,10,N)
        self.X, self.Y = np.meshgrid(x,x)
        self.R = np.sqrt(np.power(self.X,2)+np.power(self.Y,2))
        self.theta = np.arctan(self.X/self.Y)
        
        # Initial field
        if type=='axicon':
            self.initial_field = np.exp(-0.5*np.power(self.R,2))*\
                np.exp(-axicon_number*1j*self.R)
                
        if type=='axiconrev':
            self.initial_field = np.exp(-0.5*np.power(self.R,2))*\
                np.exp(axicon_number*1j*self.R)
                       
        elif type=='spiralaxicon':
            self.initial_field = np.exp(-0.5*np.power(self.R,2))*\
                np.exp(spiralaxicon_number*1j*self.R + vortex_number*1j*self.theta)
                
        elif type=='pseudo-LG':
            self.initial_field = np.exp(-0.5*np.power(self.R,2))*\
                np.exp(vortex_number*1j*self.theta)
                
        self.Uz = self.initial_field
        self.Iz = np.abs(self.initial_field)
        
    def transfer_function(self, z):
        H = np.exp(1j*2*math.pi*(z/self.wavelength) * \
            np.sqrt(1-np.power((self.wavelength*self.freqx),2) - \
            np.power((self.wavelength*self.freqy),2)+0j))
        return H

    def propagate(self, z=10):
        
        ## Calculate the transfer function
        H = self.transfer_function(z)
        ## Propagate
        A0 = np.fft.fftshift(np.fft.fft2(self.initial_field))
        Az = A0*H
        self.Uz = np.fft.ifft2(np.fft.ifftshift(Az))
        self.Iz = np.abs(self.Uz)
        
    def propagate_after_mask(self, z=10):
        
        ## Calculate the transfer function
        H = self.transfer_function(z)
        ## Propagate
        A0 = np.fft.ifft2(np.fft.ifftshift(self.initial_field))
        Az = A0*H
        self.Uz = np.fft.fftshift(np.fft.fft2(Az))
        self.Iz = np.abs(self.Uz)
        
    def propagate_after_lens(self, z=10):
        
        ## Distance from axicon to lens
        z01 = 0.1;
        ## Focal length of lens
        f = 3;
        ## Distance from lens to z=0 image plane origin
        z23 = f+z
        
        ## Calculate the transfer functions
        H01 = self.transfer_function(z01)
        H23 = self.transfer_function(z)
        
        ## Go to fourier space
        A0 = np.fft.fftshift(np.fft.fft2(self.initial_field))
    
        ## Propagate to lens
        A1 = A0*H01   
        
        ## Apply Fourier lens
        A2 = np.fft.ifft2(np.fft.ifftshift(A1))
        
        ## Propagate away from lens
        A3 = A2*H23
            
        ## Go back to real space 
        self.Uz = np.fft.fftshift(np.fft.fft2(A3))
        self.Iz = np.abs(self.Uz)
    
    def propagate_after_lens_2(self, z=10):
        
        ## Focal length of lens
        f = 0.5;
        
        ## Distance from axicon to lens
        z01 = f;
        ## Distance from lens to z=0 image plane origin
        z45 = z
        
        ## Calculate the transfer functions
        H01 = self.transfer_function(z01)
        H45 = self.transfer_function(z45)
        
        ## Go to fourier space
        A0 = np.fft.fftshift(np.fft.fft2(self.initial_field))
    
        ## Propagate to lens
        A1 = A0*H01   
        
        ## Apply lens
        A2 = np.fft.ifft2(np.fft.ifftshift(A1))
        A3 = A2*np.exp(-1j*(2*math.pi/self.wavelength)*np.power(self.R,2)/(2*f))
        A4 = np.fft.fftshift(np.fft.fft2(A3))

        ## Propagate from lens to output
        A5 = A4*H45
            
        ## Go back to real space 
        self.Uz = np.fft.ifft2(np.fft.ifftshift(A3))
        self.Uz = A4
        self.Iz = np.abs(self.Uz)
        
    def show_intensity(self):
        fig, ax = plt.subplots(figsize=(5,5))
        cax = ax.imshow(self.Iz,cmap='gray')
        plt.title('Intensity')
        plt.show()
    
    def show_phase(self):
        fig, ax = plt.subplots(figsize=(5,5))
        phase = np.angle(self.Uz)
        plt.title('Phase')
        ax.imshow(phase, cmap='hsv')
        plt.show()
        
    def show_zplot(self, zlim=20, npoints=10):
        zprofile = []
        for z in np.linspace(0,zlim,npoints):
            self.propagate_after_lens(z)
            cut = self.Iz[math.floor(self.resolution/2),:]
            zprofile.append(cut)
            
        return np.vstack(zprofile)
    
    def show_zplot_2(self, zlim=20, npoints=10):
        zprofile = []
        for z in np.linspace(0,zlim,npoints):
            self.propagate_after_lens_2(z)
            cut = self.Iz[math.floor(self.resolution/2),:]
            zprofile.append(cut)
            
        return np.vstack(zprofile)
    
    def show_zplot_nolens(self, zlim=20, npoints=10):
        zprofile = []
        for z in np.linspace(0,zlim,npoints):
            self.propagate(z)
            cut = self.Iz[math.floor(self.resolution/2),:]
            zprofile.append(cut)
            
        return np.vstack(zprofile)

    def zplot_after_mask(self, zlim=20, npoints=10):
        zprofile = []
        for z in np.linspace(0,zlim,npoints):
            self.propagate_after_mask(z)
            cut = self.Iz[math.floor(self.resolution/2),:]
            zprofile.append(cut)
            
        return np.vstack(zprofile)

    def apply_mask(self, maskrad=1):
        ## Prepare the mask
        mask = maskrad-self.R
        mask[mask>0]=0
        mask[mask<0]=1
        self.initial_field = self.Uz*mask
        self.Uz = self.initial_field
        self.Iz = np.abs(self.initial_field)



