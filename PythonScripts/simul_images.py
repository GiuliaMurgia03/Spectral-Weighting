import matplotlib.pyplot as plt
import numpy
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
from matplotlib import ticker
from matplotlib.colors import PowerNorm
from matplotlib.colors import LogNorm
import math


plt.style.use(astropy_mpl_style)

path_sim="/Users/giuliamurgia/Desktop/simul_RA_DEC/"
files_sim=["splat_merge_norfi","splat_smooth_sources_noise_rfi1", "splat_smooth_sources_noise_rfi2", "splat_merge", "splat_weighted_merge"]

channel_noise=0.00853032
splat_noise=channel_noise/math.sqrt(1024)

for image_name in files_sim:

    #Open file
    image_file = path_sim+image_name+".fits"
    image_data = fits.getdata(image_file, ext=0)
    image_data=image_data/splat_noise
    plane_image_data = image_data[0, :, :] # Remove third axis
    print(image_name)


    #Plot 
    plt.figure()
    plt.imshow(plane_image_data, cmap='afmhot',norm=PowerNorm(0.4, vmin=-0.15, vmax=15, clip=True))
    plt.grid(False) #Remove Grid
    font = {'family':'sans','color':'black','size':12}

    #Colorbar
    cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.set_label('Relative Intensity', size=10) #Title of the Colorbar
    cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
    #Axis
    plt.xlabel("RA [pixel]", fontdict = font)
    plt.ylabel("DEC [pixel]", fontdict = font)
    plt.gca().invert_yaxis() #Inverte y axis
    plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)
    tick_locator = ticker.MaxNLocator(nbins=6) #Number of ticks
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+image_name+".pdf", format="pdf", bbox_inches="tight")



#Open file for Example Cube 1 Channel 510
image_file = path_sim+"smooth_sources_noise_rfi1"+".fits"
image_data = fits.getdata(image_file, ext=0)
image_data=image_data/channel_noise
plane_image_data = image_data[510,:,:] # Take channel 510
print("smooth_sources_noise_rfi1")

#Plot 
plt.figure()
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.7, vmin=-0.25, vmax=3, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [pixel]", fontdict = font)
plt.ylabel("DEC [pixel]", fontdict = font)
plt.gca().invert_yaxis() #Inverte y axis
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+"Cube1Channel510.pdf", format="pdf", bbox_inches="tight")




#Open file for Example Cube 2 Channel 950
image_file = path_sim+"smooth_sources_noise_rfi2"+".fits"
image_data = fits.getdata(image_file, ext=0)
image_data=image_data/channel_noise
plane_image_data = image_data[950,:,:] # Take channel 950
print("smooth_sources_noise_rfi2")


#Plot 
plt.figure()
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.7, vmin=-0.25, vmax=3, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [pixel]", fontdict = font)
plt.ylabel("DEC [pixel]", fontdict = font)
plt.gca().invert_yaxis() #Inverte y axis
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+"Cube2Channel950.pdf", format="pdf", bbox_inches="tight")



#Open file for Example Weights Cube 1 Channel 510
image_file = path_sim+"weights_smooth_sources_noise_rfi1"+".fits"
image_data = fits.getdata(image_file, ext=0)
plane_image_data = image_data[510,:,:] # Remove third axis
biggest = numpy.amax(plane_image_data)
plane_image_data=plane_image_data/biggest
print("weights_smooth_sources_noise_rfi1")


#Plot 
plt.figure()
plt.imshow(plane_image_data, cmap='viridis', norm=LogNorm(vmin=0.001, vmax=0.25, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Weight', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [pixel]", fontdict = font)
plt.ylabel("DEC [pixel]", fontdict = font)
plt.gca().invert_yaxis() #Inverte y axis
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+"WeightsCube1Channel510.pdf", format="pdf", bbox_inches="tight")




#Open file for Example Weights Cube 2 Channel 950
image_file = path_sim+"weights_smooth_sources_noise_rfi2"+".fits"
image_data = fits.getdata(image_file, ext=0)
plane_image_data = image_data[950,:,:] # Remove third axis
biggest = numpy.amax(plane_image_data)
plane_image_data=plane_image_data/biggest

print(plane_image_data.shape, "weights_smooth_sources_noise_rfi2")


#Plot 
plt.figure()
plt.imshow(plane_image_data, cmap='viridis', norm=LogNorm(vmin=0.001, vmax=0.25, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Weight', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [pixel]", fontdict = font)
plt.ylabel("DEC [pixel]", fontdict = font)
plt.gca().invert_yaxis() #Inverte y axis
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+"WeightsCube2Channel950.pdf", format="pdf", bbox_inches="tight")


# Residuals Image
image_file = path_sim+"residuals.fits"
image_data = fits.getdata(image_file, ext=0)
image_data=image_data/splat_noise
plane_image_data = image_data[0, :, :] # Remove third axis
print(plane_image_data.shape, "residuals.fits")


#Plot 
plt.figure()
#plt.imshow(plane_image_data, cmap='afmhot',norm=PowerNorm(0.4, vmin=-0.0005, vmax=0.02, clip=True))
#plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(1, vmin=-0.0008, vmax=0.0008, clip=True))
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(1.15, vmin=None, vmax=None, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Intensity', size=10) #Title of the Colorbar
#cbar.set_ticks([-0.00008,0,0.00008])
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
tick_locator = ticker.MaxNLocator(nbins=3) #Number of ticks
cbar.locator = tick_locator
cbar.update_ticks()
#Axis
plt.xlabel("RA [pixel]", fontdict = font)
plt.ylabel("DEC [pixel]", fontdict = font)
plt.gca().invert_yaxis() #Inverte y axis
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+"residuals.pdf", format="pdf", bbox_inches="tight")

