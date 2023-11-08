import matplotlib.pyplot as plt
import numpy
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
from matplotlib import ticker
from matplotlib.colors import PowerNorm
from matplotlib.colors import LogNorm
import math

def createList(r1, r2):
    return [item for item in range(r1, r2+1)]
channels=createList(0,1023)

channel_noise=0.00853032
splat_noise=channel_noise/math.sqrt(0.00853032)

path_sim="/Users/giuliamurgia/Desktop/simul_RA_DEC/"
files_sim=["merge", "weighted_merge"]




 # Open file
image_file = path_sim+"smooth_sources.fits"
image_data = fits.getdata(image_file, ext=0)
image_data=image_data/splat_noise
spectrum_sky=[]

for c in channels:
    sum=0    
    for x in range(255):
        for y in range(255):
            val=image_data[c,y,x]
            sum=sum+val
    spectrum_sky.append(sum/(256*256))


for image_name in files_sim:

    # Open file
    image_file = path_sim+image_name+".fits"
    image_data = fits.getdata(image_file, ext=0)
    image_data=image_data/splat_noise

    spectrum=[]

    for c in channels:
        sum=0    
        for x in range(255):
            for y in range(255):
                val=image_data[c,y,x]
                sum=sum+val
        spectrum.append(sum/(256*256))

    # Plot
    font = {'family':'sans','color':'black','size':8}
    fig, ax = plt.subplots()
    if(image_name=='merge'):
        plt.plot(channels, spectrum,color='darkmagenta', label='Merge', linewidth=0.7)
    if(image_name=='weighted_merge'):
        plt.plot(channels, spectrum,color='darkmagenta', label='Weighted Merge', linewidth=0.7)    
    plt.plot(channels, spectrum_sky,color='limegreen', label='Sky Model', linewidth=0.7)
    leg = plt.legend(fontsize='8')
    plt.tick_params(labelsize=8)
    ax.set_xlim(0,1023)
    ax.set_ylim(0.000,0.015)
    ax.set_box_aspect(1/2.5)

    plt.xlabel("Channels", fontdict = font)
    plt.ylabel("Relative Intensity", fontdict = font)
    plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/spectrum_"+image_name+".pdf", format="pdf", bbox_inches="tight")
    plt.clf()


    
