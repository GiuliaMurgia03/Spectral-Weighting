import matplotlib.pyplot as plt
import numpy
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
from matplotlib import ticker
from matplotlib.colors import PowerNorm
from matplotlib.colors import LogNorm
from astropy.wcs import WCS
import math

def createList(r1, r2):
    return [item for item in range(r1, r2+1)]
channels=createList(2000,3000)

path_sim="/Users/giuliamurgia/Desktop/RPOL_29JUN16/"
files_sim=["29JUN16_M31_RA1_RR", "29JUN16_M31_DEC3_RR","weighted_merge"]


for image_name in files_sim:

    # Open file
    image_file = path_sim+image_name+".fits"
    image_data = fits.getdata(image_file, ext=0)
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)
    # Drop the axis
    wcs = wcs.dropaxis(3)
    cube_image_data = image_data[0,:, :, :]
    spectrum=[]

    for c in channels:
        sum=0    
        for x in range(16):
            for y in range(16):
                val=cube_image_data[c,97+(y-8),101+(x-8)]
                sum=sum+val
        spectrum.append(sum/(16*16))

    # Plot
    font = {'family':'sans','color':'black','size':8}
    fig, ax = plt.subplots()
    if (image_name=="29JUN16_M31_RA1_RR"):
        plt.plot(channels, spectrum, color='darkmagenta', label='RA',linewidth=0.6)    
    if (image_name=="29JUN16_M31_DEC3_RR"):
        plt.plot(channels, spectrum, color='darkmagenta', label='DEC',linewidth=0.6)  
    if (image_name=="weighted_merge"):
        plt.plot(channels, spectrum, color='darkmagenta', label='Weighted Merge',linewidth=0.6)         
    plt.tick_params(labelsize=8)
    ax.set_xlim(2000,2800)
    ax.set_ylim(-0.1,1.1)
    ax.set_box_aspect(1/2.5)
    leg = plt.legend(fontsize='8')

    plt.xlabel("Channels", fontdict = font)
    plt.ylabel("JY/BEAM", fontdict = font)
    plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/srt_spectrum_"+image_name+".pdf", format="pdf", bbox_inches="tight")
    plt.clf()


    
