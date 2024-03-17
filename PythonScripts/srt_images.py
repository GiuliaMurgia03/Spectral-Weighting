import matplotlib.pyplot as plt
import numpy
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
from astropy.io import fits
from matplotlib import ticker
from matplotlib.colors import PowerNorm
from matplotlib.colors import LogNorm

plt.style.use(astropy_mpl_style)

path_srt="/Users/giuliamurgia/Desktop/RPOL_29JUN16/"
files_srt=["splat","weighted_splat","splat_flagged_merge", "weighted_splat_4L","weighted_splat_4R","weighted_splat_LL+RR", "weighted_splat_LL_all","weighted_splat_RR_all","weighted_splat_LL+RR_all"]

for image_name in files_srt:

    #Open file
    image_file = path_srt+image_name+".fits"
    image_data = fits.getdata(image_file, ext=0)
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)
    # Drop the last two axes
    wcs = wcs.dropaxis(3)
    wcs = wcs.dropaxis(2)
    plane_image_data = image_data[0,0, :, :]
    print(plane_image_data.shape)


    #Plot 
    plt.figure()
    plt.subplot(projection=wcs) # Real RA/DEC axes
    plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.5, vmin=-0.001, vmax=0.06, clip=True))
    plt.grid(False) #Remove Grid
    font = {'family':'sans','color':'black','size':12}

    #Colorbar
    cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.set_label('Intensity', size=10) #Title of the Colorbar
    cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
    #Axis
    plt.xlabel("RA [J2000]", fontdict = font)
    plt.ylabel("DEC [J2000]", fontdict = font)
    plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

    plt.savefig("/Users/giuliamurgia/Desktop/NEW_NoiseRelative_Plots/"+"srt_"+image_name+".pdf", format="pdf", bbox_inches="tight")


quit()


#Open file Channel 2342 RA
image_file = path_srt+"29JUN16_M31_RA1_RR.FITS"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,2342, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.5, vmin=-0.1, vmax=1.4, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_RA_2342.pdf", format="pdf", bbox_inches="tight")



#Open file Channel 2342 DEC
image_file = path_srt+"29JUN16_M31_DEC3_RR.FITS"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,2342, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.5, vmin=-0.1, vmax=1.4, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_DEC_2342.pdf", format="pdf", bbox_inches="tight")



#Open file Channel 2342 Merge
image_file = path_srt+"weighted_merge.fits"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,2342, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.5, vmin=-0.1, vmax=1.4, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_merge_2342.pdf", format="pdf", bbox_inches="tight")


#Open file Channel 8890 Merge
image_file = path_srt+"weighted_merge.fits"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,8899, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.5, vmin=-0.1, vmax=1.4, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_merge_8890.pdf", format="pdf", bbox_inches="tight")

#Open file Channel 8890 Flagged Merge
image_file = path_srt+"flagged_merge.fits"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,8890, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='afmhot', norm=PowerNorm(0.5, vmin=-0.1, vmax=1.4, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Intensity', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_merge_flagged_8890.pdf", format="pdf", bbox_inches="tight")


#Open file Weights Channel 2342 RA
image_file = path_srt+"weights_29JUN16_M31_RA1_RR.FITS"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,2342, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='viridis', norm=LogNorm(vmin=10, vmax=10000, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Weight', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_weights_RA_2342.pdf", format="pdf", bbox_inches="tight")



#Open file Weights Channel 2342 DEC
image_file = path_srt+"weights_29JUN16_M31_DEC3_RR.FITS"
image_data = fits.getdata(image_file, ext=0)
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
# Drop the last two axes
wcs = wcs.dropaxis(3)
wcs = wcs.dropaxis(2)
plane_image_data = image_data[0,2342, :, :]
print(plane_image_data.shape)

#Plot 
plt.figure()
plt.subplot(projection=wcs) # Real RA/DEC axes
plt.imshow(plane_image_data, cmap='viridis', norm=LogNorm( vmin=0.01, vmax=100, clip=True))
plt.grid(False) #Remove Grid
font = {'family':'sans','color':'black','size':12}

#Colorbar
cbar=plt.colorbar(shrink=0.58, pad=0.07, location='top')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label('Relative Weight', size=10) #Title of the Colorbar
cbar.ax.tick_params(direction='out', length=2, width=0.5, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=9)  #Ticks and Ticks Parameters of the Colorbar
#Axis
plt.xlabel("RA [J2000]", fontdict = font)
plt.ylabel("DEC [J2000]", fontdict = font)
plt.tick_params(direction='out', length=3, width=1, colors='black',
                grid_color='black', grid_alpha=0.5,labelsize=10)

plt.savefig("/Users/giuliamurgia/Desktop/NEW_Plots/"+"srt_weights_DEC_2342.pdf", format="pdf", bbox_inches="tight")

