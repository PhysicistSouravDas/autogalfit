# program to get the pixel values from the ra dec coordinates of an object
from astropy.io import fits, ascii
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from astropy.visualization import (imshow_norm, MinMaxInterval, SqrtStretch, PercentileInterval)
from astropy.coordinates import SkyCoord
from astropy import units as u

from locate import locate

def get_pixels(img_name, coordinate, hdu=0):
    """
    Obtains the pixel value (xpix, ypix) using coordinates
    img_name (str): name of the image file
    coordinate (float): (ra, dec) value in degrees
    hdu (int): the data unit containing the image (default 0 for cutouts)
    Returns tuple of (x, y) pixel values
    """
    hdul = fits.open(img_name)
    img_header = hdul[hdu].header
    img_data = hdul[hdu].data
    ra, dec = coordinate

    c = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

    try:
        wcs = WCS(img_header)
        pix_x, pix_y = wcs.world_to_pixel(c)
    except ValueError:
        wcs = WCS(img_header).celestial
        pix_x, pix_y = wcs.world_to_pixel(c)
    
    return pix_x, pix_y




if __name__ == '__main__':

    hdul = fits.open('g5_cutout.fits')
    img_headers = hdul[0].header
    img = hdul[0].data

    header_file = open('g5_header.txt', mode='w')
    header_file.write(repr(img_headers))
    header_file.close()

    # getting ra and dec of the g5 galaxy
    g5_ra_hourang = '01 25 28.27'  # given in research paper pg6 table3 doi:10.1088/0004-637X/811/2/132
    g5_dec_deg = '-00 06 23.03'
    g5_coord = locate(ra=g5_ra_hourang, dec=g5_dec_deg, return_obj=True)

    ### plotting (object oriented)
    # plt.imshow(img, origin='lower', )
    wcs = WCS(img_headers).celestial  # .celestial is for WCS more than 3 axes
    pix_x, pix_y = wcs.world_to_pixel(g5_coord)
    print(pix_x, pix_y)
    # print(wcs[:2])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax = plt.subplot(projection=wcs)
    im, norm = imshow_norm(img, ax, origin='lower', interval=PercentileInterval(99), cmap='rainbow',)# stretch=SqrtStretch())
    ax.grid(color='white', ls='solid')
    ax.set_xlabel('Galactic Longitude')
    ax.set_ylabel('Galactic Latitude')
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='white', ls='dotted')
    overlay[0].set_axislabel('Right Ascension (J2000)')
    overlay[1].set_axislabel('Declination (J2000)')
    plt.show()
