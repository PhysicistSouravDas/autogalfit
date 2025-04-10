"""
.. module:: cutout

:Synopsis: Create cutout FITS images from larger fits image and plot them.
:Author: Sourav Das
:Year: 2023
"""

from pathlib import Path
from astropy.io import fits
from astropy import units as u

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D

from matplotlib import pyplot as plt
import numpy as np

DIR = Path('./').absolute()  # current directory
im_path = f'{DIR}/HST/hst_13024_19_acs_wfc_f814w/hst_13024_19_acs_wfc_f814w_drz.fits.gz'

def cut_galaxy(fits_file: str, coordinate: tuple, size: int, output_fits: str) -> None:
    """
    Function which creates cutout fits image of galaxy from a large fits image
    fits_file: Input whole fits file name with path
    coordinate: Coordinate of the galaxy, ra, dec in degrees in tuple (ra, dec) in floats
    size: Size of the cutout square image
    output_fits: name of the output cutout file to be saved
    Returns `None`
    """
    hdul = fits.open(fits_file)
    img = hdul[1].data
    ra, dec = coordinate
    position = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    wcs = WCS(hdul[1].header)
    cutout = Cutout2D(img, position, (size, size), wcs=wcs)
    # print(type(cutout.data))
    # https://stackoverflow.com/a/33515855 --> cutout.wcs.to_header()
    out_hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())#, header=hdul[1].header)
    out_hdu.writeto(output_fits, overwrite=True)

if __name__ == '__main__':
    # name of qso source --> HE0153-4520
    hdul = fits.open(im_path)
    print(hdul.info())
    img = hdul[1].data
    # print(img)
    # print(np.min(img), np.max(img))
    # plt.imshow(
    #     img, vmin=-0.03, vmax=0.1, origin='lower',
    #     cmap='viridis'
    #     )
    # plt.show()

    # We have the coordinates of few galaxies stored in separate coord_ra.txt and 
    # coord_dec.txt

    ra_file = open('coord_ra.txt', 'r')
    ra = []

    for line in ra_file.readlines():
        ra_val = line.lstrip().rstrip()
        ra.append(float(ra_val))

    ra_file.close()

    dec_file = open('coord_dec.txt', 'r')
    dec = []

    for line in dec_file.readlines():
        dec_val = line.lstrip().rstrip()
        dec.append(float(dec_val))

    dec_file.close()

    # initiate preparing cutout one by one using the coordinates
    position = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    wcs = WCS(hdul[1].header)
    plt.subplot(projection=wcs)
    plt.imshow(
        img, vmin=-0.03, vmax=0.1, origin='lower',
        cmap='viridis'
        )
    plt.grid(color='white', ls='solid')
    plt.show()

    cutout = Cutout2D(img, position[6], (200, 200), wcs=wcs)
    plt.imshow(
        cutout.data, 
        vmin=-0.03, vmax=0.1, origin='lower',
        cmap='viridis'
        )
    plt.show()
