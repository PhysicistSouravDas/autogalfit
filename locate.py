from astropy import units as u
from astropy.coordinates import SkyCoord

def locate(ra, dec, return_obj=False):
    """
    ra: in hours (eg 01 25 28.27)
    dec: in deg (eg -00 06 23.03)
    return_obj: Return SkyCoord object if True, else return ra, dec value in tuple
    frame is 'icrs'
    returns the SkyCoord object with coordinates in degrees
    eg: <SkyCoord (ICRS): (ra, dec) in deg (345., -90.)>
    """

    c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
    if return_obj:
        return c
    return c.ra.deg, c.dec.deg


if __name__ == "__main__":
    
    g5_coord = locate(ra='01 25 28.27', dec='-00 06 23.03')

    ds9_head = '# Region file format: DS9 version 4.1\n\
    global color=black dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n\
    icrs\n'
    ra_deg = g5_coord[0]
    dec_deg = g5_coord[1]

    ds9_file = open('ds9.reg', mode='w')
    ds9_file.write(ds9_head)
    ds9_file.write(f'circle({ra_deg}, {dec_deg} 1") # text={{g5}}\n')
    ds9_file.close()
