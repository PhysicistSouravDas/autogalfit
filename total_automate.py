# this program will create the cutouts of all the galaxies
# obtain the pixel position of the galaxies in that cutout
# and prepare the GALFIT input files to fit the galaxies in each cutout

# the coordinates and other data is stored in coordinates.txt file

from astropy.io import ascii
import numpy as np
from matplotlib import pyplot as plt
import subprocess, argparse, time

# local imports
from cutout import cut_galaxy
from get_pixels import get_pixels
from autogalfit import run_galfit, fetch_params
from plot_dist import plot_dist


### PART I : Creating the cutouts

init_time = time.time()
# parsing argument to catch the obs_id from the terminal
parser = argparse.ArgumentParser()
parser.add_argument(
    'fname', type=str, help='Enter the file name (obs_id) with extension'
    )
parser.add_argument('psf_file', type=str, help='Enter psf file name without extension')
args = parser.parse_args()

# obs_id = 'hst_14660_02_acs_wfc_f814w_jd8e02_drc.fits.gz'
obs_id = str(args.fname)
table = ascii.read('coordinates.txt', format='rst')

c = float(table[0]['g_ra']), float(table[0]['g_dec'])
# """
count = 0
print('Preparing the cutouts ...')

# checking whether a folder already exists or not
# for succesful command returncode is 0
if subprocess.run(['ls', "results"]).returncode != 0:
    subprocess.run(['mkdir', "results"])
if subprocess.run(['ls', "cutouts"]).returncode != 0:
    subprocess.run(['mkdir', "cutouts"])
# galfit_output is already created by run_galfit(); but still no harm in below code
if subprocess.run(['ls', "galfit_output"]).returncode != 0:
    subprocess.run(['mkdir', "galfit_output"])
if subprocess.run(['ls', "results/plots"]).returncode != 0:
    subprocess.run(['mkdir', "results/plots"])
if subprocess.run(['ls', "results/data"]).returncode != 0:
    subprocess.run(['mkdir', "results/data"])
if subprocess.run(['ls', "results/raw_data"]).returncode != 0:
    subprocess.run(['mkdir', "results/raw_data"])

# storing all the cutouts in cutouts/ folder
for row in table:
    coord = float(row['g_ra']), float(row['g_dec'])
    gal_name = str(row['gal'])
    cut_galaxy(obs_id, coord, 200, f'cutouts/{gal_name}.fits')
    count += 1

print(f'{count} galaxies have been cut out.')
# """
### PART II : Getting pixel values of galaxies
# getting pixel values for feeding in GALFIT
print('Extracting pixel values ...')

pix_x = []
pix_y = []

for row in table:
    coord = float(row['g_ra']), float(row['g_dec'])
    xpix, ypix = get_pixels(f"cutouts/{row['gal']}.fits", coord)
    pix_x.append(float(xpix))
    pix_y.append(float(ypix))

print('Extraction done')

### PART III : Auto GALFIT Execution

# program to run galfit by changing slightly few parameters
# we want to change the following parameters
# TODO: Using autogalfit's function to fit a galaxy, and obtain its parameters.
# """
count = 0
psf_file = str(args.psf_file)
for row in table:
    # running GALFIT n no. of times over all galaxy cutouts
    galname = row['gal']

    print(f"Fitting for {galname} ...")

    run_galfit(f'{galname}.fits', (pix_x[count], pix_y[count]), psf_file=psf_file)
    count += 1
    print(f"Fitting for {galname} done")

print("")
print(f"FITTING FOR {count} GALAXIES COMPLETED")
# """
### PART IV : Obtaining the parameters for each gal cutout 
# We will obtain the parameters using fetch_params()
# and we will prepare {galname}_result.fits with the median
# parameters for each galaxy

print("Preparing the model using median value parameters ...")
for row in table:
    galname = row['gal']
    dir = f"galfit_output/{galname}"

    # fetching all parameter arrays of particular galaxies
    pos_x, pos_y, int_mag, eff_rad, sersic, ax_ratio, PA = fetch_params(dir, 400)

    # plotting the distributions of the parameters in corresponding dirs
    color = "#4287f5"
    ec = "#4287f5"

    # creating a folder with galname inside plots directory
    # if the folder already exist (for example, in a second run), then don't create directory
    if subprocess.run(['ls', f"results/plots/{galname}"]).returncode != 0:
        # the returncode is 0 if the command succesfully 'ls' in the directory
        subprocess.run(["mkdir", f"results/plots/{galname}"])
	plot_dist(pos_x, color, ec, "Position-x", 'x', 'rho_x', f'results/plots/{galname}', 'pos_x')
	plot_dist(pos_y, color, ec, "Position-y", 'y', 'rho_y', f'results/plots/{galname}', 'pos_y')
    plot_dist(int_mag, color, ec, "Integrated magnitude", 'm', 'rho_m', f'results/plots/{galname}', 'int_mag')
    plot_dist(eff_rad, color, ec, "Effective radius", 'r', 'rho_r', f'results/plots/{galname}', 'eff_rad')
    plot_dist(sersic, color, ec, "Sersic index", 'n', 'rho_n', f'results/plots/{galname}', 'sersic')
    plot_dist(ax_ratio, color, ec, "Axis ratio", 'b/a', 'rho_(b/a)', f'results/plots/{galname}', 'ax_ratio')
    plot_dist(PA, color, ec, "Position angle", 'PA', 'rho_PA', f'results/plots/{galname}', 'PA')
    # plotting done

    # taking input from one of the galfit output file of the galaxy
    gal_inpt = open(file=f'{dir}/galfit.01', mode='r')

    gal_mod_inpt = open(file='galfit_mod.01', mode='w')

    # writing the params to file row-wise in different columns
    raw_params_file = open(f"results/raw_data/{galname}.txt", mode="w")
    for i in range(len(int_mag)):
        raw_params_file.write(f"{pos_x[i]}    {pos_y[i]}    {int_mag[i]}    {eff_rad[i]}    {sersic[i]}    {ax_ratio[i]}    {PA[i]}\n")

    # writing the median values to appropriate txt files for each galaxy
    params_file = open(f"results/data/{galname}.txt", mode='w')
    params_file.write(
        f'{np.median(pos_x)}  {np.median(pos_y)}  {np.median(int_mag)}  {np.median(eff_rad)}  {np.median(sersic)}  {np.median(ax_ratio)}  {np.median(PA)}'
        )
    params_file.close()

    # preparing model of a cutout using galfit -o2 to only create the model image

    print(f"Preparing model of {galname} ...")
    # notice that line starts from 0 here
    for i, line in enumerate(gal_inpt):
        if i == 7:
            line = f' A) ../cutouts/{galname}.fits          # Input data image (FITS file)\n'
        if i == 8:
            line = f' B) {galname}_result.fits         # Output data image block\n'
        if i == 10:
            line = f' D) ../{psf_file}.fits               # Input PSF image and (optional) diffusion kernel\n'
		if i == 47:
			line = f' 1) {np.median(pos_x)} {np.median(pos_y)}  1 1  #  Position x, y'
        if i == 48:
            line = f' 3) {np.median(int_mag)}     1          #  Integrated magnitude \n'
        if i == 49:
            line = f' 4) {np.median(eff_rad)}      1          #  R_e (effective radius)   [pix]\n'
        if i == 50:
            line = f' 5) {np.median(sersic)}      1          #  Sersic index n (de Vaucouleurs n=4)\n'
        if i == 54:
            line = f' 9) {np.median(ax_ratio)}      1          #  Axis ratio (b/a)\n'
        if i == 56:
            line = f'10) {np.median(PA)}      1          #  Position angle (PA) [deg: Up=0, Left=90]\n'
        
        gal_mod_inpt.write(str(line))

    gal_inpt.close()
    gal_mod_inpt.close()
    # estimating that the actual data is stored in cols 5 to 20


    subprocess.run(["galfit", "-o2","../galfit_mod.01"], cwd=f"results")
    
    print(f"Model for {galname} prepared\n")
    
print("All processes completed")
print(f"Time taken(s): {time.time() - init_time}")
