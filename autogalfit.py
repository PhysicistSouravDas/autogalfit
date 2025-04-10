"""
.. module:: autogalfit

:Synopsis: Run GALFIT multiple times over a field and compute median parameter values.
:Author: Sourav Das
:Year: 2023

This module is to run galfit for a field multiple times and 
obtain median value of the parameters
The idea is to run GALFIT every time, which will create two
types of output txt files, and the result.fits file
which will be stored in corresponding folders with gal names
then the median values of the parameters will be returned
by the function
"""

import os, random, subprocess

def run_galfit(cutout: str, pix_val: tuple, psf_file: str='psfacs00', n: int=400):
    """
    Runs galfit n no. of times on a cutout and obtain the median
    values of the parameters
    cutout(str): Name of the cutout file
    pix_val(tuple): (xpix, ypix) value of the galaxy in the cutout
    n(int): No. of times to iterate (default 400)

    """
    """
    Integrated magnitude (line 47), Effective radius (48) 
    Sersic index (49), Axis ratio (53), Position angle (54)

    Range of values:
    Integrated magnitude: (15, 30)
    Effective radius: (2, 20)
    Sersic index: (0.5, 5)
    Axis ratio: (0.1, 1)
    Position angle: (-90, 90)
    """


    # getting current directory
    curr_dir = os.path.abspath(os.getcwd())
    out_dir = f'{curr_dir}/galfit_output/'
    subprocess.run(["mkdir", f"galfit_output/{cutout[:-5]}"])
    ITER = n
    num_runs = 0
    actualrun = 0


    # for _ in range(ITER):
    while (num_runs < ITER):
    	actualrun += 1
        # To scale random(0,1) between (a,b): => random(a,b) = a + (b-a)*random(0,1)

        # choosing random numbers for the parameters
        int_mag = 15 + (30-15)*random.random()
        eff_rad = 2 + (20-2)*random.random()
        sersic = 0.5 + (5-0.5)*random.random()
        ax_ratio = 0.1 + (1-0.1)*random.random()
        PA = -90 + (180)*random.random()
        
        file_names = []
        for filename in os.listdir(f"galfit_output/{cutout[:-5]}/"):
        	if 'galfit' in filename: file_names.append(filename)
        
        file_runs = [i.split('.')[-1] for i in file_names]
        file_runs = [int(i) for i in file_runs]
        num_runs = max(file_runs)
        	

        gal_inpt = open(file='galfit_input.feedme', mode='r')

        gal_mod_inpt = open(file='galfit_mod.feedme', mode='w')

        # notice that line starts from 0 here
        for i, line in enumerate(gal_inpt):
            if i == 5:
                line = f' A) ../../cutouts/{cutout}          # Input data image (FITS file)\n'
            if i == 8:
                line = f' D) ../../{psf_file}.fits               # Input PSF image and (optional) diffusion kernel\n'
            if i == 45:
                line = f' 1) {pix_val[0]}  {pix_val[1]}  1 1  #  Position x, y\n'
            if i == 46:
                line = f' 3) {int_mag}     1          #  Integrated magnitude \n'
            if i == 47:
                line = f' 4) {eff_rad}      1          #  R_e (effective radius)   [pix]\n'
            if i == 48:
                line = f' 5) {sersic}      1          #  Sersic index n (de Vaucouleurs n=4)\n'
            if i == 52:
                line = f' 9) {ax_ratio}      1          #  Axis ratio (b/a)\n'
            if i == 53:
                line = f'10) {PA}      1          #  Position angle (PA) [deg: Up=0, Left=90]\n'
            
            gal_mod_inpt.write(str(line))

        gal_inpt.close()
        gal_mod_inpt.close()
        # estimating that the actual data is stored in cols 5 to 20


        while True :
            errcode = subprocess.run(["galfit", "../../galfit_mod.feedme"], cwd=f"galfit_output/{cutout[:-5]}").returncode
            # the above code executes the command and stores the returncode
            if errcode == 0: break
        # this command runs galfit from './galfit_output/cutout'
        # {cutout[:-5]} trims .fits from gal_name.fits and only retains gal_name 
        
        logfile = open(f'galfit_output/{cutout[:-5]}/fit.log')
        for i, logline in enumerate:
        	if i == 15*(actualrun-1) + (10-1):
        		# line containing the converged or diverged params
        		if '*' in logline:
        			# then delete the last file
        			if num_runs <= 8:
        				os.remove(f'galfit_output/{cutout[:-5]}/galfit.0{num_runs}')
        			else:
        				os.remove(f'galfit_output/{cutout[:-5]}/galfit.{num_runs}')
		logfile.close()
        				
        


def fetch_params(dir_name: str, num_files: int):
    """
    Fetch the required parameters from the galfit output txt files
    dir_name(str): Name of directory in which many galfit txt files stored
    num_files(int): Number of galfit txt files in the directory

    Returns array of the parameters
    """
    # the name of galfit output files are like galfit.01, galfit.02 and so on
    # in the output files obtained via a small test run, the target parameters
    # are located on the following lines
    # int_mag: 49
    # eff_rad: 50
    # sersic:  51
    # ax_ratio:55
    # PA:      56
    pos_x = []
    pos_y = []
    int_mag = []
    eff_rad = []
    sersic = []
    ax_ratio = []
    PA = []

    # galfit file names until 99 files are like galfit.01, ..., galfit.99, then,
    # galfit.100, galfit.101, ..., galfit.400

    file_list = os.listdir(f'{dir_name}')

    for filename in file_list:
    # for i in range(num_files):
    #     i += 1  # since galfit files start from 01 instead of 00
    #     if i < 10:
    #         f = open(f'{dir_name}/galfit.0{i}', 'r')
    #     else:
    #         f = open(f'{dir_name}/galfit.{i}', 'r')

        if 'galfit' in filename:
            f = open(f'{dir_name}/{filename}', 'r')
        
            for j, line in enumerate(f):
                # note that j starts from 0 instead of 1
                # so line 1 becomes line 0 here
                if j == 47:
                	pos_x.append(float(line.split()[1]))
                	pos_y.append(float(line.split()[2]))
                if j == 48:
                    int_mag.append(float(line.split()[1]))
                if j == 49:
                    eff_rad.append(float(line.split()[1]))
                if j == 50:
                    sersic.append(float(line.split()[1]))
                if j == 54:
                    ax_ratio.append(float(line.split()[1]))
                if j == 55:
                    PA.append(float(line.split()[1]))
            
            f.close()

    return pos_x, pos_y, int_mag, eff_rad, sersic, ax_ratio, PA
