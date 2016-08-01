import json
import numpy as np

def Initialise():
    with open('../params.json') as param_file:    
		param = json.load(param_file)
    return param

def main():
    # Read all parameters from params.json file.
    param           = Initialise()
    
    # Define the data directory
    dat_directory   = param["directories"]["workdir"]

    W, F, E         = np.genfromtxt(dat_directory+"N_no2014_2016_07_27.txt"
                      ,unpack=True,skip_header=500,skip_footer= 500)

    f = open(dat_directory+'NI_no2014_owens.dat', 'w+')
    for i in range(len(W)):
        if not np.isnan(F[i]) and F[i] > 0:
            print >> f, " ","{: 1.10e}".format(W[i])," "+"{: 1.10e}".format(F[i])," "+"{: 1.10e}".format(E[i])
    f.close()
    
    print "done"

if __name__ == '__main__':
    main()
