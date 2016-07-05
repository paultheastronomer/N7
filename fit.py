import numpy as np
import matplotlib.pyplot as plt
import json

def Initialise():
    with open('params.json') as param_file:    
		param = json.load(param_file)
    return param

def BasicPlot(W,F):
    plt.step(W,F)
    plt.show()

def main():    
    param = Initialise()
    dat_directory   = param["directories"]["workdir"]

    W, F, E         = np.genfromtxt(dat_directory+param["files"]["datafile"],unpack=True)
    BasicPlot(W,F)



if __name__ == '__main__':
    main()

