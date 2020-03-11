import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

if __name__ == "__main__":
    filename = None
    if len(sys.argv) < 1:
        filename = "contact_matrix.dat"
    else:
        filename = sys.argv[1]
        
    data = np.loadtxt(filename, comments=["#"])
    print(data.shape)

    sumx = []
    sumy = []
    for threshold in range(1000):
        sumx.append((np.sum(data, axis=0) > 100*threshold ).sum())
        sumy.append((np.sum(data, axis=1) > 100*threshold ).sum())
        
    print(' max : {} , min {} '.format(data.max(), data.min()))

    # plt.plot(sumx)
    # plt.plot(sumy)   
    
    plt.imshow(data, interpolation='nearest', cmap=plt.cm.viridis)
               #               norm=colors.PowerNorm(gamma=0.5))
    plt.colorbar()
    plt.show()
    
