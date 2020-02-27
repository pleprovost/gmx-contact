
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

if __name__ == "__main__":
    data = np.loadtxt('contact_matrix.dat')
    print(data.shape)

    plt.imshow(data, interpolation='nearest', cmap=plt.cm.viridis,\
               norm=colors.PowerNorm(gamma=0.5))
    plt.colorbar()
    plt.show()
    
