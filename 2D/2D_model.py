import numpy as np
import scipy as sp
from scipy.special import jv,jvp,hankel1,h1vp, hankel1e
import cmath, time
import random

#create the geometry by randomly distributing circle inside a box. The position of the center, their radii and polar coordinates with respect to on edge of the box are returned
#"fill" is the filling factor, "a" the length of one size of the box, "mean_radius" the average radius and #std_radius# the standard deviation of the radius
def random_particles(fill,a, mean_radius, std_radius):
    
    # Initialize lists
    positions = [] #position of the centers of the rods
    rad = [] #radius of the rods
    
    #define the mean and std of the half gaussian to be the same as the values passed to the function
    std_half = np.sqrt((std_radius**2)/(1-2/np.pi))
    mean_half = mean_radius-np.sqrt(2*(std_half**2)/np.pi)
    
    # Generate random positions for particles until their total area occupies the desired filling factor
    while np.sum((np.pi*np.abs(rad)**2))/(a**2)<fill:
            rad.append(sp.stats.truncnorm.rvs(0, np.inf,loc=mean_half,scale=std_half))

    #places the cicles inside the box. Note that generating the circles while placing them deforms the distribution, that's why they were generated outside the loop
    for j in range(len(rad)):
        inte=0
        while True:
            # Generate random cartesian coordinates for the particle with origin at 0,0
            x = np.random.uniform(0, a)
            y = np.random.uniform(0, a)
            inte=inte+1
            
            # Check if particle overlaps with existing particles
            overlap = False

            #it is not guaranteed that the particle will fit in the box, so after too many tries we generate a new radius. This minimizes the distortion of the distribution
            if inte>10**6:
                rad[j]=sp.stats.halfnorm.rvs(loc=mean_half,scale=std_half)

            #condition for two circle to overlap
            for i in range(len(positions)):
                if np.linalg.norm(np.array(positions[i]) - np.array([x+1j*y])) < rad[i]+rad[j]:
                    overlap = True
                    break

            # If no overlap, add particle position to list
            if not overlap:
                positions.append([x+1j*y])
                break
                
    x, y = np.meshgrid(positions, positions) #define a vector with the positions
    z = x - y #define a vector with the relative distances
    r = np.abs(z) #radial distance between each fiber
    theta = np.angle(z) #angular distance between each fiber
    nrod = len(positions) #number of particles generated
    
    print("geometry done")
    
    return np.array(positions), np.array(rad), nrod, r, theta


#creates the T matrix for the cylinder (see eq. B5 of the paper)
#"k0" is the wave vector in vacuum, "radius" is a list with the radii of each rod, "epsilon" is a list with the electric constant of each rod, "nfour" is the number of harmoncis as from -nfour to nfour, "loss" is the amount of dielectric loss in units of conductivity  
def cylinder_T_matrix(i,k0, radius, epsilon, nfour,loss):
    
    lambd=2*np.pi/k0 #wavelength in vacuum
    P = 2 * nfour + 1 #total number of harmonics
    nu = np.sqrt(epsilon) + 1j*(loss/100)*lambd #refractive index

    d = np.arange(-nfour, nfour + 1) #list that runs through the harmonics
    ar1 = k0*radius #argument for the bessel functions outside the cylinder
    ar2 = k0*nu*radius.reshape(len(nu),1) #argument for the bessel functions inside the cylinder
    d2, ar2 = np.meshgrid(d, ar2) #the T matrix is matrix in harmonic space for each rod
    d1, ar1 = np.meshgrid(d, ar1) #the T matrix is matrix in harmonic space for each rod

    Dem = -nu*np.ones((1,P))*jvp(d2,ar2,1)*hankel1(d1,ar1) + jv(d2,ar2)*h1vp(d1,ar1,1) #denominator
    Num = nu*np.ones((1,P))*jv(d1,ar1)*jvp(d2,ar2,1)-jv(d2, ar2)*jvp(d1, ar1,1) #numerator

    #T matrix
    T = np.diag((Num/Dem).ravel())

    return T

#Axion source matrix (see eq B11 in the paper)
#"k0" is the wave vector in vacuum, "radius" is a list with the radii of each rod, "epsilon" is a list with the electric constant of each rod, "nfour" is the number of harmoncis as from -nfour to nfour 
def axion(k0, radius, epsilon, nfour):
    
    P = 2 * nfour + 1  #total number of harmonics
    nu = np.sqrt(epsilon) #refractive index (note that there is no loss term here)
    E0=((1/nu**2)-1) #induce electric field
    d = np.arange(-nfour, nfour + 1) #list that runs through the harmonics
    ar1 = k0*radius.reshape(len(radius),1) #argument for the bessel functions outside the cylinder
    ar2 = k0*nu*radius.reshape(len(nu),1) #argument for the bessel functions inside the cylinder
    
    Dem = -nu*jvp(0,ar2,1)*hankel1(0,ar1) + jv(0,ar2)*h1vp(0,ar1,1) #denominator
    Num = nu*E0*jvp(0,ar2,1) #numerator
    
    SSE = (Num/Dem).transpose() #axion source matrix

    #Kronecker  delta at the zeroth harmonic
    delta=np.zeros((2*nfour+1,1));
    delta[nfour]=1;

    SSE=delta*SSE #final axion source matrix

    #reshaping to be consistent with the T-matrix 
    SSE=SSE.transpose()
    SSE=SSE.ravel()
    
    return SSE

#Geometrical matrx coupling the rods (see eq B19 of the paper)
#"k0" is the wave vector in vacuum, "r" is a list with the radial distance between each fiber, "theta" is a list with the radial distance between each fiber, "nfour" is the number of harmoncis as from -nfour to nfour, "nrod" is the number of rods
def coupling_matrix(k0, r, theta, nrod, nfour):

    if ntige == 1: #no coupling matrix if only 1 rod
        s = []
    else:
        P = 2*nfour+1  #total number of harmonics
        nmat = nrod*(nrod-1)//2 #size of the matrix, pairing the rods
        inx=np.zeros(nmat,int) #index runs through the size of the matrix
        iny=np.zeros(nmat,int) #index runs through the size of the matrix
        
        sp1 = np.zeros((nmat,2*nfour+1,2*nfour+1),complex) #matrix size is number of pairs per harmonic per harmonic
        sp2 = np.zeros((nmat,2*nfour+1,2*nfour+1),complex) #matrix size is number of pairs per harmonic per harmonic
        spp = np.zeros((nrod,nrod,2*nfour+1,2*nfour+1),complex) #matrix size is number of pairs per harmonic per harmonic, will be the G matrix after reshaping

        k, l = np.meshgrid(np.arange(-nfour, nfour + 1), np.arange(-nfour, nfour + 1)) #contain the harmonics pairs
        ind=0 #will run through the matrix

        #fill the paired indeces to correctly account for their respectave rod, as given by the r and theta in the fucntion defined above
        for m in range(nrod):
            for nn in range(m+1,nrod):
                ind += 1
                inx[ind-1] = m
                iny[ind-1] = nn

        #runs over the pairs
        for m in range(nmat):
            ans = np.exp(1j*(k-l)*theta[inx[m],iny[m]])*hankel1(l-k, k0*r[inx[m],iny[m]])
            sp1[m] = ans #matrix for the pair in the position inx[m],iny[m]
            sp2[m] =(-1.0)**(k-l)*ans #matrix for the pair in the position iny[m],iny[m] this is the transpose

        #store the values in the spp matrix
        for n in range(len(inx)):
            spp[inx[n],iny[n]] = sp1[n]
            spp[iny[n],inx[n]] = sp2[n]
            
        #the diagonals are zero
        for n in range(nrod):
            spp[n,n]=np.zeros((P,P))
            
        #reshape to be consistent with the T matrix and the Axion source matrix
        G_matrix=np.vstack([np.hstack(c) for c in spp])
        
    return G_matrix


#function to compute the electric field inside the cylinders given the solution of the scattering problem (see eq B17 in the paper)
#"ycarte" is a list with the [x+j*y] positions we want to compute the electric field, "D_coeff" are the solutions of the system,"affixe" is a list with the position of the center of each fiber, "k0" is the wave vector in vacuum, "radius" is a list with the radii of each rod, "nfour" is the number of harmoncis as from -nfour to nfour, "nrod" is the number of rods
def mapfield(ycarte, radius, k0, affixe, nfour, nrod, D_coeff):
    mf = np.arange(-nfour, nfour + 1) #run through the harmonics
    m, n = ycarte.shape #run through the points we want the electric field at
    ZE = np.empty((m, n), dtype=np.complex128) #will store the electric field

    #fill the list
    for l in range(m):
        for p in range(n):
            R = []
            indice = np.where(np.abs(ycarte[l, p] - affixe) - radius.reshape(nrod,1) < 0)[0] #pick the values outside the cylinder
            if len(indice) == 0:
                pc = ycarte[l, p] - affixe #points we want the electric field at with the rods as the origin
                mmf, pc = np.meshgrid(mf, pc)
                R = (hankel1(mmf, k0 * np.abs(pc)) * np.exp(1j * np.angle(pc) * mmf)).flatten() #matrix of the hankel functions in harmonic space
                ZE[l, p] = np.dot(R, D_coeff) #electric field

    return ZE

#integrate the poynting vector over a detector that occupies the whole angular domain encompasing the configuration
#"a" is the box size, "steptheta" and "stepr" are the discretization of the angular/radial coordinate for integration, "L" is the distance of the detector to the center of the box
def integrate_over_angle(a, steptheta, stepr, L, radius, k0, epsilon, affixe, nfour, nrod, D_coeff):

    #to find the magnetic field, we take a derivative of the electric field, so we compute it at L+dr and L-dr then substract
    first_circle = (L-stepr)*np.cos(np.arange(0,2*np.pi,steptheta)) + 1j*( (L-stepr)*np.sin(np.arange(0,2*np.pi,steptheta))) #cartesian coordinates at L-dr
    second_circle= (L+stepr)*np.cos(np.arange(0,2*np.pi,steptheta)) + 1j*( (L+stepr)*np.sin(np.arange(0,2*np.pi,steptheta))) #cartesian coordinates at L-dr
    
    ppoyn=np.array([first_circle,second_circle]) #array with all the points we want to compute the electric field at
    ZpE = mapfield(ppoyn, radius, k0, affixe, nfour, nrod, D_coeff) #computes the electric field at the points above
    uE = np.diff(ZpE, axis=0) / (2 * stepr) #find the magnetic field by taking a numerical derivative
    FluxPE = 1j * ZpE * np.conj(uE) / (2 * k0) #compute the poynting vector
    PE = L*np.sum(FluxPE[0, :] * steptheta) #numerical integration of the poynting vector over the whole angular domain
    
    return PE

#function that computes the power emitted by the configuration, build the matrix and call integrate_over_angle() to find the power
#"radius0" is the mean radius of the particles, "sigma" is the standard deviation of particles, "a" is the box size, "fill" is the filling factor from 0 to 1, "lam" is a list with all the wavelengths we want to compute the power at, "L" is the radius of the detector that occupies the whole angular domain
def power(radius0, sigma, a, fill, lam, L):
    
    nfour = 5 #number of harmonics
    
    stepr = 2e-10 #step to compute the magnetic field at integrate_over_angle()
    steptheta=np.pi/20 #step to integrate the poynting vector at integrate_over_angle()

    #builds the geometry and save the position, radius, number of rods, and relative radial/angular position
    affixe, radius, nrod, r, theta = random_particles(fill, a, radius0, sigma)

    #all cylinders will have the same dielectric constant in our case
    epsilon0 = 1.77 * 1.77 #dielectric constant inside the cylinder
    epsilon = epsilon0 * np.ones((nrod, 1)) #all rods have the same dielectric constant
    
    loss=0 #dielectric loss as defined in the paper
    
    poynting = np.ones_like(lam) #intialize list

    #for each wavelength, build the matrices and solves the system
    for i in range(len(lam)):
        
        k0 = 2 * np.pi / lam[i] #wave number in vacuum

        #build the G and T matrix
        G = coupling_matrix(k0, r, theta, nrod, nfour) #G in the paper
        T_matrix = cylinder_T_matrix(i,k0, radius, epsilon, nfour,loss) #T in the paper

        #matrix that will be inverted in eq B25 of the paper
        if T_matrix.size == 0: #in case there are no cylinders we get identity
            HE = np.eye(T_matrix.shape[0])
        else:
            HE = np.eye(G.shape[0]) - T_matrix @ G #build the matrix required for inversion 

        #axion source matrix
        S = axion(k0, radius, epsilon, nfour)

        #solves the linear system in eq B25 of the paper
        D_coeff = np.linalg.solve(HE, S) #D in the paper

        #finds the power emitted by the configuration for these wavelength
        PE = integrate_over_angle(a, steptheta, stepr, L, radius, k0, epsilon, affixe, nfour, nrod, D_coeff)
        
        poynting[i] = np.abs(PE) #stores the power
        
    return poynting


#worker function, ready to be parallelized in the lam loop
def worker(lam):
    
    L = 100 #detector radius
    radius0 = 1 #average radius of the rods
    sigma = 0.5 #standard deviation of the rods
    a = 20*np.sqrt(10) #size of the box, just make sure it is reasonable with the filling factor, if it is too large the code becomes inefficient
    fill = 0.5 #filling factor
    
    poynting = power(radius0, sigma, a, fill, lam, D, L) #call the function for a particular wavelength
    
    return poynting

def main():
    
    lam = np.logspace(-2,1,1000) #wavelengths we want to probe

    iterations=np.arange(0,100,1) #how many times we want to run for each lambda, the result of the figures are averaged over 100 runs
    data=np.zeros((100,len(lam))) #store the final data

    #loop in the iterations
    for i in iterations:
        data[i]=worker(lam)

    #saves the data
    np.save("data_test.npy",data)

main()
