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
    r = np.abs(z) #radial position in polar coordinates
    theta = np.angle(z) #angular position in polar coordinates
    nrod = len(positions) #number of particles generated
    
    print("geometry done")
    
    return np.array(positions), np.array(rad), nrod, r, theta


def cylinder_T_matrix(i,k0, radius, epsilon, nfour,loss):
    
    lambd=2*np.pi/k0
    P = 2 * nfour + 1
    nu = np.sqrt(epsilon) + 1j*(loss/100)*lambd #loss is included

    d = np.arange(-nfour, nfour + 1)
    ar1 = k0*radius
    ar2 = k0*nu*radius.reshape(len(nu),1)
    d2, ar2 = np.meshgrid(d, ar2)
    d1, ar1 = np.meshgrid(d, ar1)

    N = -nu*np.ones((1,P))*jvp(d2,ar2,1)*hankel1(d1,ar1) + jv(d2,ar2)*h1vp(d1,ar1,1)    
    D = nu*np.ones((1,P))*jv(d1,ar1)*jvp(d2,ar2,1)-jv(d2, ar2)*jvp(d1, ar1,1)
    
    T = np.diag((D/N).ravel())

    return T


def axion(k0, radius, epsilon, nfour):
    
    P = 2 * nfour + 1
    nu = np.sqrt(epsilon)
    E0=((1/nu**2)-1)
    d = np.arange(-nfour, nfour + 1)
    ar1 = k0*radius.reshape(len(radius),1)
    ar2 = k0*nu*radius.reshape(len(nu),1)
    
    D = -nu*jvp(0,ar2,1)*hankel1(0,ar1) + jv(0,ar2)*h1vp(0,ar1,1) 
    N = nu*E0*jvp(0,ar2,1)
    
    SSE = (N/D).transpose()

    delta=np.zeros((2*nfour+1,1));
    delta[nfour]=1;

    SSE=delta*SSE
    SSE=SSE.transpose()
    
    SSE=SSE.ravel()
    
    return SSE


def coupling_matrix(k0, r, theta, ntige, nfour):

    if ntige == 1:
        s = []
    else:
        P = 2*nfour+1
        nmat = ntige*(ntige-1)//2
        inx=np.zeros(nmat,int)
        iny=np.zeros(nmat,int)
        
        sp1 = np.zeros((nmat,2*nfour+1,2*nfour+1),complex)
        spp = np.zeros((ntige,ntige,2*nfour+1,2*nfour+1),complex)
        sp2 = np.zeros((nmat,2*nfour+1,2*nfour+1),complex)
        k, l = np.meshgrid(np.arange(-nfour, nfour + 1), np.arange(-nfour, nfour + 1))
        ind=0
        
        for m in range(ntige):
            for nn in range(m+1,ntige):
                ind += 1
                inx[ind-1] = m
                iny[ind-1] = nn

        for m in range(nmat):
            ans = np.exp(1j*(k-l)*theta[inx[m],iny[m]])*hankel1(l-k, k0*r[inx[m],iny[m]])
            sp1[m] = ans
            sp2[m] =(-1.0)**(k-l)*ans

        for n in range(len(inx)):
            spp[inx[n],iny[n]] = sp1[n]
            spp[iny[n],inx[n]] = sp2[n]
            

        for n in range(ntige):
            spp[n,n]=np.zeros((P,P))
            
        
        s=np.vstack([np.hstack(c) for c in spp])
        
    return s


def mapfield(ycarte, radius, k0, affixe, nfour, nrod, BE):
    mf = np.arange(-nfour, nfour + 1)
    m, n = ycarte.shape
    ZE = np.empty((m, n), dtype=np.complex128)
    for l in range(m):
        for p in range(n):
            R = []
            indice = np.where(np.abs(ycarte[l, p] - affixe) - radius.reshape(nrod,1) < 0)[0]
            if len(indice) == 0:
                pc = ycarte[l, p] - affixe
                mmf, pc = np.meshgrid(mf, pc)
                R = (hankel1(mmf, k0 * np.abs(pc)) * np.exp(1j * np.angle(pc) * mmf)).flatten()
                ZE[l, p] = np.dot(R, BE)

    return ZE


def integrate_over_det(a, xmin, xmax, stepx, stepy, gpoyn, radius, k0, epsilon, affixe, nfour, nrod, BE):
    xpoyn, ypoyn = np.meshgrid(np.arange(xmin, xmax + stepx, stepx), [-stepy, stepy])
    ppoyn = xpoyn + a / 2 + 1j * ypoyn + gpoyn
    ZpE = mapfield(ppoyn, radius, k0, affixe, nfour, nrod, BE)
    uE = np.diff(ZpE, axis=0) / (2 * stepy)
    FluxPE = 1j * ZpE * np.conj(uE) / (2 * k0)
    PE = np.sum(FluxPE[0, :] * stepx) / (xmax - xmin)
    return PE

def integrate_over_angle(a, steptheta, stepr, L, radius, k0, epsilon, affixe, nfour, nrod, BE):
    
    first_circle = (L-stepr)*np.cos(np.arange(0,2*np.pi,steptheta)) + 1j*( (L-stepr)*np.sin(np.arange(0,2*np.pi,steptheta)))
    second_circle= (L+stepr)*np.cos(np.arange(0,2*np.pi,steptheta)) + 1j*( (L+stepr)*np.sin(np.arange(0,2*np.pi,steptheta)))
    
    ppoyn=np.array([first_circle,second_circle])
    ZpE = mapfield(ppoyn, radius, k0, affixe, nfour, nrod, BE)
    uE = np.diff(ZpE, axis=0) / (2 * stepr)
    FluxPE = 1j * ZpE * np.conj(uE) / (2 * k0)
    PE = L*np.sum(FluxPE[0, :] * steptheta)
    
    return PE


def power(radius0, sigma, a, fill, lam, D, L):
    
    nfour = 5 #number of harmonics
    
    xmin = -D #size of the dector
    xmax = D
    gpoyn = -L * 1j 
    stepx = D / 10 
    stepy = 2e-10 
    
    stepr = 2e-10 
    steptheta=np.pi/20 

    affixe, radius, nrod, r, theta = random_particles(fill, a, radius0, sigma)
    
    epsilon0 = 1.77 * 1.77 #dielectric epsilon
    epsilon = epsilon0 * np.ones((nrod, 1)) #all rods have the same dielectric constant
    loss=0 #dielectric loss as defined in the paper
    
    poynting = np.ones_like(lam)
    
    for i in range(len(lam)):
        
        k0 = 2 * np.pi / lam[i]
        if nfour == 0:
            rreg = r + np.eye(r.shape[0])
            G = hankel1(0, k0 * np.sqrt(epsilon0) * rreg)
            np.fill_diagonal(G, 0)
        else:
            G = coupling_matrix(k0, r, theta, nrod, nfour) #G in the paper

        T_matrix = cylinder_T_matrix(i,k0, radius, epsilon, nfour,loss) #T in the paper

        if T_matrix.size == 0:
            HE = np.eye(T_matrix.shape[0])
        else:
            HE = np.eye(G.shape[0]) - T_matrix @ G

        S = axion(k0, radius, epsilon, nfour)
        
        D_coeff = np.linalg.solve(HE, S) #D in the paper
        
        PE = integrate_over_angle(a, steptheta, stepr, L, radius, k0, epsilon, affixe, nfour, nrod, D_coeff)
        
        poynting[i] = np.abs(PE)
        
        if i%10==0:
            print(i)
        
    return poynting


def worker(lam):
    
    L = 100
    D = 50
    radius0 = 1
    sigma = 0.1
    a = 20*np.sqrt(10) 
    fill = 0.5
    
    poynting = power(radius0, sigma, a, fill, lam, D, L)
    
    return poynting

def main():
    
    lam = np.logspace(-2,1,100)
    
    data=np.zeros((100,len(lam)))
    iterations=np.arange(0,5,1)
    
    for i in iterations:
        data[i]=worker(lam)
        print("one done")
    
    np.save("data_test.npy",data)

main()
