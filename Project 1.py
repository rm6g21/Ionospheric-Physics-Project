import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import i0 
from scipy.integrate import simps
from PIL import Image

#import brightness images --> Im sorry if you run this code and it doesnt give you access
#im not sure how to change that. sorry
#if code doesnt work just ignore the images and remove them from brightness plot, they arent needed just nice visuals
img1 = np.asarray(Image.open('c:\\Users\\robbi\\OneDrive\\Pictures\\12_seconds_barium.PNG'))
img2 = np.asarray(Image.open('c:\\Users\\robbi\\OneDrive\\Pictures\\27_seconds_barium.PNG'))
img3 = np.asarray(Image.open('c:\\Users\\robbi\\OneDrive\\Pictures\\64_seconds_barium.PNG'))

def Brightness(C,ds):


   

    lim = len(C[:][0][0])
    B = np.zeros((lim,lim))
    #for now only do positive x and z
    for n_x in range(lim):
        for n_z in range(lim):
            for n_y in range(lim):
                B[n_x][n_z] += 2*ds*C[n_x][n_y][n_z]
                
    # Create a new array that is twice the size of B in each dimension
    B_full = np.zeros((2*lim, 2*lim))

    # Fill the new array with the values of B in each quadrant
    B_full[lim:, lim:] = B  # B(x, y)
    B_full[:lim, lim:] = np.flip(B, axis=0)  # B(-x, y)
    B_full[lim:, :lim] = np.flip(B, axis=1)  # B(x, -y)
    B_full[:lim, :lim] = np.flip(B)  # B(-x, -y)

    # Plot the new array
    """
    plt.imshow(B_full, cmap='PuBu', interpolation='nearest')
    plt.xlabel('x / length')
    plt.ylabel('z / length')
    plt.title('Relative brightness when looking up at sphere at n_t=4')
    plt.colorbar()
    """
    return B_full
    
def theoretical_Brightness(d_x,d_t,R,max_dist,n_t):
    
    #For some reason numpy wouldn't recognise zeros here so we import it directly
    from numpy import zeros

    
    B = zeros((max_dist,max_dist))

    #Define D

    D = d_x**2 / (6*d_t)

    t = n_t*d_t

  
    

    

    
    
    
    #loop over to find B(x,y) rather than B(rho)
    #Harder to visualise on 2D grid when only in terms of rho
    for n_x in range(max_dist):
        for n_y in range(max_dist):
            x = n_x * d_x
            y = n_y * d_x
            rho = np.sqrt(x**2 + y**2)

            #code below for integration using Simspon's rule
            def integrand(theta):
                return np.exp(((-R**2 * (np.sin(theta))**2)/(4*D*t))) *i0(R*rho*np.sin(theta)/(2*D*t))*np.sin(theta)*(np.cos(theta)**2)
            
            range_list = np.linspace(0,np.pi/2,1000)

            y = integrand(range_list)
    #Calculate the integral
            B[n_x][n_y] = (R**3/(D*t))*np.exp((-(rho)**2/(4*D*t)))*simps(y,range_list,axis = -1)
    
    #Make a full circle
    B_full = np.zeros((2*max_dist, 2*max_dist))
    B_full[max_dist:, max_dist:] = B  # B(x, y)
    B_full[:max_dist, max_dist:] = np.flip(B, axis=0)  # B(-x, y)
    B_full[max_dist:, :max_dist] = np.flip(B, axis=1)  # B(x, -y)
    B_full[:max_dist, :max_dist] = np.flip(B)  # B(-x, -y)
    return B_full





   

def new_erfc(x):
    if x<0:
        return 2-math.erfc(-x)
    else :
        return math.erfc(x)

def theoretical_density(R, max_dist, max_t, dt, dr):
    #Calculate D using relation
    D = dr**2 / (6*dt)

    C = np.zeros((max_dist,max_t),dtype = float)
    for n_t in range(max_t):
        for n_r in range(max_dist):
            #define r and t
            r = (n_r * dr) +10e-30
            t = (n_t * dt)+10e-30
            
            #Define r_plus and r_minus
            r_plus = (R+r)/(2*np.sqrt(D*t))
            r_minus = (R-r)/(2*np.sqrt(D*t))

            #Calculate C
            
            C[n_r][n_t] = 0.5*(2-new_erfc(r_plus) - new_erfc(r_minus))+((1/r)*np.sqrt(D*t/np.pi))*(np.exp(-(r_plus**2)) - np.exp(-(r_minus**2)))
           
    return C

#Define 3D array as set of zeros

#We are measuring 10 minutes of diffusion with 1 second intervals
time_steps, dt = np.linspace(0,200,200,retstep=True)
#We are measuring across 300m of space with 1m intervals
space_steps, ds = np.linspace(0,6000,50,retstep=True)
max_t = len(time_steps)
max_dist = len(space_steps)








#make a 3D array: cubic spatial
C = np.zeros((max_dist,max_dist,max_dist),dtype=float)





#initial conditions
# Set radius for initial conditions
initial_radius = max_dist / 2




#initial number density
N0 = 3.89e13

# Set initial conditions
for n_x in range(max_dist):
    for n_y in range(max_dist):
        for n_z in range(max_dist):
            # Check if point is within the initial radius
            if np.sqrt(n_x**2 + n_y**2 + n_z**2) < initial_radius:
                # Set initial conditions with Gaussian decay
                C[n_x][n_y][n_z] = 1
            elif initial_radius -0.5 <=np.sqrt(n_x**2 + n_y**2 + n_z**2) < initial_radius + 0.5:
                C[n_x][n_y][n_z] = 0.5
                




newC = C[:][:][:]


#error_without_boundary = [0.041362931081663855, 0.039615176343294284, 0.04181115470594442, 0.049983722129009894, 0.05567421592429207, 0.05946333423543956, 0.0617756170736944, 0.06291914507740265, 0.06324789420195438, 0.06302255328748958, 0.06242755513907333, 0.06159166827324028, 0.06231960725465234, 0.06291559553713662, 0.0632484998894215, 0.0633717444245428, 0.06332695307183994, 0.06314685059813405, 0.06285740119982122, 0.062479384603944754, 0.062029564199310536, 0.061521560323752306, 0.060966509731398694, 0.0603735688660888, 0.060262174000482144, 0.060333917911331025, 0.060337516904713395, 0.06028110549128912, 0.06017181859321692, 0.060015925060654396, 0.05981894193259911, 0.05958573238013634, 0.05932058981339704, 0.05902731024299021, 0.05870925465613394, 0.05836940288963946, 0.058010400248735966, 0.05763459792551809, 0.057244088107456084, 0.05684073452966447, 0.05679074270140694, 0.05672395219937387, 0.05663347820574563, 0.0565212388695564, 0.05638899448560612, 0.0562383615708367, 0.056070825610634326, 0.055887752601360405, 0.055690399504747995, 0.055479923719553695, 0.05525739166617437, 0.05502378657090903, 0.054780015528200124, 0.05452691591153362, 0.0542652611966869, 0.05399576625466479, 0.05371909216591603, 0.053435850602229135, 0.053164439619856174, 0.05308155952570755, 0.05298750186146728, 0.052882918161405876, 0.05276842160879973, 0.052644589495936894, 0.05251196551566914, 0.05237106189614192, 0.052222361389636554, 0.05206631912576985, 0.05190336433862608, 0.051733901976759714, 0.0515583142043924, 0.05137696180154826, 0.051190185470326416, 0.050998307053994496, 0.05080163067510213, 0.050600443798369754, 0.05039501822368006, 0.05018561101411617, 0.04997246536362149, 0.04975581140852364, 0.049577192849266496, 0.04948669392788446, 0.04939052086695796, 0.04928892363734952, 0.049182141246073355, 0.04907040225957199, 0.04895392530023138, 0.04883291951752114, 0.04870758503508421, 0.04857811337503888, 0.04844468786069078, 0.04830748399879917, 0.04816666984247761, 0.04802240633576017, 0.04787484764080517, 0.04772414144866123, 0.04757042927447093, 0.04741384673793812, 0.0472545238298446, 0.047092585165355694, 0.04692815022481603, 0.04676133358269907, 0.04659224512533569, 0.04642099025801566, 0.04624767010201894, 0.04616100482251441, 0.0460724862557824, 0.0459807789330676, 0.045885992985592164, 0.045788234781092604, 0.04568760706232035, 0.04558420908008641, 0.04547813672106911, 0.04536948263059379, 0.04525833633058919, 0.04514478433291709, 0.04502891024825978, 0.04491079489075117, 0.044790516378521966, 0.04466815023032859, 0.04454376945842571, 0.044417444657836355, 0.044289244092168265, 0.04415923377611801, 0.044027477554799785, 0.04389403718002792, 0.04375897238368008, 0.04362234094826023, 0.043484198774774774, 0.043344599948035356, 0.04320359679948965, 0.043086480176903884, 0.043006526114709724, 0.04292452045333968, 0.042840519433510506, 0.04275457773122504, 0.042666748504168975, 0.04257708343666755, 0.04248563278324693, 0.04239244541084495, 0.04229756883971107, 0.042201049283039076, 0.0421029316853697, 0.04200325975980466, 0.04190207602406855, 0.04179942183545372, 0.04169533742468742, 0.04158986192875054, 0.04148303342268726, 0.041374888950431524, 0.04126546455468763, 0.041154795305891584, 0.04104291533028333, 0.040929857837120276, 0.040815655145056445, 0.04070033870771746, 0.04058393913849449, 0.04046648623458296, 0.04034800900029212, 0.04022853566964564, 0.04013627477705772, 0.040061576268343745, 0.03998548016287793, 0.03990801661200985, 0.03982921507646503, 0.039749104343598254, 0.03966771254420516, 0.039585067168900624, 0.039501195084076304, 0.03941612254744702, 0.039329875223196295, 0.039242478196730574, 0.03915395598905416, 0.0390643325707708, 0.03897363137572459, 0.03888187531428858, 0.03878908678630891, 0.038695287693714414, 0.038600499452800384, 0.038504743006194615, 0.03840803883451323, 0.03831040696771645, 0.03821186699616896, 0.038112438081416365, 0.038012138966680924, 0.03791098798708879, 0.03780900307963017, 0.03770620179286495, 0.03760260129637687, 0.03749821838998342, 0.03739306951270936, 0.03729120658668592, 0.037219295376618194, 0.03714635387431986, 0.037072399673146364, 0.03699745004278389, 0.036921521936128486, 0.036844631996013094, 0.036766796561783464, 0.03668803167572548, 0.03660835308935014, 0.03652777626953605, 0.03644631640453381, 0.03636398840983614, 0.036280806933915485, 0.036196786363831966, 0.036111940830715505, 0.03602628421512441, 0.03593983015228282]






"""
#time evoluiton
fig, axs = plt.subplots(2,2)
axs[0,0].plot(C[:][0][0],ls='-',label='Finite Difference')
axs[0,0].plot(C_theoretical[:,0],ls='--',label='Theoretical Solution')
"""

#Plotting the theoretical brightness using crude animation
fig = plt.figure()

C_theoretical = theoretical_density(initial_radius*ds,max_dist,max_t,dt,ds)

    
error_list = []
arg_list = []
value_list = []

central_B = [Brightness(C,ds)[50,50]]
B_error = [Brightness(C,ds)[50,50]-theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,0)[50,50]]

radius_list = [25*ds]
radius_error = [60]

fig, axs = plt.subplots(2,2)

for n_t in range(1,20): 
    #qualitative tracker for how long code is running for
    print(n_t)

    ##newC[max_dist-1][:][:]= 0.5*C[max_dist-2][:][:] 
    #newC[:][max_dist-1][:]= 0.5*C[:][max_dist-2][:] 
    #newC[:][:][max_dist-1]= 0.5*C[:][:][max_dist-2]
    

    for n_x in range(0,max_dist-1):
        for n_y in range(0,max_dist-1):
            for n_z in range(0,max_dist-1):
                if n_x == 0 and n_y>0 and n_z>0:  # Check if the point is on the boundary
                    # Use different time evolution equation for boundary points
                    newC[n_x][n_y][n_z] =  (1/6)*(2*C[n_x+1][n_y][n_z] + C[n_x][n_y+1][n_z] + C[n_x][n_y-1][n_z] + C[n_x][n_y][n_z+1] + C[n_x][n_y][n_z-1])
                elif n_y == 0 and n_x>0 and n_z>0:
                    newC[n_x][n_y][n_z] =  (1/6)*(C[n_x+1][n_y][n_z] + C[n_x-1][n_y][n_z] + 2*C[n_x][n_y+1][n_z] + C[n_x][n_y][n_z+1] + C[n_x][n_y][n_z-1])
                elif n_z == 0 and n_x>0 and n_y>0:
                    newC[n_x][n_y][n_z] =  (1/6)*(C[n_x+1][n_y][n_z] + C[n_x-1][n_y][n_z] + C[n_x][n_y+1][n_z] + C[n_x][n_y-1][n_z] + 2*C[n_x][n_y][n_z+1])
                elif n_x==0 and n_y == 0 and n_z>0:
                    newC[n_x][n_y][n_z] =  (1/6)*(2*C[n_x+1][n_y][n_z] + 2*C[n_x][n_y+1][n_z] + C[n_x][n_y][n_z+1] + C[n_x][n_y][n_z-1])
                elif n_x==0 and n_z == 0 and n_y>0:
                    newC[n_x][n_y][n_z] =  (1/6)*(2*C[n_x+1][n_y][n_z] + C[n_x][n_y+1][n_z] + C[n_x][n_y-1][n_z] + 2*C[n_x][n_y][n_z+1])
                elif n_y==0 and n_z == 0 and n_x>0:
                    newC[n_x][n_y][n_z] =  (1/6)*(C[n_x+1][n_y][n_z] + C[n_x-1][n_y][n_z] + 2*C[n_x][n_y+1][n_z] + 2*C[n_x][n_y][n_z+1])
                elif n_x==0 and n_y == 0 and n_z == 0:
                    newC[n_x][n_y][n_z] =  (1/3)*(C[n_x+1][n_y][n_z] + C[n_x][n_y+1][n_z] + C[n_x][n_y][n_z+1])
                
     
                
                #boundary condition
               
                
                else:
                    # Use the original time evolution equation for non-boundary points
                    newC[n_x][n_y][n_z] =  (1/6) * (C[n_x+1][n_y][n_z] + C[n_x-1][n_y][n_z] + C[n_x][n_y+1][n_z] + C[n_x][n_y-1][n_z] + C[n_x][n_y][n_z+1] + C[n_x][n_y][n_z-1])
                
                C = newC[:][:][:]
"""  
    #Central Brightness and Error
    if n_t%10 ==0:
        B = Brightness(C,ds)[50,50]
        central_B.append(B)
        theoretical_B = theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50]
        error = abs((B-theoretical_B))
        B_error.append(error)

        #Calculating Radius
        # Calculating Radius
        radius = None
        for i in range(max_dist):
            if np.round(C[i,0,0],2) ==0.0:
                radius = (i-1)*ds
                radius_list.append(radius)
                error = np.max(np.abs(((C[:,0,0]-C_theoretical[:,n_t]))))
                r_error = radius*error
                radius_error.append(r_error)
                break
   
    
fig, ax1 = plt.subplots()
t = np.linspace(0,max_t-1,int(max_t/10))
color = 'tab:blue'
ax1.set_ylim(2000,10000)
ax1.set_xlabel('Time/s')
ax1.set_ylabel('Radius/m', color=color)
ax1.plot(t, radius_list, color=color,ls='--',label='radius')
ax1.tick_params(axis='y', labelcolor=color)
ax1.errorbar(t, radius_list, yerr=radius_error, fmt='o',capsize=5, color=color)


ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylim(0,7000)
ax2.set_ylabel(r'Brightness /$m^{-2}$', color=color)  
ax2.plot(t, central_B, color=color,ls='--',label='brightness')
ax2.tick_params(axis='y', labelcolor=color)
ax2.errorbar(t, central_B, yerr=B_error, fmt='o',capsize = 5, color=color)

# Adjust the position of the legends
ax1.legend(loc='upper right', bbox_to_anchor=(1, 0.9), fontsize = 'x-small')
ax2.legend(loc='upper right', bbox_to_anchor=(1, 1.0), fontsize = 'x-small')
#plt.title('Central Brightness and Radius Evolution Without Wind')
plt.show()

   

   



    

    #Calculating Error
    
    

    #Radius and Error



   

plt.clf()
central_B.append(Brightness(C,ds)[max_dist][max_dist])
plt.plot(central_B)
plt.xlabel('Time/s')
plt.ylabel(r'Brightness $m^{-2}$')
plt.title('Brightness at the centre of the sphere')
plt.pause(0.01)

"""
   

"""

    if n_t == 12:
        axs[0,0].plot(space_steps,C[:][0][0],label='Finite Difference',color='blue',ls='--')
        axs[0,0].plot(space_steps,C_theoretical[:,n_t],ls='--',label='Theoretical Solution',color='red')
        axs[0,0].legend(loc='upper right',fontsize = 'x-small')
        axs[0,0].set_title(r'$t=12s$')
        axs[0,0].set_xlabel(r'Distance from Centre / m')
        axs[0,0].set_ylabel(r'Number Density / $m^{-3}$')
        error_list.append(np.max(np.abs(((C[:,0,0]-C_theoretical[:,n_t])))))
        arg = (np.argmax((np.abs(((C[:,0,0]-C_theoretical[:,n_t]))))))
        arg_list.append(arg)
        value_list.append(C[arg,0,0])

    elif n_t==27:
        axs[0,1].plot(space_steps,C[:][0][0],label='Finite Difference',color='blue',ls='--')
        axs[0,1].plot(space_steps,C_theoretical[:,n_t],ls='--',label='Theoretical Solution',color='red')
        axs[0,1].legend(loc='upper right',fontsize = 'x-small')
        axs[0,1].set_title(r'$t=27s$')
        axs[0,1].set_xlabel(r'Distance from Centre / m')
        axs[0,1].set_ylabel(r'Number Density / $m^{-3}$')
        
        error_list.append(np.max(np.abs(((C[:,0,0]-C_theoretical[:,n_t])))))
        arg = (np.argmax((np.abs(((C[:,0,0]-C_theoretical[:,n_t]))))))
        arg_list.append(arg)
        value_list.append(C[arg,0,0])

       
    elif n_t == 64:
        axs[1,0].plot(space_steps,C[:][0][0],label='Finite Difference',color='blue',ls='--')
        axs[1,0].plot(space_steps,C_theoretical[:,n_t],ls='--',label='Theoretical Solution',color='red')
        axs[1,0].legend(loc='upper right',fontsize = 'x-small')
        axs[1,0].set_title(r'$t=64s$')
        axs[1,0].set_xlabel(r'Distance from Centre / m')
        axs[1,0].set_ylabel(r'Number Density / $m^{-3}$')
        error_list.append(np.max(np.abs(((C[:,0,0]-C_theoretical[:,n_t])))))
        arg = (np.argmax((np.abs(((C[:,0,0]-C_theoretical[:,n_t]))))))
        arg_list.append(arg)
        value_list.append(C[arg,0,0])

       
    elif n_t == 200:
        axs[1,1].plot(space_steps,C[:][0][0],label='Finite Difference',color='blue',ls='--')
        axs[1,1].plot(space_steps,C_theoretical[:,n_t],ls = '--',label='Theoretical Solution',color='red')
        axs[1,1].legend(loc = 'upper right',fontsize = 'x-small')
        axs[1,1].set_title(r'$t=200s$')
        axs[1,1].set_xlabel(r'Distance from Centre / m')
        axs[1,1].set_ylabel(r'Number Density / $m^{-3}$')
        error_list.append(np.max(np.abs(((C[:,0,0]-C_theoretical[:,n_t])))))
        arg = (np.argmax((np.abs(((C[:,0,0]-C_theoretical[:,n_t]))))))
        value_list.append(C[arg,0,0])
        arg_list.append(arg)
    plt.tight_layout()
    plt.show()
"""




#plotting the brightness
    if n_t ==12:
        img=axs[0,0].imshow(Brightness(C,ds))        
        axs[0,0].set_title(r'$t=12s$')
        axs[0,0].set_xticklabels((-6,0,6))
        axs[0,0].invert_yaxis()  # This line is added to flip the y-axis
        axs[0,0].set_yticklabels((-6,0,6))
        axs[0,0].set_xlabel('X/km')
        axs[0,0].set_ylabel('Y/km')
        colorbar = fig.colorbar(img, ax=axs[0,0])
        

        axs[2,0].plot(space_steps,(Brightness(C,ds)[50,50:]-theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50:])/(theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50]),ls='--',label='Finite Difference',color='blue')
        axs[1,0].plot(space_steps,theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50:],ls='--',label='Theoretical Solution',color='red')
        axs[1,0].plot(space_steps,Brightness(C,ds)[50,50:],ls='--',label='Finite Difference',color='blue')
    
        axs[1,0].legend(loc='upper right',fontsize = 'xx-small')
        axs[1,0].set_xlabel(r'Distance from Centre / m')
        axs[1,0].set_ylabel(r'Brightness $m^{-2}$')
        
        axs[2,0].set_xlabel(r'Distance from Centre / m')
        axs[2,0].set_ylabel(r'Relative Error')
        

    elif n_t ==27:
        axs[0,1].imshow(Brightness(C,ds))
        colorbar = fig.colorbar(axs[0,1].imshow(Brightness(C,ds)))
        axs[0,1].set_title(r'$t=27s$')
        axs[0,1].set_xticklabels((-6,0,6))
        axs[0,1].invert_yaxis()  # This line is added to flip the y-axis
        axs[0,1].set_yticklabels((-6,0,6))
    

    
        axs[2,1].plot(space_steps,(Brightness(C,ds)[50,50:]-theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50:])/theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50],ls='--',label='Finite Difference',color='blue')
        axs[1,1].plot(space_steps,theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50:],ls='--',label='Theoretical Solution',color='red')
        axs[1,1].plot(space_steps,Brightness(C,ds)[50,50:],ls='--',label='Finite Difference',color='blue')
        
        axs[1,1].legend(loc='upper right',fontsize = 'xx-small')
        
        axs[1,1].set_xlabel(r'Distance from Centre / m')
    
        axs[2,1].set_xlabel(r'Distance from Centre / m')
            


        elif n_t == 64:
            axs[0,2].imshow(Brightness(C,ds))
            colorbar = fig.colorbar(axs[0,2].imshow(Brightness(C,ds)))
            axs[0,2].set_title(r'$t=64s$')
            axs[2,2].plot(space_steps,(Brightness(C,ds)[50,50:]-theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50:])/theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50],ls='--',label='Finite Difference',color='blue')
            axs[1,2].plot(space_steps,theoretical_Brightness(ds,dt,initial_radius*ds,max_dist,n_t)[50,50:],ls='--',label='Theoretical Solution',color='red')
            axs[1,2].plot(space_steps,Brightness(C,ds)[50,50:],ls='--',label='Finite Difference',color='blue')
            axs[0,2].set_xticklabels((-6,0,6))
            axs[0,2].invert_yaxis()  # This line is added to flip the y-axis
            axs[0,2].set_yticklabels((-6,0,6))
            
            axs[1,2].legend(loc='upper right',fontsize = 'xx-small')
            axs[1,2].set_xlabel(r'Distance from Centre / m')
        
            axs[2,2].set_xlabel(r'Distance from Centre / m')

            
        #elif n_t == 200:
            #axs[2,0].imshow(img3)
            #axs[2,1].imshow(Brightness(C,ds))
            #axs[2,2].imshow(theoretical_Brightness(ds,dt,initial_radius,max_dist,n_t))
    #plt.subplots_adjust(wspace=0.5,hspace=0.5)
   # plt.show()
    
    




    
    #plot density distribution
    




        
        
#plot brightness



    
    

    











