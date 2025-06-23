import numpy as np

def random_particles(box_length, box_width, box_height ,mean_radius, std_radius, fill_factor):
    # Calculate the number of particles based on the fill factor
    vol_box = box_length * box_width * box_height
    vol_particle = (4/3)*np.pi * (mean_radius)**3
    num_particles = int(vol_box * fill_factor / vol_particle)
    
    # Initialize list to store particle positions
    positions = []
    rad = []
    count = 0

    # Generate random positions for particles
    for _ in range(num_particles):
        while True:
            # Generate random coordinates for particle
            x = np.random.uniform(0, box_length)
            y = np.random.uniform(0, box_width)
            z = np.random.uniform(0, box_height)
            
            #Generate random radius for particle
            particle_radius=np.random.normal(loc=mean_radius,scale=std_radius)
            
            # Check if particle radius is valid
            if particle_radius < 0:
                break
            
            # Check if particle overlaps with existing particles
            overlap = False
            for pos in positions:
                if np.linalg.norm(np.array(pos) - np.array([x, y, z])) < 2*particle_radius + 2*particle_radius:
                    overlap = True
                    break
            
            # If no overlap, add particle position to list
            if not overlap:
                positions.append([x, y, z])
                rad.append(particle_radius)
                count = count + 1
                print(count)
                break

    return np.array(positions), np.array(rad)

#%%
# Define box dimensions and particle radius
box_length = 2e-6  
box_width = 2e-6
box_height = 2e-6
mean_radius = 200e-9
std_radius = 50e-9

# Define fill factor
fill_factor = 0.8
rand_pos_list, rand_rad_list = random_particles(box_length, box_width, box_height, mean_radius, std_radius,fill_factor)

rand_pos_list = rand_pos_list * 1e6  # Convert position from meters to micrometers
rand_rad_list = rand_rad_list * 1e6

# Export data to a CSV file
np.savetxt("mtlcad_3D_spheres_FF75_final.csv", np.column_stack((rand_pos_list, rand_rad_list)), delimiter=",", header="x,y,z,radius", comments='')