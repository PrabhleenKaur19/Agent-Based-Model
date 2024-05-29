#Set the parameter values in this file
seed_grid <- 12

#Specify number of rows and columns of land grid
grid_size_x <- 40
grid_size_y <- 40

#We initialise a 3 dimensional grid where the third dimension provides information about the 
#current resource quality index and number of animals present in each cell
array_names <- c("land_quality", "n_animals")

#Specify number of agents
N <- 500

#Total time steps
total_time = 4380

#Specify the probability vector for different states of "resting", "foraging"
state_names <- c("resting", "foraging")
state_initial_probability = c(0.3, 0.7)

movement_span <- 1
replenish_constant <- 0.25 *(N/(grid_size_x*grid_size_y)) #This value depends on the ratio of N and grid size
social_thresh <- 0.5
energy_deplete_moving <- 0.1
energy_deplete_resting <- 0.05

#state change parametres
max_forage_thresh <- 1
min_resting_thresh <- 0.75 #animal wont rest once it's energy declines to 0.75 or below

#If want to initialise a homogenous land with equvivalent total resource
#homogenous_grid_quality <- (sum(threeD_grid[,,1]))/(grid_size_x*grid_size_y)

search_span = 2
food_span = 1

#Agent Groups
x_partitions <- 4
y_partitions <- 4
n_groups <- x_partitions * y_partitions #Total 16 groups required
probability_edge <- 0.9

#number of cores for parallel processing
n_cores <- 8



