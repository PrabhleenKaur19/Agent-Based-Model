#Main simulations code
set.seed(seed_grid)
threeD_grid <- initialize_3Dgrid(grid_size_x = grid_size_x, grid_size_y = grid_size_y, 
                                 array_names <- c("land_quality", "n_animals"),
                                 homogeneous_land = FALSE,
                                 n_patches = round((grid_size_x*grid_size_y)/120))
initial_grid_quality <- threeD_grid[,,1]

grid_plot_function(threeD_grid, 1, plot_agents = FALSE)

set.seed(seed_grid)
agents_initial <- initialize_agents(N = N, 
                                    state_names = state_names, 
                                    state_initial_probability = state_initial_probability, 
                                    threeD_grid = threeD_grid,
                                    max_forage_thresh = max_forage_thresh)
agents_generated <- add_connections_to_agents(agents_initial, 
                                                      threeD_grid, 
                                                      x_partitions = 4, 
                                                      y_partitions = 4,
                                                      probability_connection = probability_edge)

agents <- agents_generated[[1]]
agents_pre_network <- agents_generated[[2]]
col_vec <- randomcoloR::distinctColorPalette(n_groups)

#Example visualisation of a single group
g=subgraph(graph=agents_pre_network, vids=which(V(agents_pre_network)$group_number==4))
plot(g)
aniSNA::get_network_summary(agents_pre_network)

matrix_agents_pre_network <- igraph::as_adjacency_matrix(agents_pre_network)
assort_list <- assortnet::assortment.discrete(matrix_agents_pre_network, types = V(agents_pre_network)$group_number
                                              , weighted = FALSE, SE = FALSE, M = 1)
assort_list$r #initial assortativity is 1

#Update grid
threeD_grid <- add_animals_to_gridcells(agents = agents, threeD_grid = threeD_grid)

#Plot grid along with animals
grid_plot_function(threeD_grid, 1, plot_agents = TRUE, agents = agents, n_groups = n_groups)

animal_movement_list <- animal_movements_record_list_initialize(agents, total_time)

#start with an empty network
network <- matrix(0, nrow=N, ncol=N)

for (t in 2:total_time) {
  
  #the whole area is replenished at a certain rate at each time step
  threeD_grid[,,1] <- threeD_grid[,,1] + (replenish_constant*threeD_grid[,,1])
  
  #Ensure that the resource index doesnt exceed initially assigned values
  threeD_grid[,,1][threeD_grid[,,1] >initial_grid_quality] <- initial_grid_quality[threeD_grid[,,1] > initial_grid_quality]
  
  #filter out the coordinates from threeDgrid that are changed due to foraging
  coordinates_foraging_agents <- as.matrix(agents[,c("pos_x", "pos_y")][agents$energy_level < agents$land_QI & 
                                                                          agents$energy_level < social_thresh & 
                                                                          agents$state == "foraging",])
  
  #subset agents whose state == "foraging" and energy_level > social threshold
  agent_two_step_index <- which(agents$energy_level >= social_thresh & agents$state == "foraging")
  agents_decision_step1 <- agents[agent_two_step_index, c("agent_no", "group_number")]
  
  #Change in agent properties due to movement/resting
  agents[,c("energy_level", "pos_x", "pos_y")] <- list.rbind(mclapply(1:N, function(i) 
    agent_values_manipulate_food_and_social(i, agents, threeD_grid, social_thresh), mc.cores = n_cores))
  
  #agents final movements will change again based on collective movements aggregated in step 2
  agents_decision_step1[,c("x_coordinate", "y_coordinate")] <- agents[agent_two_step_index, c("pos_x", "pos_y")]
  agents_decision_step1[,c("x_movement", "y_movement")] <- agents[agent_two_step_index, c("pos_x", "pos_y")]
  agents_decision_final_movement <- final_movement_with_food_and_social(agents_decision_step1, agents, distance_thresh = 1.5)
  agents[agent_two_step_index,c("pos_x", "pos_y")] <- agents_decision_final_movement[, c("x_movement", "y_movement")]
  
  #If agents foraging and decide to eat, decrease in threeDgrid values
  threeD_grid_food <- threeD_grid[,,1]
  threeD_grid_food[coordinates_foraging_agents] <- threeD_grid_food[coordinates_foraging_agents] - depletion_function(threeD_grid_food[coordinates_foraging_agents])
  threeD_grid[,,1] <- threeD_grid_food
  
  agents[,c("land_QI", "state")] <- list.rbind(mclapply(1:N, function(i){
    land_QI <- threeD_grid[agents$pos_x[i], agents$pos_y[i], 1]
    state <- new_state(i, t, agents$state[i], agents$energy_level[i], animal_movement_list, min_resting_thresh)
    return(c(land_QI, state))}, mc.cores = n_cores))
  
  animal_movement_list <- mclapply(seq_along(animal_movement_list), function(i) 
    {animal_movement_list[[i]][t,] <- c(paste0(agents$agent_no[i]), 
                                        t, 
                                        agents$pos_x[i], 
                                        agents$pos_y[i], 
                                        agents$state[i], 
                                        agents$energy_level[i])
  return(animal_movement_list[[i]])}, mc.cores = n_cores)
  
  # restrict min and max values for grid quality
  threeD_grid[, , 1][threeD_grid[, , 1] < 0] <- 0
  threeD_grid[, , 1][threeD_grid[, , 1] > 1] <- 1
  
  # Plot the grid and the animals at each time point
  # grid_plot_function(threeD_grid, t, plot_agents = TRUE, agents = agents, n_groups = n_groups)
  
  #Update network values
  coords <- data.frame(X=agents$pos_x, Y=agents$pos_y)
  dists <- as.matrix(dist(coords))
  network[which(dists<1,arr.ind=T)] <- network[which(dists<1,arr.ind=T)] + 1
  
  #Check assortativity values every 50 timepoints
  if(t%%50 == 0){
    network_after_tsteps <- network/t
    assort_list_network <- assortnet::assortment.discrete(network_after_tsteps, 
                                                          types = agents$group_number, 
                                                          weighted = TRUE, 
                                                          SE = FALSE, 
                                                          M = 1)
    cat("Assortivity value at time point",t, "-----",assort_list_network$r, "\n")
  }
}
