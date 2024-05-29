#**************************************************************************************************#
#**************************************************************************************************#

initialize_3Dgrid <- function (grid_size_x, 
                               grid_size_y, 
                               array_names, 
                               homogeneous_land = FALSE, 
                               n_patches, 
                               homogenous_grid_quality){
  x_names <- as.character(seq(1, grid_size_x, 1))
  y_names <- as.character(seq(1, grid_size_y, 1))
  if (homogeneous_land) {
    land_quality <- rep(homogenous_grid_quality, grid_size_x * 
                          grid_size_y)
  }
  else {
    land_quality <- green_land(n_patches, grid_size_x, grid_size_y)
  }
  threeD_grid <- array(c(land_quality, rep(0, grid_size_x * grid_size_y)), 
                       dim = c(grid_size_x, grid_size_y, 2), 
                       dimnames = list(x_names, y_names, array_names))
  threeD_grid[, , 1][threeD_grid[, , 1] < 0] <- 0
  threeD_grid[, , 1][threeD_grid[, , 1] > 1] <- 1
  return(threeD_grid)
}

#**************************************************************************************************#
#**************************************************************************************************#
#function to get the horizontal and vertical lines to plot for distinguishing the cells in the grid
polygon_coordinates <- function(x,y){
  x_coordinates <- c(x-0.5, x+0.5, x+0.5, x-0.5)
  y_coordinates <- c(y-0.5, y-0.5, y+0.5, y+0.5)
  return(list(x_coordinates, y_coordinates))
}


#**************************************************************************************************#
#**************************************************************************************************#
initialize_agents <- function(N, 
                               state_names, 
                               state_initial_probability,
                               threeD_grid, 
                               max_forage_thresh = 1){
  agents <- data.frame(agent_no = 1, 
                       pos_x = sample(1:dim(threeD_grid)[1], size = 1),
                       pos_y = sample(1:dim(threeD_grid)[2], size = 1), 
                       state = sample(state_names, size = 1, prob = state_initial_probability)
                       )
  agents$land_QI[1] <- threeD_grid[agents$pos_x[1], agents$pos_y[1], 1]
  
  if(agents$state[1] == "foraging"){
    agents$energy_level[1] <- stats::runif(1, min = 0, max = max_forage_thresh)
  }else{
    agents$energy_level[1] <- 1
  }
  for(i in 2:N) {
    agents_more <- data.frame(agent_no = i, 
                              pos_x = sample(1:dim(threeD_grid)[1], size = 1), 
                              pos_y = sample(1:dim(threeD_grid)[2], size = 1), 
                              state = sample(state_names, size = 1, prob = state_initial_probability))
    agents_more$land_QI[1] <- threeD_grid[agents_more$pos_x[1], 
                                          agents_more$pos_y[1], 1]
    
    if (agents_more$state[1] == "foraging"){
      agents_more$energy_level[1] <- stats::runif(1, min = 0, max = max_forage_thresh)
    }else{
      agents_more$energy_level[1] <- 1
    }
    agents <- rbind(agents, agents_more)
  }
  return(agents)
}

#**************************************************************************************************#
#**************************************************************************************************#
add_connections_to_agents <- function(agents, 
                                      threeD_grid, 
                                      x_partitions, 
                                      y_partitions, 
                                      probability_connection){
  
  x_breaks <- seq(0, dim(threeD_grid)[1], by = dim(threeD_grid)[1]/x_partitions)
  y_breaks <- seq(0, dim(threeD_grid)[2], by = dim(threeD_grid)[2]/y_partitions)
  agents_group_wise <- NULL
  
  k <- 1
  n_groups_vector <- seq(1, (length(x_breaks) - 1) * (length(y_breaks) - 1), by = 1)
  
  for (i in 1:(length(x_breaks) - 1)) {
    agents_filtered <- agents[agents$pos_x > x_breaks[i] & 
                                agents$pos_x <= x_breaks[i + 1], ]
    for (j in 1:(length(y_breaks) - 1)) {
      group <- agents_filtered[agents_filtered$pos_y > 
                                 y_breaks[j] & agents_filtered$pos_y <= y_breaks[j + 
                                                                                   1], ]
      pre_connection_network <- igraph::erdos.renyi.game(nrow(group), 
                                                         probability_connection)
      igraph::V(pre_connection_network)$name <- as.character(group$agent_no)
      pre_connection_network <- igraph::set_vertex_attr(pre_connection_network, 
                                                        name = "group_number", value = n_groups_vector[k])
      pre_connection_network_matrix <- as.matrix(igraph::as_adjacency_matrix(pre_connection_network))
      group$connections <- sapply(group$agent_no, function(i) {
        names(which(pre_connection_network_matrix[as.character(i), 
        ] == 1))
      })
      group$group_number <- n_groups_vector[k]
      agents_group_wise <- rbind(agents_group_wise, group)
      if (k != 1) {
        complete_network <- igraph::disjoint_union(complete_network, 
                                                   pre_connection_network)
      }
      else {
        complete_network <- pre_connection_network
      }
      k <- k + 1
    }
  }
  agents_group_wise <- agents_group_wise[order(agents_group_wise$agent_no), 
  ]
  return(list(agents_group_wise, complete_network))
}


#**************************************************************************************************#
#**************************************************************************************************#
add_animals_to_gridcells <- function(agents, threeD_grid) 
{for (i in 1:nrow(agents)) {
    threeD_grid[agents$pos_x[i], agents$pos_y[i], 2] <- threeD_grid[agents$pos_x[i], 
                                                                    agents$pos_y[i], 2] + 1}
  return(threeD_grid)
}

#**************************************************************************************************#
#**************************************************************************************************#
#Function to plot the grid
grid_plot_function <- function(threeD_grid, 
                               time_step, 
                               plot_agents = TRUE, 
                               agents, 
                               n_groups = n_groups){
  
  plot(1, 2, type = "n", xlim = c(1, dim(threeD_grid)[1]), 
       ylim = c(1, dim(threeD_grid)[2]), col = 1:10, pch = 16, 
       xlab = "x coordinate", ylab = "y coordinate", main = paste(
         "Agents at time = ", time_step, sep = ""), xaxt = "n", 
       yaxt = "n")
  graphics::axis(side = 1, at = seq(1, dim(threeD_grid)[1], 
                                    1), labels = TRUE, tick = TRUE)
  graphics::axis(side = 2, at = seq(1, dim(threeD_grid)[2], 
                                    1), labels = TRUE, tick = TRUE)
  grid_col = grDevices::terrain.colors(100, alpha = 0.8, rev = TRUE)
  for (i in 1:nrow(threeD_grid[, , 1])) {
    for (j in 1:ncol(threeD_grid[, , 1])) {
      pc <- polygon_coordinates(i, j)
      graphics::polygon(pc[[1]], pc[[2]], 
                        col = grid_col[round(threeD_grid[i,j, 1] * 100, 0)], 
                        border = NA)
    }
  }
  set.seed(1)
  col_vec <- randomcoloR::distinctColorPalette(n_groups)
  if (plot_agents) {
    
    graphics::points(jitter(agents$pos_x, 0.5), 
                     jitter(agents$pos_y, 0.5), 
                     pch = 21, 
                     cex = 0.7, 
                     col = col_vec[as.numeric(as.factor(agents$group_number))],
                     bg = c("green", "red")[as.numeric(as.factor(agents$state))],
                     lwd = 2)
  }
}

#**************************************************************************************************#
#**************************************************************************************************#
animal_movements_record_list_initialize <- function(agents, total_time) 
{
  animal_movement_list <- list()
  animal_movement_list <- as.list(agents$agent_no)
  animal_movement_list <- lapply(animal_movement_list, function(i) {
    data.frame(animal_id = character(total_time), 
               time_point = numeric(total_time), 
               x_coordinate = numeric(total_time), 
               y_coordinate = numeric(total_time), 
               state = character(total_time), 
               energy_level = numeric(total_time))
  })
  names(animal_movement_list) <- as.character(agents$agent_no)
  for (i in 1:nrow(agents)) {
    animal_movement_list[[i]][1, ] <- c(paste0(agents$agent_no[i]), 
                                        1, agents$pos_x[i], agents$pos_y[i], agents$state[i], 
                                        agents$energy_level[i])
  }
  return(animal_movement_list)
}

#**************************************************************************************************#
#**************************************************************************************************#
depletion_function <- function(x){
  exponential_constant_times_x <- exp(3*x)
  return(exponential_constant_times_x/(exp(3)*2))
}

#**************************************************************************************************#
#**************************************************************************************************#
#Function that returns next choice based on food need
food_direction_probabilistic <- function(i, agents, threeD_grid, food_span){
  
  current_pos_x <- agents$pos_x[i]
  current_pos_y <- agents$pos_y[i]
  current_state <- agents$state[i]
  
  #changes in position
  x_range <- seq(current_pos_x - food_span, current_pos_x + food_span, 1)
  y_range <- seq(current_pos_y - food_span, current_pos_y + food_span, 1)
  x_range <- x_range[x_range>0 & x_range <= grid_size_x]
  y_range <- y_range[y_range>0 & y_range <= grid_size_y]
  
  possible_cells <- expand.grid(x_range, y_range)
  
  #Take the exponential of each cell
  exponential_food <- apply(possible_cells, 1, function(j) {
    exp(threeD_grid[j[1], j[2], 1])})
  exponential_food_probablity_vector <- exponential_food/sum(exponential_food)
  
  #Choose the cell based on above probability
  next_move_cell <- sample(nrow(possible_cells), size = 1, prob = exponential_food_probablity_vector)
  return(c(possible_cells[next_move_cell, 1], possible_cells[next_move_cell, 2]))
  
}

#**************************************************************************************************#
#**************************************************************************************************#
choose_random_cell <- function(current_pos_x, current_pos_y, grid, movement_span = 1){
  x_range <- seq(current_pos_x - movement_span, current_pos_x + movement_span, 1)
  y_range <- seq(current_pos_y - movement_span, current_pos_y + movement_span, 1)
  x_range <- x_range[x_range>0 & x_range <= grid_size_x]
  y_range <- y_range[y_range>0 & y_range <= grid_size_y]
  
  possible_cells <- expand.grid(x_range, y_range)
  random_row <- sample(1:nrow(possible_cells), size = 1)
  
  x_cell <- possible_cells$Var1[random_row]
  y_cell <- possible_cells$Var2[random_row]
  
  return(c(x_cell, y_cell))
}

#**************************************************************************************************#
#**************************************************************************************************#
#Function returns the coordinates of the nearest connections in the neighborhood so an agent can make a move towards them
social_direction <- function(i, agents, grid, search_span, step_length){
  
  #get agent's all connections from the population
  friends <- agents$connections[i][[1]]
  
  #get agent's current position
  current_pos_x <- agents$pos_x[i]
  current_pos_y <- agents$pos_y[i]
  
  #if no connections, choose random cell to make a move
  if(length(friends) == 0){
    next_cell <- choose_random_cell(current_pos_x, current_pos_y, grid, step_length)
    return(next_cell)
    
  }else{
    
    #look into the neighbourhood
    x_range <- seq(current_pos_x - search_span, current_pos_x + search_span, 1)
    y_range <- seq(current_pos_y - search_span, current_pos_y + search_span, 1)
    x_range <- x_range[x_range>0 & x_range <= grid_size_x]
    y_range <- y_range[y_range>0 & y_range <= grid_size_y]
    possible_cells <- expand.grid(x_range, y_range)
    
    #remove the current position of agent from possible values, so that if a friend is already present, it should move
    remove_ind <- which(possible_cells$Var1 == current_pos_x & possible_cells$Var2 == current_pos_y)
    possible_cells <- possible_cells[-remove_ind,]
    
    #returns the name of animals in neighboring cells
    neighbours_all <- c()
    for(k in 1:nrow(possible_cells)){
      neighbours <- which(apply(agents[, c("pos_x", "pos_y")], 1, function(x) all(x==c(possible_cells$Var1[k], possible_cells$Var2[k]))))
      neighbours_all <- c(neighbours_all, neighbours)
    }
    
    #Find out if any friend among the neighbors
    nearby_friends <- as.integer(friends[which(friends %in% neighbours_all)])
    
    #if friend/s present in neighborhood
    if(length(nearby_friends) != 0 ){
      
      #get euclidean distance between nearby friends
      distances <- sapply(nearby_friends, function(j){distance_compute(agents, i, j)})
      
      #Make a move towards your nearest friend
      move_towards <- nearby_friends[which.min(distances)]
      
      #Make a move in the direction of friend
      #get position of friend
      friend_pos_x <- agents$pos_x[agents$agent_no == move_towards]
      friend_pos_y <- agents$pos_y[agents$agent_no == move_towards]
      
      #check how far is the friend in terms of number of cells
      diff_x <- friend_pos_x - current_pos_x
      diff_y <- friend_pos_y - current_pos_y
      
      if(diff_x != 0){
        current_pos_x <- current_pos_x + (sign(diff_x)*step_length)
      }
      
      if(diff_y != 0){
        current_pos_y <- current_pos_y + (sign(diff_y)*step_length)
      }
      
      return(c(current_pos_x, current_pos_y))
      
    }else{#if no friend present in neighborhood
      next_cell <- choose_random_cell(current_pos_x, current_pos_y, grid, step_length)
      return(next_cell)
    }
  }
}


#**************************************************************************************************#
#**************************************************************************************************#
agent_values_manipulate_food_and_social <- function(i, agents, threeD_grid, social_thresh){
  
  if(agents$state[i] == "foraging" & agents$energy_level[i] < min_resting_thresh){
    
    if(agents$energy_level[i] < social_thresh & agents$energy_level[i] < agents$land_QI[i]){
      energy_level <- agents$energy_level[i] + depletion_function(threeD_grid[agents$pos_x[i], agents$pos_y[i], 1])
      next_steps <- food_direction_probabilistic(i, agents, threeD_grid, food_span = food_span)
    }else if(agents$energy_level[i] < social_thresh & agents$energy_level[i] >= agents$land_QI[i]){
      
      energy_level <- agents$energy_level[i] - energy_deplete_moving # Decrease in energy levels due to exploring and no eating
      next_steps <- food_direction_probabilistic(i, agents, threeD_grid, food_span = food_span)
    }else{
      energy_level <- agents$energy_level[i] - energy_deplete_moving # Decrease in energy levels due to exploring and no eating
      next_steps <- social_direction(i, agents, threeD_grid, search_span = search_span, step_length = 1)}
  }else{#resting phase
    # In resting phase, the position of animal remains same, just energy depletes
    energy_level <- agents$energy_level[i] - energy_deplete_resting
    next_steps <- c(agents$pos_x[i], agents$pos_y[i])
  }
  energy_level <- min(energy_level, 1)
  energy_level <- max(energy_level, 0)
  return(c(energy_level, next_steps[1], next_steps[2]))
}

#**************************************************************************************************#
#**************************************************************************************************#
new_coordinates_with_food_and_social <- function(agents_grouped, agents_filtered_grouped, mat_dist, i, distance_thresh){
  
  #select the name/ID of agent
  col_name <- as.character(agents_filtered_grouped$agent_no[i])
  
  #Take the mean of all agents that are in neighbourhood of the current position of selected agent
  new_x <- round(mean(agents_grouped$pos_x[agents_grouped$agent_no %in% names(which(mat_dist[,col_name] < distance_thresh))]))
  new_y <- round(mean(agents_grouped$pos_y[agents_grouped$agent_no %in% names(which(mat_dist[,col_name] < distance_thresh))]))
  
  return(c(new_x, new_y))
}


#**************************************************************************************************#
#**************************************************************************************************#
#Decision is changed only of those agents that are foraging and following friends. Although The decision is influenced by all of their friends irrespective of their states/energy level
final_movement_with_food_and_social <- function(agents_decision_step1, agents, distance_thresh){
  
  l <- mclapply(unique(agents_decision_step1$group_number),
                function(g){
                  agents_filtered_grouped <- agents_decision_step1[agents_decision_step1$group_number == g,]
                  agents_grouped <- agents[agents$group_number == g,]
                  mat <- as.matrix(cbind(agents_grouped$pos_x, agents_grouped$pos_y))
                  rownames(mat) <- agents_grouped$agent_no
                  distances <- dist(mat, method = "euclidean")
                  mat_dist <-  as.matrix(distances, labels = TRUE)
                  
                  for(ag in 1:nrow(agents_filtered_grouped)){
                    agents_filtered_grouped[ag, c("x_movement", "y_movement")] <- new_coordinates_with_food_and_social(agents_grouped, agents_filtered_grouped, mat_dist, ag, distance_thresh = distance_thresh)
                  }
                  return(agents_filtered_grouped)
                }, mc.cores = n_cores)
  
  agents_decision_step1 <- list.rbind(l)
  agents_decision_step1 <- agents_decision_step1[order(agents_decision_step1$agent_no),]
  return(agents_decision_step1)
}




