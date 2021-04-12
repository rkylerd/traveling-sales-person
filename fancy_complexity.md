### Overall
```
# initialize variables - time & space O(n); some variables are lists of length n

# get initial solution from greedy algorithm - Time O(n^3), Space O(n) to store the best route 

# perform ant colony min-max algorithm - Time O(n^3)

```

### Greedy Complexity
``` 
# Time Complexity O(n^3) or O(r * n^3)
for i in range(n): # check every possible greedy solution (each starting city) - O(n) time to iterate through n
    
    # initialize the cost (0) and route ([i])
    for _ in range(n - 1): # Build new greedy solution from i - O(n) time to iterate through n-1 cities
        
        for j in range(n): # Find the next closest city and append to route - O(n) time & space
            # 'if not j in route' could add another n or constant r to the overall complexity
    
        # add the city j with shortest distance
        # set i equal to the city with the shortest distance
    
    # update best greedy solution if a better one was found 
```

### Ant Colony Min-Max Algorithm
```
while the timer has not ended and iterations_count < iteration_limit: # time O(1)
    
    for ant k in ants: # iterate through all ants - O(k)
        # start from random city i
        # add i to the current route

        for _ in (n - 1): # Build the route city by city - O(n) time

            for j in n: # Calculate & store the probability of traveling to each city - Time & Space O(n) (possibly time O(n^2) because of the use of 'j in route')
                # the probability of visiting city j from i is based on a combination of the pheromone level and the distance between

            # randomly use one of the following methods to determine the next city to visit:
                - select the city with the highest probability - time O(1)
                - consider the probablities of all possible cities (calculated above) - time O(n)
            
            # add the selected city to the current route - time O(1)
            i = the selected city
        
        # Update iteration_best_route if a better route was found - time O(1)

    # if the best_iteration_ant was better than the global_best_solution, update the global_best_solution - time O(1)

    # Reduce the pheromone levels along all edges (but not below the minimum level) - time O(n^2)

    # Randomly choose to increase the pheromone levels of each edge in the iteration_best_route or the global_best_solution - O(n) time to increase pheromone levels
```

### Create TSP solution from best route of cities
```
route = []
for i in BSSF_route: # O(n) time & space
    route.append(cities[i])

solution = TSPSolution(route)
```