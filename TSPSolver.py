#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
import random


class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: 1. cost of best solution, 
		2. time spent to find best solution, 3. total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	# stringToIndex() turns a city name of 'A' to index 0, for example
	# changes 'AA' to 26
	# This function is useful because the cities are ordered by city name
	def stringToIndex(self, s):
		returnIdx = 0
		# For all leading characters 
		# (in the string 'ABY', 'AB' are the leading characters)
		for idx in range(len(s)-1):
			# We don't want to multiply by 0, so use 64 NOT 65
			# 1 = 'A' - 64 
			pos = ord(s[idx]) - 64 
			returnIdx += 26 * pos

		return returnIdx + ord(s[-1])-65

	def getIndexOfClosestCity( self, all_cities, remaining_cities, city_idx):
		city = remaining_cities[0]
		closest_city_idx = self.stringToIndex(city._name)
		city_idx_to_pop = 0

		if len(remaining_cities) == 1:
			return closest_city_idx, city_idx_to_pop

		shortest_path = all_cities[city_idx].costTo(remaining_cities[0])

		for idx in range(1, len(remaining_cities)):
			
			city = remaining_cities[idx]
			if shortest_path > all_cities[city_idx].costTo(city):
				shortest_path = all_cities[city_idx].costTo(city)
				# index of the city in the complete list, NOT in the remaining list of cities
				closest_city_idx = self.stringToIndex(city._name)
				city_idx_to_pop = idx

		return closest_city_idx, city_idx_to_pop

	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		btsf_time = None
		start_cities = [False] * ncities

		while (time.time() - start_time) < time_allowance and count < ncities:
			new_timer = time.time()
			new_solution = None
			# start from a random city that hasn't been started from yet
			city_idx = random.randint(0, len(cities)-1 )
			while start_cities[city_idx]:
				city_idx = random.randint(0, len(cities) - 1)
			start_cities[city_idx] = True

			route = [ cities[city_idx] ]
			
			remaining_cities = [ *cities[:city_idx], *cities[city_idx+1:] ]
			
			while len(remaining_cities) > 0:
				# find the city with the smallest cost
				closest_city_idx, city_idx_to_pop = self.getIndexOfClosestCity( cities, remaining_cities, city_idx )
				
				# remove the smallest-cost city (so we don't return to it)
				remaining_cities.pop(city_idx_to_pop)
				route.append( cities[closest_city_idx] )
				city_idx = closest_city_idx

			new_solution = TSPSolution(route)
			count += 1
			if new_solution.cost < math.inf and (not(bssf) or bssf.cost > new_solution.cost):
				# Found a valid route
				btsf_time = time.time() - new_timer
				foundTour = True
				bssf = new_solution

		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = btsf_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results
	
	
	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		pass



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):

		cities = self._scenario.getCities()
		n = len(cities)
		a = 1
		b = 1
		p = 0.2
		ants = n
		max_iterations = 1000
		q = 0.5
		best_frequency = 0.25
		distances = [[cities[i].costTo(cities[j]) for j in range(n)] for i in range(n)]
		BSSF = math.inf
		BSSF_route = []

		# Set initial tour cost using the greedy algorithm
		for i in range(n): # check every possible greedy solution (each starting city)
			route = [i]
			cost = 0
			current = i
			for _ in range(n - 1): # Build path
				next = -1
				next_cost = math.inf
				for j in range(n):
					if not j in route and distances[current][j] < next_cost: # Get shortest path
						next = j
						next_cost = distances[current][j]
				if next == -1: # No path exists
					cost = math.inf
					break
				route.append(next)
				current = next
				cost += next_cost
			cost += distances[route[-1]][route[0]]
			if cost < math.inf and cost < BSSF:
				BSSF_route = route
				BSSF = cost

		t_max = 1 / (p * BSSF)
		t_min = t_max / (2 * n)
		pheromones = [[t_max] * n for _ in range(n)]

		count = 0

		# ant_cost = [0.0] * ants
		# ant_history = [[] for _ in range(ants)]
		iteration = 0
		prob = [0.0]*n
		start_time = time.time()
		last_improvement = 0 # iteration number of last improvement
		# Begin iterating
		while time.time() - start_time < time_allowance and iteration - last_improvement < max_iterations:
			iteration += 1
			iterationBestCost = math.inf
			iterationBestRoute = []
			# iterate through each ant
			for k in range(ants):
				# random starting city
				i = random.randrange(0, n)
				route = [i]
				ant_cost = 0

				# Build route city by city
				for _ in range(n - 1):
					sum = 0
					max_index = -1
					max_value = 0
					# Calculate the probability of traveling to each city
					for j in range(n):
						if j in route: # Can't revisit a city
							prob[j] = 0
						else:
							dist = distances[i][j]
							if dist == 0: # Avoid division by zero
								prob[j] = math.inf
								sum += 1
								max_value = math.inf
								max_index = j
							elif dist < math.inf:
								# Set the probability with a combination of the pheromone level and the distance
								prob[j] = (pheromones[i][j] ** a) * ((1/dist) ** b)
								sum += prob[j]
								if prob[j] > max_value:
									max_index = j
									max_value = prob[j]
							else: # Don't travel infinite distances
								prob[j] = 0

					if sum == 0: # No valid paths were found
						ant_cost = math.inf
						break

					if random.random() < q: # Randomly use the best option or choose probabilistically
						j = max_index
					else:
						choice = random.random()
						j = -1
						while choice >= 0:
							j += 1
							choice -= prob[j] / sum

					# Add to route
					route.append(j)
					ant_cost += distances[i][j]
					i = j

				# Finish route and update iteration best if it is a better route
				if ant_cost < math.inf:
					ant_cost += distances[route[-1]][route[0]]
				if ant_cost < iterationBestCost:
					iterationBestRoute = route
					iterationBestCost = ant_cost
			# if the best ant was better than the global best solution, update the global best solution
			if iterationBestCost < BSSF:
				BSSF = iterationBestCost
				BSSF_route = iterationBestRoute
				print("Iteration {}: {}".format(iteration, BSSF))
				last_improvement = iteration
				t_max = 1 / (p * BSSF)
				t_min = t_max / (2 * n)
				count += 1
			# Reduce all edges
			for i in range(n):
				for j in range(n):
					pheromones[i][j] *= (1 - p)
					if pheromones[i][j] < t_min:
						pheromones[i][j] = t_min
			# Add pheromones to each edge with either the iteration best or global best solution
			if random.random() < best_frequency:
				if BSSF < math.inf:
					for i in range(n):
						pheromones[BSSF_route[i]][BSSF_route[(i + 1) % n]] += 1 / BSSF
						if pheromones[BSSF_route[i]][BSSF_route[(i + 1) % n]] > t_max:
							pheromones[BSSF_route[i]][BSSF_route[(i + 1) % n]] = t_max
			elif iterationBestCost < math.inf:
				for i in range(n):
					pheromones[iterationBestRoute[i]][iterationBestRoute[(i + 1) % n]] += 1 / iterationBestCost
					if pheromones[iterationBestRoute[i]][iterationBestRoute[(i + 1) % n]] > t_max:
						pheromones[iterationBestRoute[i]][iterationBestRoute[(i + 1) % n]] = t_max
		end_time = time.time()

		route = []
		for i in BSSF_route:
			route.append(cities[i])
		solution = TSPSolution(route)

		results = {}
		results['cost'] = BSSF
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = solution
		results['max'] = iteration
		results['total'] = None
		results['pruned'] = None

		return results


