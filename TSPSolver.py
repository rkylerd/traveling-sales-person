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
		bssf_time = None
		start_cities = [False] * ncities
		starting_cities_remaining = ncities

		while (time.time() - start_time) < time_allowance:
			new_timer = time.time()
			new_solution = None
			if starting_cities_remaining == 0:
				break
			# start from a random city that hasn't been started from yet
			city_idx = random.randint(0, len(cities)-1 )
			while start_cities[city_idx]:
				city_idx = random.randint(0, len(cities) - 1)
			starting_cities_remaining -= 1
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
				bssf_time = time.time() - new_timer
				foundTour = True
				bssf = new_solution

		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = bssf_time
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
		# parameters (they usually have values between 0 and 1)
		a = 1
		b = 2
		p = 0.2
		q_not = .5
		xi = p
		t_not = .6

		ants = 25
		max_iterations = 1000
		pheromones = [[1] * n for _ in range(n)]
		
		ant_cost = [0.0]*ants
		ant_history = [[] for _ in range(ants)]
		distances = [[cities[i].costTo(cities[j]) for j in range(n)] for i in range(n)]
		iteration = 0
		prob = [0.0]*n
		BSSF = math.inf
		BSSF_route = []
		BSSF_time = None
		count = 0
		start_time = time.time()
		while time.time() - start_time < time_allowance and iteration < max_iterations:
			solution_start_time = time.time()
			iteration += 1
			bestAnt = -1
			bestCost = math.inf
			# for every ant k
			for k in range(ants):
				# random starting city i
				i = random.randrange(0, n)
				ant_history[k] = [i]
				ant_cost[k] = 0
				
				# visit every city ( '_' is just a counter to ensure i get to be every city n )
				for _ in range(n - 1):
					# sum of probabilities (used for randomly selecting the next city)
					sum = 0
					# look at each city to determin if reachable from i
					for j in range(n):
						# if ant k has already visited city j (or is at city j), don't visit j again
						if j in ant_history[k]:
							prob[j] = 0
						else:
							dist = distances[i][j]
							
							if dist == 0:
								prob[j] = math.inf
								sum += 1
							elif dist < math.inf:
								prob[j] = (pheromones[i][j] ** a) * ((1/dist) ** b)
								sum += prob[j]
							else:
								prob[j] = 0
					# try a new ant if no other cities are visitable from city i
					if sum == 0:
						ant_cost[k] = math.inf
						break
					
					# ACS: test whether or not to use tour construction (the most probable path) 
					# 	or normal probabilistic path selection
					if random.random() <=  q_not:
						most_probable_path = 0
						for j in range(1,n):
							if prob[j] > prob[most_probable_path]:
								most_probable_path = j
						j = most_probable_path
					else :
						# cities with higher probability are more likely to cause choice to go below 0 
						choice = random.random()
						j = -1
						while choice >= 0:
							j += 1
							choice -= prob[j] / sum

					ant_history[k].append(j)
					ant_cost[k] += distances[i][j]

					# ACS: evaporate pheromone level immediately after taking path i to j
					pheromones[i][j] = (1 - xi) * pheromones[i][j] + xi * t_not

					# set next starting city
					i = j

				if ant_cost[k] < math.inf:
					# Is the original city reachable from the last city?
					cost = distances[ant_history[k][n-1]][ant_history[k][0]]
					if cost < math.inf:
						ant_cost[k] += cost
					else:
						ant_cost[k] = math.inf
				if ant_cost[k] < bestCost:
					bestAnt = k
					bestCost = ant_cost[k]
			# if the best ant of all ant iterations found better than global BSSF 
			if bestAnt != -1:
				count += 1
				if bestCost < BSSF:
					BSSF_time = time.time() - solution_start_time
					BSSF = bestCost
					BSSF_route = ant_history[bestAnt]
					print("Iteration {}: {}".format(iteration, BSSF))
			
			# ACS: only the global best ant is allowed
			# 	to add pheromone after each iteration.
			for i in BSSF_route:
				current = pheromones[ BSSF_route[i] ][ BSSF_route[ (i+1) % len(BSSF_route) ] ] 
				pheromones[ BSSF_route[i] ][ BSSF_route[ (i+1) % len(BSSF_route) ] ] = (1 - p) * current + p * (1 / BSSF)
		end_time = time.time()

		route = []
		for i in BSSF_route:
			route.append(cities[i])
		solution = TSPSolution(route)

		results = {}
		results['cost'] = BSSF
		results['time'] = BSSF_time
		results['count'] = count
		results['soln'] = solution
		results['max'] = iteration
		results['total'] = None
		results['pruned'] = None

		return results


