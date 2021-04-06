#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))



import copy
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

	def constructResults( self, cost, solutionTime, count, soln, maxQueueSize = None, total = None, pruned = None):
		results = {}
		results['cost'] = cost
		results['time'] = solutionTime if solutionTime != None else time.time() 
		results['count'] = count
		results['soln'] = soln
		results['max'] = maxQueueSize
		results['total'] = total
		results['pruned'] = pruned
		return results

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
		closest_city_idx = city._index
		city_idx_to_pop = 0

		if len(remaining_cities) == 1:
			return closest_city_idx, city_idx_to_pop

		shortest_path = all_cities[city_idx].costTo(remaining_cities[0])

		for idx in range(1, len(remaining_cities)):
			
			city = remaining_cities[idx]
			if shortest_path > all_cities[city_idx].costTo(city):
				shortest_path = all_cities[city_idx].costTo(city)
				# index of the city in the complete list, NOT in the remaining list of cities
				closest_city_idx = city._index
				city_idx_to_pop = idx

		return closest_city_idx, city_idx_to_pop

	def greedy( self,time_allowance=60.0, quitAfterFirstSolution=False ):
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
			solution_timer = time.time()
			new_solution = None
			# quiteAfterFirstSolution is to give the branchAndBound algorithm an initial bssf
			if bssf and quitAfterFirstSolution:
				break 
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
				bssf_time = time.time() - solution_timer
				foundTour = True
				bssf = new_solution

		return self.constructResults(
			bssf.cost if foundTour else math.inf, 
			bssf_time, 
			count, 
			bssf)
	
	
	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		cities = self._scenario.getCities()
		ncities = len(cities)
		bbStartTime = time.time()

		greedySolution = self.greedy(time_allowance, True)
		# Start over with new random seed
		if not(greedySolution['soln']):
			return self.constructResults(math.inf, None, 0, None)

		self.bpsf = [city._index for city in greedySolution['soln'].route]
		self.bssfCost = greedySolution['cost']
		self.bbBstsf = time_allowance
		self.bbSolutionCount = 0
		self.bbSolutionTimer = time.time()
		self.pruneCount = 0
		self.maxQueueSize = 0
		self.totalStates = 0
		
		# Get the initial state
		matrix = self.createCostMatrix(cities, ncities)
		initialCost = self.reduceMatrix(matrix, ncities)
		
		# Create a priority queue of the first level in the tree of cities
		Q, bounds, matrices = self.makeQueue(cities, copy.deepcopy(matrix), ncities, initialCost)
		# while first-level states still exist
		while len(Q):
			if (time.time() - bbStartTime) >= time_allowance:
				break
			self.totalStates += 1
			nextCity = Q.pop(0)
			lBoundOfNextCity = bounds.pop(0)

			# always check bound before trying a path (ordered by lowest bounds)
			if lBoundOfNextCity < self.bssfCost:
				indexOfCity = nextCity._index
				remainingCityIndeces = [ *[n for n in range(1,indexOfCity)], *[n for n in range(indexOfCity+1,ncities)]]
				
				# Always start from city 0
				pathSoFar = [0, indexOfCity]
				localMatrix = copy.deepcopy(matrices[indexOfCity-1])
				self.findSolution(localMatrix, cities, ncities, [*pathSoFar], [*remainingCityIndeces], lBoundOfNextCity, bbStartTime, time_allowance)
		route = []		
		for i in self.bpsf:
			route.append(cities[i])

		bssf = TSPSolution(route)

		return self.constructResults(
			self.bssfCost, 
			self.bbBstsf, 
			self.bbSolutionCount, 
			bssf, 
			self.maxQueueSize, 
			self.totalStates, 
			self.pruneCount)


	def markPathVisited(self, matrix, n, row, col):
		self.markDirectionVisited(matrix, n, row)
		self.markDirectionVisited(matrix, n, col, False)
		matrix[col][row] = None

	def findSolution(self, matrix, cities, n, pathSoFar, remainingCityIndeces, parentLowerBound, bbStartTime, timeAllowance):
		remainingCityIndecesInDepth = [*remainingCityIndeces]
		indexInOriginalRemainingCitiesIndeces = -1

		while len(remainingCityIndeces):
			if (time.time() - bbStartTime) >= timeAllowance:
				return
			self.totalStates += 1
			localMatrix = copy.deepcopy(matrix)
			indexInOriginalRemainingCitiesIndeces += 1
			index = remainingCityIndeces.pop(0)

			lb = self.getLowerBound(localMatrix, parentLowerBound, n, pathSoFar[-1], index)
			
			if lb < self.bssfCost:
				newPathSoFar = [*pathSoFar, index]

				if len(newPathSoFar) == n:
					self.bbBstsf = time.time() - self.bbSolutionTimer
					self.bbSolutionCount += 1
					self.bssfCost = lb
					self.bpsf = newPathSoFar
				else:
					self.findSolution(copy.deepcopy(localMatrix), cities, n, [*newPathSoFar], [ *remainingCityIndecesInDepth[:indexInOriginalRemainingCitiesIndeces], *remainingCityIndecesInDepth[indexInOriginalRemainingCitiesIndeces+1:]], lb, bbStartTime, timeAllowance)
			else:
				self.pruneCount += 1

	def makeQueue(self, cities, matrix, n, initialCost):
		bounds = []
		Q = []
		row = 0
		matrices = []

		# Start at column 1 because the initial state always begins at city 0
		for col in range(1,n):
			firstLevelMatrix = copy.deepcopy(matrix)
			lb = self.getLowerBound(firstLevelMatrix, initialCost, n, row, col)
			
			# If on the first city, add it to the queue and move on to the next city
			if len(bounds) == 0:
				bounds.append(lb)
				Q.append(cities[col])
				matrices.append(firstLevelMatrix)
				continue
			
			# To make the queue a priority queue, 
			# cities with smaller bounds are placed at the front of the queue
			indexOfNextSmallerThanMe = -1
			for idx, bound in enumerate(bounds):
				if bound > lb:
					break
				# the new lb needs will be placed after the current bound
				indexOfNextSmallerThanMe += 1 
			
			# Insert the new lower bound after the final element with a lower bound (the next smallest)
			i = 1 + indexOfNextSmallerThanMe
			bounds.insert(i, lb)
			Q.insert(i, cities[col])
			matrices.append(firstLevelMatrix)
		return Q, bounds, matrices

	def getLowerBound(self, matrix, parentLowerBound, n, row, col, ignoreRowCol=True, markNoneBeforeReduction=True):
		costAtCoordinate = matrix[row][col]
		if markNoneBeforeReduction:
			# set all values in the row and column to None (so we don't visit them again) 
			self.markPathVisited(matrix, n, row, col)
		optimizationCost = self.reduceMatrix(matrix, n, ignoreRowCol, row, col)
		return parentLowerBound + costAtCoordinate + optimizationCost

	def reduceMatrix(self, matrix, n, ignoreRowCol=False, ignoreRow=None, ignoreCol=None):
		costToReduce = 0

		for row in range(n):
			if ignoreRowCol and row == ignoreRow:
				continue

			rowSmallestValue = self.getSmallestElementInDirection(matrix, n, row, True, ignoreCol)

			# 0 < rowSmallestValue < infinity
			if rowSmallestValue and rowSmallestValue < math.inf:
				self.reduceElementsInDirection(matrix, n, rowSmallestValue, row)
			costToReduce += rowSmallestValue
		
		for col in range(n):
			if ignoreRowCol and col == ignoreCol:
				continue

			colSmallestValue = self.getSmallestElementInDirection(matrix, n, col, False, ignoreRow)
			# 0 < colSmallestValue < infinity
			if colSmallestValue and colSmallestValue < math.inf:
				self.reduceElementsInDirection(matrix, n, colSmallestValue, col, False)
			costToReduce += colSmallestValue
		
		return costToReduce

	def getSmallestElementInDirection(self, matrix, n, index, isRow, ignoreRowCol=None, defaultSmallestValue=0):

		smallestValue = math.inf
		flag = False
		if isRow:	
			for col in range(n):
				if ignoreRowCol != None and col == ignoreRowCol:
					continue
				if matrix[index][col] == None:
					continue
				flag = True
				if matrix[index][col] < smallestValue:
					smallestValue = matrix[index][col]
		else:
			for row in range(n):
				if ignoreRowCol != None and row == ignoreRowCol:
					continue
				if matrix[row][index] == None:
					continue
				flag = True
				if matrix[row][index] < smallestValue:
					smallestValue = matrix[row][index]

		return smallestValue if flag else defaultSmallestValue

	def reduceElementsInDirection(self, matrix, n, reducer, index, isRow = True):

		if isRow:	
			for col in range(n):
				if matrix[index][col] != None and matrix[index][col] < math.inf:
					matrix[index][col] -= reducer 
		else:
			for row in range(n):
				if matrix[row][index] != None and matrix[row][index] < math.inf:
					matrix[row][index] -= reducer

	def markDirectionVisited(self, matrix, n, index, isRow = True):
		
		if isRow:	
			for col in range(n):
				matrix[index][col] = None 
		else:
			for row in range(n):
				matrix[row][index] = None

	# Space complexity - O(n^2) for 2D matrix
	# Time complexity - O(n^2) to create the 2D matrix
	def createCostMatrix(self, cities, n):
		# O(n) time and space to setup rows
		matrix = [[] for _ in range(n)]

		# O(n^2) time and space
		for row in range(n):
			for col in range(n):
				# next two operations take constant time and space
				iToJCost = cities[row].costTo(cities[col]) if row != col else math.inf
				matrix[row].append( iToJCost )
		return matrix

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


