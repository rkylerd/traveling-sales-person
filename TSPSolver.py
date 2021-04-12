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
	
	# O(r) represents the number of remaining cities (it's indicitave of the tree depth)
	def branchAndBound( self, time_allowance=60.0 ):
		cities = self._scenario.getCities()
		ncities = len(cities)
		self.bbStartTime = time.time()
		self.timeAllowance = time_allowance

		# Calculate an initial bssf using the greedy algorithm
		# O(n^2) time complexity, O(n) space complexity 
		greedySolution = self.greedy(time_allowance)

		# If no initial solution could be found using greedy, start over with a new random seed
		if not(greedySolution['soln']):
			return self.constructResults(math.inf, None, 0, None)

		# Initialize reporting stats; use the applicable greedy stats as defaults
		self.bpsf = [city._index for city in greedySolution['soln'].route]
		self.bssfCost = greedySolution['cost']
		self.bbBstsf = time_allowance
		self.bbSolutionCount = 0
		self.pruneCount = 0
		self.runningQueueSize = 0
		self.maxQueueSize = 0
		self.totalStates = 0

		# O(n^2) for both time and space complexity to crreate the n^2 cost matrix
		matrix = self.createCostMatrix(cities, ncities)
		
		# Get the initial state
		# O(n^2) for both time and space complexity to reduce
		initialCost = self.reduceMatrix(matrix, ncities)

		# Always start from city 0
		pathSoFar = [0]
		# O(r) space complexity (basically n at this point)
		remainingCityIndeces = [n for n in range(1,ncities)]

		# O(r * n^2) time and space complexity
		self.prioritizeByBounds(matrix, cities, ncities, [*pathSoFar], [*remainingCityIndeces], initialCost)

		# Create a TSP solution using the path made from branch and bound/greedy
		route = []	
		# O(n) time and space complexity 
		for i in self.bpsf:
			route.append(cities[i])
		bssf = TSPSolution(route)

		# if the time ran out before the algorithm found a better solution,
		# report the time_allowance as the time taken along with the greedy algorithm's bssf
		bbEndTime = time.time() - self.bbStartTime
		
		return self.constructResults(
			self.bssfCost, 
			bbEndTime if self.bbBstsf > bbEndTime else self.bbBstsf,
			self.bbSolutionCount, 
			bssf, 
			self.maxQueueSize, 
			self.totalStates, 
			self.pruneCount)

	# r represents the number of remaining cities (or the depth of the tree)
	def prioritizeByBounds(self, matrix, cities, n, pathSoFar, remainingCityIndeces, parentLowerBound, depth = 1):
		# Create a priority queue of the current level in the tree of cities
		parentCityIndex = pathSoFar[-1]
		# O(r * n^2) time complexity to makeQueue() (or to reduce each of the r matrices)  
		# Q and bounds have O(r) space complexity
		# matrices has O(r * n^2)space complexity 
		Q, bounds, matrices = self.makeQueue(cities, copy.deepcopy(matrix), n, parentCityIndex, parentLowerBound, remainingCityIndeces)

		# O(r) space complexity
		# Create a copy of the remainingIndeces so that we don't lose information
		# as we remove from remainingIndeces
		remainingCityIndecesAtLevel = [*remainingCityIndeces]
		
		# keep track of the max queue size for reporting at the end
		self.runningQueueSize += len(Q)
		if self.maxQueueSize < self.runningQueueSize:
			self.maxQueueSize = self.runningQueueSize
		
		# O(r) time complexity to go through the Q * O(n^2) to drillDown or compare bounds at next level
		while len(Q):
			# stop if the timer ended
			if (time.time() - self.bbStartTime) >= self.timeAllowance:
				self.bbBstsf = self.timeAllowance
				break
			self.totalStates += 1
			nextCity = Q.pop(0)
			lBoundOfNextCity = bounds.pop(0)
			# O(n) time complexity to find the index of an element in the list
			currentCityIndex = remainingCityIndecesAtLevel.index(nextCity._index)

			# Should we prune this branch and its children?
			if lBoundOfNextCity < self.bssfCost:
				# O(r) space complexity
				newRemainingCityIndeces = [ *remainingCityIndecesAtLevel[:currentCityIndex], *remainingCityIndecesAtLevel[currentCityIndex+1:]]
				indexOfCity = nextCity._index
				
				newPathSoFar = [*pathSoFar, indexOfCity]
				localMatrix = copy.deepcopy(matrices.pop(0))
				
				# After experimenting with various numbers, 
				# .9 or slightly higher seems to be a good threshold 
				# for deciding whether to drill down or look at the breadth of the level 
				ratio = (lBoundOfNextCity / self.bssfCost)
				if ratio < .9:
					# O(n^2) time and space (don't include r because it's already accounted for)
					self.prioritizeByBounds(localMatrix, cities, n, [*newPathSoFar], [*newRemainingCityIndeces], lBoundOfNextCity, depth+1)
				else:
					# O(n^2) time and space (don't include r because it's already accounted for)
					self.drillDown(localMatrix, cities, n, [*newPathSoFar], [*newRemainingCityIndeces], lBoundOfNextCity)
			else:
				self.pruneCount += 1
			
			# consider a problem solved once all of its subproblems have been solved
			self.runningQueueSize -= 1

	# O(r * n^2) time complexity to go through the remaining cities and find the lower bound for each city
	# O(n^2) space complexity to store and work with the cost matrix 
	def drillDown(self, matrix, cities, n, pathSoFar, remainingCityIndeces, parentLowerBound):
		remainingCityIndecesAtLevel = [*remainingCityIndeces]
		currentCityIndex = -1

		# Keep track of the max queue size for reporting at the end 
		self.runningQueueSize += len(remainingCityIndeces)
		if self.maxQueueSize < self.runningQueueSize:
			self.maxQueueSize = self.runningQueueSize

		# use the first of the remaining cities as the next path to 'drill down' on

		# O(r * n^2) to get the lower bound and reduce the matrix for each remaining city
		# Use (r) instead of an extra O(n) to iterate through the remaining cities
		#  because the remaining cities reduces with each call and over the whole algorithm
		#  we never go through all n inside drillDown  
		while len(remainingCityIndeces):
			# Stop if the time is up
			if (time.time() - self.bbStartTime) >= self.timeAllowance:
				self.bbBstsf = self.timeAllowance
				return
			self.totalStates += 1

			localMatrix = copy.deepcopy(matrix)
			currentCityIndex += 1

			index = remainingCityIndeces.pop(0)
			# O(n^2) to getLowerBound()
			lb = self.getLowerBound(localMatrix, parentLowerBound, n, pathSoFar[-1], index)
			    
			# Should we prune this branch and its children?
			if lb < self.bssfCost:
				newPathSoFar = [*pathSoFar, index]

				# Have we arrived at a complete path yet?
				if len(newPathSoFar) == n:
					# Can we return to the initial city?
					if cities[index].costTo(cities[pathSoFar[0]]) < math.inf:
						self.bbBstsf = time.time() - self.bbStartTime
						self.bbSolutionCount += 1
						self.bssfCost = lb
						self.bpsf = newPathSoFar
				else:
					newRemainingCityIndeces = [ *remainingCityIndecesAtLevel[:currentCityIndex], *remainingCityIndecesAtLevel[currentCityIndex+1:]]
					# O(n^2) work and O(n^2) space complexity
					self.drillDown(copy.deepcopy(localMatrix), cities, n, [*newPathSoFar], newRemainingCityIndeces, lb)
			else:
				self.pruneCount += 1
			
			# consider a problem solved once all of its subproblems have been solved
			self.runningQueueSize -= 1
	
	# r represent the number of remaining cities, which decreases as we drill deeper in the tree
	# O(r * n^2) time complexity to get the lowerBound for each column in the row
	# O(r * n^2) space complexity to store r matrices of size n^2 
	def makeQueue(self, cities, matrix, n, row, initialCost, remainingCities):
		# the first two lists have O(n) space complexity
		bounds = []
		Q = []
		# O(r * n^2) space complexity
		matrices = []

		# skip the already-visited cities
		# O(r) to iterate through remainingCities * O(n^2) time to getLowerBound()
		for col in remainingCities:
			nLevelMatrix = copy.deepcopy(matrix)
			lb = self.getLowerBound(nLevelMatrix, initialCost, n, row, col)
			
			# If on the first city, add it to the queue and move on to the next city
			if len(bounds) == 0:
				bounds.append(lb)
				Q.append(cities[col])
				matrices.append(nLevelMatrix)
				continue
			
			# O(n) to iterate through the bounds
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
			matrices.insert(i, nLevelMatrix)

		return Q, bounds, matrices

	# O(n) time complexity - to mark all paths along the row and column as None
	# O(n^2) space complexity - to store and work with the n^2 matrix
	def markPathVisited(self, matrix, n, row, col):
		# O(n) to set each element in the row as None 
		self.markDirectionVisited(matrix, n, row)
		# O(n) to set each element in the col as None
		self.markDirectionVisited(matrix, n, col, False)
		matrix[col][row] = None

	def markDirectionVisited(self, matrix, n, index, isRow = True):
		# O(n) time complexity to go through each element in the row/col of length n
		# O(n^2) space complexity to store and work with the n^2 matrix 
		if isRow:	
			for col in range(n):
				matrix[index][col] = None 
		else:
			for row in range(n):
				matrix[row][index] = None

	# O(n^2) time complexity to reduce the matrix
	# O(n^2) space complexity to store and work with the n^2 matrix 
	def getLowerBound(self, matrix, parentLowerBound, n, row, col, ignoreRowCol=True, markVisitedBeforeReduction=True):
		costAtCoordinate = matrix[row][col]
		if markVisitedBeforeReduction:
			# O(n) to mark the paths as visited (None)
			# set all values in the row and column to None (so we don't visit them again) 
			self.markPathVisited(matrix, n, row, col)
		# O(n^2) to reduce the whole matrix
		optimizationCost = self.reduceMatrix(matrix, n, ignoreRowCol, row, col)
		return parentLowerBound + costAtCoordinate + optimizationCost

	# O(n^2) space complexity to store and work with the n^2 matrix
	def reduceMatrix(self, matrix, n, ignoreRowCol=False, ignoreRow=None, ignoreCol=None):
		costToReduce = 0

		# O(n^2) time complexity to iterate through each row (n) and do one or both of the following tasks:
		# 	1. find the smallest element in a row
		# 	2. reduce each element in the row by an amount
		for row in range(n):
			if ignoreRowCol and row == ignoreRow:
				continue
			
			# O(n) time to find the smallest element in the row
			# O(1) space to store the smallest value 
			rowSmallestValue = self.getSmallestElementInDirection(matrix, n, row, True, ignoreCol)

			# 0 < rowSmallestValue < infinity
			if rowSmallestValue and rowSmallestValue < math.inf:
				self.reduceElementsInDirection(matrix, n, rowSmallestValue, row)
			costToReduce += rowSmallestValue
		
		# O(n^2) time complexity to iterate through each col O(n) and do one or both of the following tasks:
		# 	1. find the smallest element in the col O(n)
		# 	2. reduce each element in the col by an amount O(n)
		for col in range(n):
			if ignoreRowCol and col == ignoreCol:
				continue

			# O(n) time to find the smallest element in the col
			# O(1) space to store the smallest value 
			colSmallestValue = self.getSmallestElementInDirection(matrix, n, col, False, ignoreRow)

			# 0 < colSmallestValue < infinity
			if colSmallestValue and colSmallestValue < math.inf:
				# O(n) time to reduce each element in the col
				# space complexity - no extra space allocated beyond the n^2 matrix being updated
				self.reduceElementsInDirection(matrix, n, colSmallestValue, col, False)
			costToReduce += colSmallestValue
		
		return costToReduce

	# O(n) time complexity to go through each element in the row/col of length n
	# O(n^2) space complexity to store and work with the n^2 matrix
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

	# O(n) time to go through n (either the row or the column)
	# O(n^2) space complexity to store and work with the n^2 matrix
	def reduceElementsInDirection(self, matrix, n, reducer, index, isRow = True):
		if isRow:	
			for col in range(n):
				if matrix[index][col] != None and matrix[index][col] < math.inf:
					matrix[index][col] -= reducer 
		else:
			for row in range(n):
				if matrix[row][index] != None and matrix[row][index] < math.inf:
					matrix[row][index] -= reducer

	# Space complexity - O(n^2) to store the 2D matrix
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


