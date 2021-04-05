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

		while (time.time() - start_time) < time_allowance:
			solution_timer = time.time()
			new_solution = None
			# start from city random city
			city_idx = random.randint(0, len(cities)-1 )
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
				btsf_time = time.time() - solution_timer
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
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		bbStartTime = time.time()
		
		# Get the initial state
		matrix = self.createCostMatrix(cities, ncities)
		initialCost = self.reduceMatrix(matrix, ncities)
		# matrix = [[math.inf,7,3,12],[3,math.inf,6,14],[5,8,math.inf,6],[9,3,5,math.inf]]
		
		# Create a priority queue of the first level in the tree of cities
		Q, bounds = self.makeQueue(cities, copy.deepcopy(matrix), ncities, initialCost)

		self.bpsf = None
		btsf_time = None
		bssf = None
		self.bssfCost = math.inf
		self.bbSolutionCount = 0
		self.bbSolutionTimer = time.time()

		# while first-level states still exist
		while len(Q):
			if (time.time() - bbStartTime) >= time_allowance:
				break

			nextCity = Q.pop(0)
			lBoundOfNextCity = bounds.pop(0)

			# always check bound before trying a path (ordered by lowest bounds)
			if lBoundOfNextCity < self.bssfCost:
				indexOfCity = self.stringToIndex(nextCity._name)
				remainingCityIndeces = [ *[n for n in range(1,indexOfCity)], *[n for n in range(indexOfCity+1,ncities)]]
				
				# We always start from city 0
				pathSoFar = [0, indexOfCity]
				localMatrix = copy.deepcopy(matrix)
				self.markPathVisited(localMatrix, ncities, 0, indexOfCity)
				self.findSolution(copy.deepcopy(localMatrix), cities, ncities, [*pathSoFar], [*remainingCityIndeces], lBoundOfNextCity, bbStartTime, time_allowance)
		route = []		
		for i in self.bpsf:
			route.append(cities[i])

		bssf = TSPSolution(route)

		end_time = time.time()
		results['cost'] = bssf.cost if self.bssfCost < math.inf else math.inf
		results['time'] = self.bbSolutionTimer
		results['count'] = self.bbSolutionCount
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	def markPathVisited(self, matrix, n, row, col):
		self.markDirectionVisited(matrix, n, row)
		self.markDirectionVisited(matrix, n, col, False)
		matrix[col][row] = None

	def findSolution(self, matrix, cities, n, pathSoFar, remainingCityIndeces, parentLowerBound, bbStartTime, timeAllowance):
		remainingCityIndecesInDepth = [*remainingCityIndeces]
		indexInOriginalRemainingCitiesIndeces = -1
		if len(remainingCityIndeces) == 1:
			print(pathSoFar, " --- ", remainingCityIndeces)
		while len(remainingCityIndeces):
			if (time.time() - bbStartTime) >= timeAllowance:
				return
			localMatrix = copy.deepcopy(matrix)
			indexInOriginalRemainingCitiesIndeces += 1
			index = remainingCityIndeces.pop(0)

			lb = self.getLowerBound(localMatrix, parentLowerBound, n, pathSoFar[-1], index, True, True)
			if lb == math.inf:
				continue

			newPathSoFar = [*pathSoFar, index]
			if len(newPathSoFar) == n and lb < self.bssfCost:
				self.bbBstsf = time.time() - self.bbSolutionTimer
				self.bbSolutionCount += 1
				self.bssfCost = lb
				self.bpsf = newPathSoFar
				return

			if self.bpsf == None or lb < self.bssfCost:
				self.findSolution(copy.deepcopy(localMatrix), cities, n, [*newPathSoFar], [ *remainingCityIndecesInDepth[:indexInOriginalRemainingCitiesIndeces], *remainingCityIndecesInDepth[indexInOriginalRemainingCitiesIndeces+1:]], lb, bbStartTime, timeAllowance)
			
	def makeQueue(self, cities, matrix, n, initialCost):
		bounds = []
		Q = []
		row = 0

		# Start at column 1 because the initial state always begins at city 0
		for col in range(1,n):
			lb = self.getLowerBound(copy.deepcopy(matrix), initialCost, n, row, col)
			
			# If on the first city, add it to the queue and move on to the next city
			if len(bounds) == 0:
				bounds.append(lb)
				Q.append(cities[col])
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

		return Q, bounds

	def getLowerBound(self, matrix, parentLowerBound, n, row, col, cancelBeforeReduction=True, ignoreRowCol=False):
		costAtCoordinate = matrix[row][col]
		if cancelBeforeReduction:
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

	def getSmallestElementInDirection(self, matrix, n, index, isRow, ignoreRowCol=None):

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
		return smallestValue if flag else 0

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
		pass
		



