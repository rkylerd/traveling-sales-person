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
		start_time = time.time()

		matrix = self.createCostMatrix(cities, ncities)
		initialCost = self.reduceMatrix(matrix, ncities)
		
		Q, bounds = self.makeQueue(cities, matrix, ncities, initialCost)

		foundTour = False
		count = 0
		self.bssf = None
		btsf_time = None
		self.bssfCost = math.inf

		while (time.time() - start_time) < time_allowance:
			solution_timer = time.time()
			
			while len(Q):
				nextCity = Q.pop(0)
				lBoundOfNextCity = bounds.pop(0)

				if lBoundOfNextCity < self.bssfCost:
					indexOfCity = self.stringToIndex(Q[0]._name)
					remainingCityIndeces = [ *[n for n in range(1,indexOfCity)], *[n for n in range(indexOfCity+1,ncities)]]
					pathSoFar = [0, indexOfCity]
					costOfSolution, path = self.findSolution(matrix, cities, pathSoFar, remainingCityIndeces, bounds[0])

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

	def findSolution(self, matrix, cities, pathSoFar, remainingCityIndeces, parentLowerBound):
		
		for index in remainingCityIndeces:
			lb = self.getLowerBound(matrix,parentLowerBound, pathSoFar[-1], index)
			if lBoundOfNextCity < self.bssfCost:
				self.findSolution(matrix, cities, [*pathSoFar, index], remainingCityIndeces.pop(index), lb)

		if ( len(remainingCityIndeces) == 0 and parentLowerBound < self.bssfCost ):
			self.bssfCost = parentLowerBound
			self.bssf = pathSoFar

	def makeQueue(self, cities, matrix, n, initialCost):
		bounds = []
		Q = []
		row = 0

		for col in range(1,n):
			lb = self.getLowerBound(matrix, initialCost, row, col)
			
			if len(bounds) == 0:
				bounds.append(lb)
				Q.append(cities[col])
				continue

			indexOfNextSmallerThanMe = -1
			for idx, bound in enumerate(bounds):
				if bound > lb:
					break
				# the new lb needs will be placed after the current bound
				indexOfNextSmallerThanMe += 1 
			
			i = 1 + indexOfNextSmallerThanMe
			if i == len(bounds):
				print('what happens when we try to insert here?')
			
			bounds.insert(i, bound)
			Q.insert(i, cities[col])

		return Q, bounds

	def getLowerBound(self, matrix, parentLowerBound, row, col):
		optimizationCost = self.reduceMatrix(matrix)
		return parentLowerBound + matrix[row][col] + optimizationCost

	def reduceMatrix(self, matrix, n):
		costToReduce = 0

		for row in range(n):
			rowSmallestValue = self.getSmallestElementInDirection(matrix, n, row)
			if (rowSmallestValue == math.inf):
				print('rowSmallestValue == math.inf')

			# 0 < rowSmallestValue < infinity
			if rowSmallestValue and rowSmallestValue < math.inf:
				self.reduceElementsInDirection(matrix, n, rowSmallestValue, row)
			costToReduce += rowSmallestValue
		
		for col in range(n):
			colSmallestValue = self.getSmallestElementInDirection(matrix, n, col, False)
			# 0 < colSmallestValue < infinity
			if colSmallestValue and colSmallestValue < math.inf:
				self.reduceElementsInDirection(matrix, n, colSmallestValue, col, False)
			costToReduce += colSmallestValue
		
		return costToReduce

	def getSmallestElementInDirection(self, matrix, n, index, isRow = True):
		n = len(matrix[index])
		smallestValue = math.inf

		if isRow:	
			for col in range(n):
				if matrix[index][col] < smallestValue:
					smallestValue = matrix[index][col]
		else:
			for row in range(n):
				if matrix[row][index] < smallestValue:
					smallestValue = matrix[row][index]
		return smallestValue

	def reduceElementsInDirection(self, matrix, n, reducer, index, isRow = True):

		if isRow:	
			for col in range(n):
				if matrix[index][col] < math.inf:
					matrix[index][col] -= reducer 
		else:
			for row in range(n):
				if matrix[row][index] < math.inf:
					matrix[row][index] -= reducer

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
		



