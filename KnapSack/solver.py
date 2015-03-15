#!/usr/bin/python
# -*- coding: utf-8 -*-

from pulp import *
import numpy as np
from scipy.optimize import minimize
from cvxopt import matrix
from cvxopt.modeling import op, dot, variable

from collections import namedtuple
Item = namedtuple("Item", ['index', 'value', 'weight'])

def solve_it(input_data):

    # parse the input
    lines = input_data.split('\n')

    firstLine = lines[0].split()
    item_count = int(firstLine[0])
    capacity = int(firstLine[1])

    items = []

    for i in range(1, item_count+1):
        line = lines[i]
        parts = line.split()
        items.append(Item(i-1, int(parts[0]), int(parts[1])))
    
    # The method preferred, if there are multiple:  
    method = 7
    
    # Greedy - 1
    if method == 0:
        # a trivial greedy algorithm for filling the knapsack
        # it takes items in-order until the knapsack is full
        value = 0
        weight = 0
        taken = [0]*len(items)
    
        for item in items:
            if weight + item.weight <= capacity:
                taken[item.index] = 1
                value += item.value
                weight += item.weight
                
    # Greedy - 2                
    elif method == 1:
        # Another greedy algorithm for filling the knapsack
        # Create a list based on value in descending order
        # Start picking from the top of the list which has the highest value
        # Continue until knapsack is full
        value = 0
        weight = 0
        taken = [0]*len(items)
        
        newItemList = []
        newItemList = sorted(items, key=lambda item: item.value, reverse=True)
        #print newItemList
        
        for item in newItemList:
            #print item.index, item.value, item.weight
            if weight + item.weight <= capacity:
                taken[item.index] = 1
                value += item.value
                weight += item.weight        

    # Greedy - 3 
    elif method == 2:
        # Another greedy algorithm for filling the knapsack
        # Create a list based on weight, sort them in ascending order
        # Start picking from the top of the list which has the lowest weight
        # Continue until knapsack is full
        value = 0
        weight = 0
        taken = [0]*len(items)
        
        newItemList = []
        newItemList = sorted(items, key=lambda item: item.weight, reverse=False)
        #print newItemList
        
        for item in newItemList:
            #print item.index, item.value, item.weight
            if weight + item.weight <= capacity:
                taken[item.index] = 1
                value += item.value
                weight += item.weight

    # Greedy - 4             
    elif method == 3:
        # Another greedy algorithm for filling the knapsack
        # Create a list based on value per weight, sort them in descending order
        # Start picking from the top of the list which has the highest value per unit weight
        # Continue until knapsack is full
        value = 0
        weight = 0
        taken = [0]*len(items)
        
        newItemList = []
        newItemList = sorted(items, key=lambda item: (float(item.value)/float(item.weight)), reverse=True)
        #print "Capacity", capacity
        #print newItemList
        
        for item in newItemList:
            #print item.index, item.value, item.weight
            if weight + item.weight <= capacity:
                taken[item.index] = 1
                value += item.value
                weight += item.weight

    # Branch and Bound                
    elif method == 4:
        # First we need to determine an optimistic estimate of the best solution to the subproblem in order to have an upper bound
        # We will use linear relaxation to determine the upper bound
            # Order items based on value per unit weight
            # Fill knapsack with items starting from the highest value per unit weight 
            # If some space is left at the end, take a fraction of the next item in the list and calculate the value
            # This value is the upper bound
        # Then apply depth first branch and bound
        
        value = 0
        weight = 0
        taken = [0]*len(items)
        
        newItemList = []
        newItemList = sorted(items, key=lambda item: (float(item.value)/float(item.weight)), reverse=True)
        
        count = 0
        upperBound = 0
        tempWeight = 0
        for item in newItemList:
            if tempWeight + item.weight <= capacity:
                upperBound += item.value
                tempWeight += item.weight
            else:
                if count == 0:
                    upperBound += (capacity - tempWeight) * item.value / item.weight
                count += 1
        
        for level in range(len(items)):
            print level + 1
                

    # Dynamic Programming
    elif method == 5:
        # Dynamic Programming
        # We need to create a table with row number equal to the capacity and column number equal to the item count.
        # In the largest problem the capacity is 1MM and item count is 10K. 
        # The memory requirement is approx 1MM * 10 K * 8 / 1e9 = 80 GB - which is huge. 
        
        value = 0
        weight = 0
        taken = [0]*len(items)
        
        prevDict = {}
        currDict = {}
        
        for j in xrange(len(items)+1):
            for i in range(capacity + 1):
                if j == 0:
                    currDict[i] = 0
                elif items[j-1].weight <= i:
                    currDict[i] = max(prevDict[i], items[j-1].value + prevDict[i-items[j-1].weight])
                else:
                    currDict[i] = prevDict[i]
            prevDict = currDict.copy() 
            currDict = {}
                
        print prevDict
        print currDict


    # Use Pulp-OR
    elif method == 6:

        value = 0
        weight = 0
        taken = [0]*len(items)
        #Tuning parameter
        #capacity = capacity + 1
        
        #print(items)
        
        itemList = [ str(itm) for itm in xrange(len(items))]
        valueDict = {}
        weightDict = {}
        
        for i, idx in enumerate(itemList):
            valueDict[idx] = items[i].value
            weightDict[idx] = items[i].weight
        
        #print itemList
        #print valueDict
        #print weightDict
        
        prob = LpProblem("Knapsack Problem", LpMaximize)
        
        item_vars = LpVariable.dicts("x", itemList , lowBound = 0, upBound = 1, cat = LpInteger)
        
        prob += lpSum([item_vars[itm]*valueDict[itm] for itm in itemList]), "Sum of values"
        
        prob += lpSum([item_vars[itm]*weightDict[itm] for itm in itemList]) <= capacity
        
        prob.writeLP("Knapsack.lp")
        
        prob.solve()
        
        #print "Status:", LpStatus[prob.status]
        
        #for variable in prob.variables():
        #    print variable.name, "=", variable.varValue
            
        #print "Total value of the knapsack = ", pulp.value(prob.objective)
        
        
        #value = int(pulp.value(prob.objective))
        
        for i, itm in enumerate(itemList):
            if item_vars[itm].varValue == 1.0:
                #print item_vars[itm].name, item_vars[itm].varValue
                #print i, items[i]
                taken[i] = 1
                value += items[i].value
                

    # Use CVXOPT
    elif method == 7:

        value = 0
        weight = 0
        taken = [0]*len(items)
        
        x = variable(4)
        aList = [[] for _ in xrange(len(items))]
        for i,item in enumerate(items):
            aList[i].append(float(item.weight))
            aList[i].append(float(0))
            aList[i].append(float(0))
        A = matrix(aList)
        b = matrix([capacity,0,0])
        c = matrix([ -item.value for item in items])
        ineq = ( A * x <= b )
        lp2 = op(dot(c,x), ineq)
        #lp2.solve()
        #print(lp2.objective.value())
        
        
    else:
        value = 0
        weight = 0
        taken = [0]*len(items)
    
    # prepare the solution in the specified output format
    output_data = str(value) + ' ' + str(1) + '\n'
    output_data += ' '.join(map(str, taken))
    return output_data


import sys

if __name__ == '__main__':

    if len(sys.argv) == 1:
        file_location = "./data/ks_4_0"
        input_data_file = open(file_location, 'r')
        input_data = ''.join(input_data_file.readlines())
        input_data_file.close()
        print solve_it(input_data) 
    elif len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        input_data_file = open(file_location, 'r')
        input_data = ''.join(input_data_file.readlines())
        input_data_file.close()
        print solve_it(input_data)       
    else:
        print 'This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/ks_4_0)'

