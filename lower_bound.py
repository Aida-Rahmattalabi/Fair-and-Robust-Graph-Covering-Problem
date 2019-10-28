#!/usr/bin/python2.7
from __future__ import division
from gurobipy import *
import random
import math
import operator as op


class lower_bound():
    
    def __init__(self, N, M, I, J, A, fairness_x, fairness_y, alpha_x, alpha_y, node_attribute):
        self.N = N  # Number of target nodes in the network
        self.M = M  # Number of target nodes
        self.I = I  # Total number of invited youth
        self.J = J  # Number of nodes that show up
        
        self.B = self.I - self.J
        
        self.fairness_x = fairness_x
        self.fairness_y = fairness_y
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y

        self.node_attribute = node_attribute        # this vector(list) gives the attributes of each individual
        self.categories_num = max(node_attribute)  # how many categories are there? 
       
        self.A = A
        self.tau = 0
        
        self.x_n = [0 for u in range(self.N)]
        self.y_n = [0 for u in range(self.M)]
    
        self.m = Model()
        
        # model parameters
        # self.m.Params.MIPFocus = 2
        # self.m.setParam('TimeLimit', 30*60)
        # self.m.setParam('MIPGap', stopGap)
        self.m.setParam('Presolve', 1)
        self.m.setParam('OutputFlag', 0 ) 
        self.m.setParam('LogToConsole', 0 )
        self.m.setParam('LogFile', "" ) 

        self.obj = self.createOptimizeModel()

    def addVariables(self):
        
        print('adding the variables...')
        for u in range(self.N):
            self.x_n[u] = self.m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name="xh_"+str(u))
                    
        for u in range(self.M):
            self.y_n[u] = self.m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name ='y_'+str(u))
        
        self.tau = self.m.addVar(vtype=GRB.CONTINUOUS, name = "tau")

    def addConstraints(self):
        print('adding the constraints...')
        
        self.m.addConstr(quicksum(self.x_n[u] for u in range(self.N)), GRB.LESS_EQUAL, self.B, name ="resource_constraint")
        count = 0
        for u in range(self.M):
            self.m.addConstr(self.y_n[u], GRB.LESS_EQUAL, quicksum(self.A[v][u]*self.x_n[v] for v in range(self.N)), "constr_"+str(count))
            count = count + 1
            
        self.m.addConstr(self.tau, GRB.LESS_EQUAL, quicksum(self.y_n[u] for u in range(self.M)), "constr_final")    
        

        # --- add fairness constraints ---
        count = [0 for t in range(self.categories_num)]
        
        for c in range(self.categories_num):
            for u in range(self.N):
                if(self.node_attribute[u] == (c+1)):
                    count[c] += 1  

        print count
        if(self.fairness_x):
            for c in range(self.categories_num):
                expr = LinExpr(0)
                for u in range(self.N):
                    if(self.node_attribute[u] == (c+1)):
                        expr.add(self.x_n[u])

                self.m.addConstr((1/count[c])*expr - self.I/self.N, GRB.LESS_EQUAL, self.alpha_x, name='fariness_x_1_'+str(c))
                self.m.addConstr(-self.alpha_x, GRB.LESS_EQUAL, (1/count[c])*expr - self.I/self.N, name='fariness_x_2_'+str(c))

        if(self.fairness_y):
            for c in range(self.categories_num):
                expr = LinExpr(0)
                expr0 = LinExpr(0)
                for u in range(self.N):
                    expr0.add(self.y_n[u])
                    if(self.node_attribute[u] == (c+1)):
                        expr.add(self.y_n[u])
                self.m.addConstr((1/count[c])*expr - expr0/self.N, GRB.LESS_EQUAL, self.alpha_y, name='fariness_y_1_'+str(c))
                self.m.addConstr(-self.alpha_y, GRB.LESS_EQUAL, (1/count[c])*expr - expr0/self.N, name='fariness_y_2_'+str(c))

    def addObjective(self):  
        self.m.setObjective(self.tau, GRB.MAXIMIZE)  
         
    def createOptimizeModel(self):
        self.addVariables()
        self.addConstraints()
        self.addObjective()
        self.m.update()
        self.m.optimize()

        return self.tau.x   
        
        