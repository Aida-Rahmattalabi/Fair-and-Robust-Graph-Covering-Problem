#!/usr/bin/python2.7
from __future__ import division
from gurobipy import *
import random
import math

'''
Given seed set X, it finds the worst case fairness of nodes
'''

class find_worst_case_fairness():
    
    def __init__(self, N, M, X, I, J, A, D, fairness_x, fairness_y, alpha_x, alpha_y):
        
        self.N = N  # network size
        self.M = M
        self.X = X  # seed set
        self.I = I  # budget
        self.J = J  # unavilable nodes
        self.A = A  # graph's adjacency matrix
        self.D = D  # nodes categories   rows are categories, columns are members of each category 
        self.fairness_x = fairness_x
        self.fairness_y = fairness_y
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y

        self.theta = 0
        self.x_n = [0 for u in range(self.N)]   # the (worst-case) available nodes
        self.y_n = [0 for u in range(self.M)]   # the (worst-case) covered nodes
        self.c = [0 for u in range(len(self.D))]
        self.v = [[0 for u in range(self.N)] for v in range(len(self.D))]  # linearization variables

        self.m = Model()

        self.y = [0 for u in range(self.M)]
        self.OPT = -1

        self.F_x = [0 for c in range(len(self.D))]
        self.F_y = [0 for c in range(len(self.D))]

        # model parameters
        # self.m.Params.MIPFocus = 2
        # self.m.setParam('TimeLimit', 30*60)
        # self.m.setParam('MIPGap', stopGap)
        # self.m.setParam('Presolve', 1)
        # self.m.setParam('ScaleFlag', 2)

        self.createOptimizeModel()
    

    def addVariables(self):
        print('addeing the variables ...')
        
        for u in range(self.N):
            self.x_n[u] = self.m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "x_"+str(u))  # the nodes that are available 

        for u in range(self.M):                                                              # the coverage
            self.y_n[u] = self.m.addVar(lb = 0, ub = 1, vtype = GRB.CONTINUOUS, name ='y_'+str(u))
        
        self.theta = self.m.addVar(vtype = GRB.CONTINUOUS, name = "theta")

        # a variable that decides if the fairness constraint for group c is violated. 
        for m in range(len(self.D)):
            self.c[m] = self.m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name="c_"+str(m))  # the nodes that are available 


        for m in range(len(self.D)):
            for u in range(self.N):
                if(u in self.D[m]):
                    self.v[m][u] = self.m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "v_"+str(m)+'_'+str(u)) 



    def addConstraints(self):
        
        print('adding the constraints ...')
        print(self.X, sum(self.X))
        for u in range(self.N):
            self.m.addConstr(self.x_n[u] ,GRB.LESS_EQUAL, self.X[u], name ="x_"+str(u))
        
        self.m.addConstr(self.I - self.J, GRB.LESS_EQUAL, quicksum(self.x_n[u] for u in range(self.N)))

        for u in range(self.M):
            for v in range(self.N):
                self.m.addConstr(self.A[v][u]*self.x_n[v], GRB.LESS_EQUAL, self.y_n[u], "y_"+str(u))

        self.m.addConstr(quicksum(self.y_n[u] for u in range(self.M)), GRB.EQUAL, self.theta , "obj")


        for m in range(len(self.D)):
            expr1 = LinExpr(0)
            for u in range(self.N):
                if(u in self.D[m]):
                    expr1.add(self.v[m][u])     
            self.m.addConstr(expr1/len(self.D[m]), GRB.LESS_EQUAL, self.alpha_y - 0.001)

        # add linearization constraints 
        for m in range(len(self.D)):
            for u in range(self.N):
                if(u in self.D[m]):
                    self.m.addConstr(self.v[m][u], GRB.LESS_EQUAL, self.c[m])
                    self.m.addConstr(self.v[m][u], GRB.LESS_EQUAL, self.y_n[u])
                    self.m.addConstr(self.y_n[u]+ self.c[m] - 1, GRB.LESS_EQUAL, self.v[m][u])


        self.m.addConstr(1, GRB.LESS_EQUAL, quicksum(self.c[cc] for cc in range(len(self.D))))


    def addObjective(self):  
        self.m.setObjective(self.theta, GRB.MINIMIZE)  
        
    def createOptimizeModel(self):
        self.addVariables()
        self.addConstraints()
        self.addObjective()
        self.m.update()
        
        '''  
        a = self.m.computeIIS()
        
        for c in self.m.getConstrs():
            if (c.IISConstr > 0):
                print('here:',  c.ConstrName)
    
        '''
        
        self.m.optimize()
        
        if(self.m.status !=3 and self.m.status !=4):

            self.OPT = self.theta.X
            
            for u in range(self.N):
                self.y[u] = round(self.y_n[u].x)


            for c in range(len(self.D)):
                count_y = 0
                count_x = 0
                for u in range(self.N):
                    if(u in self.D[c] and self.y[u]):
                        self.F_y[c] += 1

                    if(u in self.D[c]):
                        count_y += 1

                    if(u in self.D[c] and self.X[u]):
                        self.F_x[c] += 1

                    if(u in self.D[c]):
                        count_x += 1

                self.F_y[c] =  self.F_y[c]/(len(self.D[c]))       
                self.F_x[c] =  self.F_x[c]/count_x       

            for cc in range(len(self.D)):
                print(self.c[cc].X, self.F_y[cc])

    

        print('Done')

