#!/usr/bin/python2.7
from __future__ import division
from gurobipy import *
import random
import math

'''
Given a solution X (selected set of nodes), it finds the worst case coverage of nodes
'''

class find_worst_case():
    
    def __init__(self, N, M, J, gamma, A, fairness_x, fairness_y, node_attribute, alpha_x, alpha_y):
        
        self.N = N  # Number of nodes in the network
        self.M = M
        self.X = [1 for u in range(self.N)]  # Input solution
        self.I = self.N  # Budget
        self.J = J  # Number of nodes that are unavailable
        self.gamma = gamma # Proportion that do not show up
        if(gamma>=0):
            self.T = math.floor(self.gamma*self.I) # Number of nodes that are unavailable
        
        self.A = [row[0:M] for row in A[0:N]]  # The graph adjacency matrix
        
        self.theta = 0
        self.x_n = [0 for u in range(self.N)]   # The (worst-case) available nodes
        self.y_n = [0 for u in range(self.M)]   # The (worst-case) covered nodes
        self.OPT = -1
        self.categories_num = max(node_attribute)
        self.fairness_x = fairness_x
        self.fairness_y = fairness_y
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y

        self.node_attribute = node_attribute        # this vector(list) gives the attributes of each individual
        self.m = Model()
        self.Y = [0 for u in range(self.M)]
        # model parameters
        self.m.Params.MIPFocus = 2
        self.m.setParam('TimeLimit', 30*60)
        #self.m.setParam('MIPGap', stopGap)
        #self.m.setParam('Presolve', 1)
        #self.m.setParam('ScaleFlag', 2)

        self.createOptimizeModel()
    
    def random_node_failure(self):
        
        count = 0    
        x_temp = [i for i,x in enumerate(self.X) if(x == 1)]
        while (count <= self.J):
            u = random.randint(0, len(x_temp)-1)
            print u, count, x_temp, self.X
            
            if(self.X[x_temp[u]] == 1.0):
                self.X[x_temp[u]] = 0
                count = count + 1 
        print('failed nodes: ', self.J, count)


    def addVariables(self):
        print('added the variables')
        for u in range(self.N):
            self.x_n[u] = self.m.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name="x_"+str(u))
        self.theta = self.m.addVar(vtype=GRB.CONTINUOUS, name="theta")

        for u in range(self.M):
            self.y_n[u] = self.m.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name ='y_'+str(u))
        
    def addConstraints(self):
        
        print('added the constraints')
        
        for u in range(self.N):
            self.m.addConstr(self.x_n[u] ,GRB.LESS_EQUAL, self.X[u], name ="x_"+str(u))
        
        self.m.addConstr(self.I - self.J, GRB.LESS_EQUAL, quicksum(self.x_n[u] for u in range(self.N)))
                         
        '''                 
        for u in range(self.N):
            for v in range(self.N):
                for w in range(self.N):
                    self.m.addConstr(self.y_n[u], GRB.LESS_EQUAL, 0.5*(self.A[w][u]*self.x_n[w]+self.A[v][u]*self.x_n[v]), "y_"+str(u))
        '''
        
        for u in range(self.M):
            for v in range(self.N):
                self.m.addConstr(self.A[v][u]*self.x_n[v], GRB.LESS_EQUAL, self.y_n[u], "y_"+str(u))

        self.m.addConstr(quicksum(self.y_n[u] for u in range(self.M)), GRB.EQUAL, self.theta , "obj")
        
        '''

        if(1):
            expr = LinExpr(0)
            for u in range(self.N):
                if(u in [12, 17, 18, 19]):
                    print('u', u)
                    expr.add(self.y_n[u])
            self.m.addConstr(2, GRB.LESS_EQUAL, expr)

        '''
        # --- add fairness constraints ---
        count = [0 for t in range(self.categories_num)]
        
        for c in range(self.categories_num):
            for u in range(self.N):
                if(self.node_attribute[u] == (c+1)):
                    count[c] += 1  

        '''
        print count
        if(self.fairness_x):
            for c in range(self.categories_num):
                expr = LinExpr(0)
                for u in range(self.N):
                    if(self.node_attribute[u] == (c+1)):
                        expr.add(self.x_n[u])

                self.m.addConstr((1/count[c])*expr - self.I/self.N, GRB.LESS_EQUAL, self.alpha_x, name='fariness_x_1_'+str(c))
                self.m.addConstr(-self.alpha_x, GRB.LESS_EQUAL, (1/count[c])*expr - self.I/self.N, name='fariness_x_2_'+str(c))
        '''
        if(self.fairness_y):
            for c in range(self.categories_num):
                expr = LinExpr(0)
                expr0 = LinExpr(0)
                for u in range(self.N):
                    expr0.add(self.y_n[u])
                    if(self.node_attribute[u] == (c+1)):
                        expr.add(self.y_n[u])

                #self.m.addConstr((1/count[c])*expr - expr0/self.N, GRB.LESS_EQUAL, self.alpha_y, name='fariness_y_1_'+str(c))
                self.m.addConstr(self.alpha_y, GRB.LESS_EQUAL, (1/count[c])*expr, name='fariness_y_2_'+str(c))

    def addObjective(self):  
        self.m.setObjective(self.theta, GRB.MINIMIZE)  
        
    def createOptimizeModel(self):
        self.addVariables()
        self.addConstraints()
        self.addObjective()
        self.m.update()
        self.m.optimize()

        self.OPT = self.theta.X
        
        for u in range(self.N):
            self.Y[u] = self.y_n[u].x
        '''
        Final = [0 for u in range(self.N)]
        for u in range(self.N):
            print(u, self.x_n[u].X)
            if(self.x_n[u]):
                Final = [x + y for x, y in zip(Final, self.A[u])]  
        print(Final)  
        '''
        
        print('Done')

