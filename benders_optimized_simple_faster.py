#!/usr/bin/python2.7

'''
********************************************************
* The model shows a Benders Decomposition for a the K-adaptability MILP in the network monitoring problem.
* The original MILP is decomposed into two problems. 


Additional comment: 
1) Lexicographic constraints are used to speed-up the solution
#********************************************************
'''

from __future__ import division

# Standard packages
from gurobipy import *
import copy
import random
import time
import itertools
from itertools import permutations

class network_monitoring:
    
    def __init__(self, N, M, I, J, A, K, LB, X_init, Y_init, givenlimit, communities, fairness_x, fairness_y, node_attribute, alpha_x, alpha_y):

        self.N = N  # number of nodes in the network
        self.M = M  # number of second stage binary variables (i.e. number of target nodes)
        self.I = I  # total number of invited youth
        self.J = J  # number of people unavailable
        self.A = [row[0:self.M] for row in A[0:self.N]]  # the graph adjacency matrix
        self.K  = K
        self.TAU = LB    
        self.LB = LB                                     # lower Bound
        self.X = X_init
        self.Y = Y_init
        self.communities = communities
        self.fairness_x = fairness_x
        self.fairness_y = fairness_y
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y
        self.node_attribute = node_attribute        # this vector(list) gives the attributes of each individual
        self.categories_num = max(node_attribute)   # how many categories are there? 
        self.X_init = X_init
        self.L  = self.M                                 # number of second stage constraints
        self.elements = [i for i in range(self.L+1)]     # number of constraints + 1 -> to form the set L={0, ..., L}^{K}
        self.Q  = self.N + 1                             # number of uncertain parameters
        self.R  = 2*self.N + 3                           # number of constraints that define the uncertainty set
        self.epsilon = 1
        self.BM = 10000
        self.sm = Model('SubProblem')
        #self.sm.setParam('Method', 1.0)           
        self.sm.setParam( 'OutputFlag', 0 ) 
        self.sm.setParam( 'LogToConsole', 0 )
        self.sm.setParam( 'LogFile', "" )  
        #self.sm.setParam('FeasibilityTol', 10**(-2))
        #self.sm.params.threads = 1
        #self.sm.params.NodefileStart = 0.5
        #self.sm.params.timeLimit = 30*60
        self.sm.params.DualReductions = 0                # turn off to avoid the optimization status of INF_OR_UNBD
        self.sm.params.InfUnbdInfo = 1                   # additional info for infeasible/unbounded models
        #self.sm.ModelSense = -1                         # sub problem is a maximization problem
        #self.sm.setParam('IntFeasTol', 10**(-2))

        self.mm = Model('MasterProblem')
        #self.mm.setParam( 'OutputFlag', 0 ) 
        #self.mm.setParam( 'LogToConsole', 0 )
        #self.mm.setParam( 'LogFile', "" )  
        self.mm.setParam('MIPFocus', 3)
        #self.mm.setParam('IntFeasTol', 0.0001)
        self.mm.setParam('MIPGap', 0.05)
        #self.mm.setParam('timeLimit', 2*60)
        #self.mm.setParam('FeasibilityTol', 10**(-4))
        #self.mm.setParam('IntFeasTol', 10**(-2))
        #self.mm.params.threads = 1
        #self.mm.params.NodefileStart = 0.5
        #self.mm.params.timeLimit = 1800
        self.mm.params.DualReductions = 0                # turn off to avoid the optimization status of INF_OR_UNBD
        self.mm.params.InfUnbdInfo = 1                   # additional info for infeasible/unbounded models
        
        # Master problem variables
        self.tau = [0]  
        self.x_n   = [0 for u in range(self.N)] 
        self.y_k_n = [[0 for u in range(self.M)] for v in range(self.K)]

        self.create_master_problem()
        self.mm.optimize()  
         
        c = self.mm.getVarByName('tau')
        self.TAU = c.x 
        for u in range(M):
            for k in range(K):
                c = self.mm.getVarByName('y_%s_%s' %(k,u))
                self.Y[k][u] = c.x         
        for u in range(N):
            c = self.mm.getVarByName('x_%s' %(u))
            self.X[u] = c.x

        # Matrices 
        self.Qc = [[0 for u in range(self.M)] for v in range(self.Q)]
        self.W =  [[0 for u in range(self.M)] for v in range(self.L)]
        self.Ae = [[0 for u in range(self.Q)] for v in range(self.R)]
        self.b =  [0 for u in range(self.R)]
        self.H_new = [[0 for u in range(self.Q)] for v in range(self.L)]

        self.populate_matrices()

        self.constr_list = list(itertools.combinations_with_replacement(self.elements, K))   
        #self.constr_list = list(tuple(itertools.product(self.elements, repeat = K))) 
        self.perms =  list(permutations([i for i in range(self.K)]))
        self.l_p = 1 #len(self.perms)

    def populate_matrices(self):

        for u in range(self.R):
            for v in range(self.N+1):
                if(u==0):
                    if(v!=self.N):
                        self.Ae[u][v] = -1  
                if(u>=1 and u<=(self.N+1)):
                    self.Ae[u][u-1] = 1
                if(u==self.N+2): 
                    if(v==self.N):
                        self.Ae[u][v] = -1
                if(u>self.N+2):
                    self.Ae[u][u-(self.N+3)] = -1
        for u in range(self.R):
            if(u==0):
                self.b[u] = -(self.N-self.J)
            if(u>=1 and u<= (self.N+1)):
                self.b[u] = 1
            if(u==self.N+2): 
                self.b[u] = -1
            if(u>(self.N+2)):
                self.b[u] = 0
        for u in range(self.M):
            self.Qc[self.N][u] = -1
            self.W[u][u] = 1
        for u in range(self.M):
            for v in range(self.N+1):
                if(v != self.N):
                    self.H_new[u][v] = self.A[v][u]     

    def Sub_problem(self, lk, X, Y, TAU):
    
        self.sm.reset()
        for v in self.sm.getVars():
            self.sm.remove(v)
        for c in self.sm.getConstrs():
            self.sm.remove(c) 
            
        epsilon_x = [[[[0 for u in range(self.N)] for k in range(self.K)] for l in range(self.L)] for t in range(4)]
        epsilon_y = [[[[0 for u in range(self.M)] for k in range(self.K)] for l in range(self.L)] for t in range(4)]
        epsilon_x_gamma = [[[0 for u in range(self.N)] for k in range(self.K)] for t in range(4)]
        epsilon_y_gamma = [[[0 for u in range(self.M)] for k in range(self.K)] for t in range(4)]
        epsilon_y_lambda = [[[0 for u in range(self.M)] for k in range(self.K)] for t in range(4)]
            
        #------------------------------------ dual variables for linearization constraints -----------------------------------------------
        if(0 in self.constr_list[lk]):    
            for l in range(self.L):
                for k in range(self.K):
                    for u in range(self.N):                                                                             
                        epsilon_x[0][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_x_0_"+str(u)+'_'+str(k)+'_'+str(l))
                        epsilon_x[1][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = (1-X[u])*self.BM,vtype=GRB.CONTINUOUS, name="epsilon_x_1_"+str(u)+'_'+str(k)+'_'+str(l))
                        epsilon_x[2][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = X[u]*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_x_2_"+str(u)+'_'+str(k)+'_'+str(l))
                        epsilon_x[3][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_x_3_"+str(u)+'_'+str(k)+'_'+str(l))

                    for u in range(self.M):
                        epsilon_y[0][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_y_0_"+str(u)+'_'+str(k)+'_'+str(l))
                        epsilon_y[1][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = (1 - Y[k][u])*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_y_1_"+str(u)+'_'+str(k)+'_'+str(l))
                        epsilon_y[2][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = Y[k][u]*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_y_2_"+str(u)+'_'+str(k)+'_'+str(l))
                        epsilon_y[3][l][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_y_3_"+str(u)+'_'+str(k)+'_'+str(l))

        for k in range(self.K):
            for u in range(self.N): 
                epsilon_x_gamma[0][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_x_gamma_0_"+str(u)+"_"+str(k))
                epsilon_x_gamma[1][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = (1-X[u])*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_x_gamma_1_"+str(u)+"_"+str(k))
                epsilon_x_gamma[2][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj =  X[u]*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_x_gamma_2_"+str(u)+"_"+str(k))
                epsilon_x_gamma[3][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_x_gamma_3_"+str(u)+"_"+str(k))

            for u in range(self.M):
                epsilon_y_gamma[0][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_y_gamma_0_"+str(u)+"_"+str(k))
                epsilon_y_gamma[1][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = (1 - Y[k][u])*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_y_gamma_1_"+str(u)+"_"+str(k))
                epsilon_y_gamma[2][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = Y[k][u]*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_y_gamma_2_"+str(u)+"_"+str(k))
                epsilon_y_gamma[3][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_y_gamma_3_"+str(u)+"_"+str(k)) 
                         
            if(0 in self.constr_list[lk]):     
                for u in range(self.M):
                    epsilon_y_lambda[0][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_y_lambda_0_"+str(u)+"_"+str(k))
                    epsilon_y_lambda[1][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = (1 - Y[k][u])*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_y_lambda_1_"+str(u)+"_"+str(k))
                    epsilon_y_lambda[2][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = Y[k][u]*self.BM, vtype=GRB.CONTINUOUS, name="epsilon_y_lambda_2_"+str(u)+"_"+str(k))
                    epsilon_y_lambda[3][k][u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_y_lambda_3_"+str(u)+"_"+str(k)) 

        # Add constraints
        if(1):
            if (0 in self.constr_list[lk]):
                self.epsilon_u = [0 for u in range(self.Q+self.K+2)]       
                #------------------------------------ Other dual variables (other constraints) -----------------------------------------------
                for u in range(self.Q+self.K+2): 
                    if(u == self.Q):
                        self.epsilon_u[u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = TAU, vtype=GRB.CONTINUOUS, name="epsilon_"+str(u))
                    else:
                        if(u == self.Q+self.K+1):
                            self.epsilon_u[u] = self.sm.addVar(lb = -GRB.INFINITY, ub = GRB.INFINITY, obj = 1, vtype=GRB.CONTINUOUS, name="epsilon_"+str(u))
                        else:
                            self.epsilon_u[u] = self.sm.addVar(lb = -GRB.INFINITY, ub = GRB.INFINITY, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_"+str(u))
                
                self.sm.update()

                # Gamma k 
                for k in range(self.K):
                    expr1 = LinExpr(0)
                    expr2 = LinExpr(0)
                    expr3 = LinExpr(0)
                    
                    expr1.add(quicksum((epsilon_x_gamma[1][k][u] - epsilon_x_gamma[3][k][u]) for u in range(self.N)))
                    expr2.add(quicksum((epsilon_y_gamma[1][k][u] - epsilon_y_gamma[3][k][u]) for u in range(self.N)))
                    if(self.constr_list[lk][k] != 0):
                        expr3.add(-self.epsilon*self.epsilon_u[self.Q])                
                    self.sm.addConstr(expr1+expr2+expr3 , GRB.LESS_EQUAL, 0, name = 'gamma_k_'+str(k))
                    
                # Gamma_x k n
                for k in range(self.K):
                    for u in range(self.N):
                        expr1 = LinExpr(0)
                        expr1.add(-epsilon_x_gamma[0][k][u]-epsilon_x_gamma[1][k][u]+epsilon_x_gamma[2][k][u]+epsilon_x_gamma[3][k][u]) 
                        if(self.constr_list[lk][k] != 0):
                            l = self.constr_list[lk][k]
                            l = l-1          
                            if(self.H_new[l][u] != 0):
                                if(u != (self.Q-1)):
                                    expr1.add(self.H_new[l][u]*self.epsilon_u[u])   
                        self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0, name = 'gamma_x_k_u_'+str(k)+'_'+str(u))

                # Gamma_y k n
                for k in range(self.K):
                    for u in range(self.M):
                        expr1 = LinExpr(0)
                        expr1.add(-epsilon_y_gamma[0][k][u]-epsilon_y_gamma[1][k][u]+epsilon_y_gamma[2][k][u]+epsilon_y_gamma[3][k][u]) 
                        if(self.constr_list[lk][k] != 0):
                            l = self.constr_list[lk][k]
                            l = l-1          
                            if(self.W[l][u] != 0):
                                if(u != (self.Q-1)):
                                    expr1.add(self.W[l][u]*self.epsilon_u[self.Q])   
                        self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0, name = 'gamma_y_k_u_'+str(k)+'_'+str(u))

                # Lambda k
                for k in range(self.K):
                    expr1 = LinExpr(0)
                    expr1.add(quicksum((epsilon_y_lambda[1][k][u] - epsilon_y_lambda[3][k][u]) for u in range(self.M)))
                    expr1.add(self.epsilon_u[self.Q+self.K+1])
                    if(self.constr_list[lk][k] != 0):
                        expr1.add(self.epsilon_u[self.Q+k+1])        
                    self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0, name = 'lambda_'+str(k))

                # y_k_lambda k n 
                for k in range(self.K):
                    for u in range(self.M):           
                        expr1 = LinExpr(0)
                        expr2 = LinExpr(0)
                        expr1.add(-epsilon_y_lambda[0][k][u]-epsilon_y_lambda[1][k][u]+epsilon_y_lambda[2][k][u]+epsilon_y_lambda[3][k][u])
                        expr2.add(quicksum(self.Qc[v][u]*self.epsilon_u[v] for v in range(self.Q)))
                        self.sm.addConstr(expr1-expr2, GRB.LESS_EQUAL, 0, name = 'lambda_y_'+str(k)+'_'+str(u))
    
                # beta_l_l_k    
                for k in range(self.K):
                    for l in range(self.L):
                        expr1 = LinExpr(0)
                        expr1.add(quicksum((epsilon_y[1][l][k][u] - epsilon_y[3][l][k][u]) for u in range(self.M)))
                        expr1.add(quicksum((epsilon_x[1][l][k][u] - epsilon_x[3][l][k][u]) for u in range(self.N))) 
                        self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0, name = 'beta_l_k'+str(k)+'_'+str(l))
                        
                # beta_x
                for k in range(self.K):
                    for l in range(self.L):
                        for u in range(self.N):
                            expr1 = LinExpr(0)
                            expr1.add(-epsilon_x[0][l][k][u]-epsilon_x[1][l][k][u]+epsilon_x[2][l][k][u]+epsilon_x[3][l][k][u])
                            if(self.constr_list[lk][k] == 0):
                                expr1.add(-self.H_new[l][u]*self.epsilon_u[u])
                                
                            self.sm.addConstr(expr1 , GRB.LESS_EQUAL, 0)
                        
                # beta_y_k
                for k in range(self.K):
                    for l in range(self.L):
                        for u in range(self.M):
                            expr1 = LinExpr(0)
                            expr1.add(-epsilon_y[0][l][k][u]-epsilon_y[1][l][k][u]+epsilon_y[2][l][k][u]+epsilon_y[3][l][k][u])
                            if(self.constr_list[lk][k] == 0):
                                expr1.add(-self.W[l][u]*self.epsilon_u[self.Q])
                                
                            self.sm.addConstr(expr1 , GRB.LESS_EQUAL, 0)
                            
                # alpha_r 
                for r in range(self.R): 
                    expr1 = LinExpr(0)
                    expr1.add(self.b[r]*self.epsilon_u[self.Q])
                    expr1.add(quicksum((self.Ae[r][u]*self.epsilon_u[u]) for u in range(self.Q)))
                    
                    self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0)
                    
                self.sm.update()
                
            else:
                self.epsilon_u = [0 for u in range(self.Q+1)]       
                #------------------------------------ Other dual variables (other constraints) -----------------------------------------------
                for u in range(self.Q+1): 
                    if(u == self.Q):
                        self.epsilon_u[u] = self.sm.addVar(lb = -GRB.INFINITY, ub = 0, obj = -1, vtype=GRB.CONTINUOUS, name="epsilon_"+str(u))
                    else:
                        self.epsilon_u[u] = self.sm.addVar(lb = -GRB.INFINITY, ub = GRB.INFINITY, obj = 0, vtype=GRB.CONTINUOUS, name="epsilon_"+str(u))
		self.sm.update()
                # Gamma k 
                for k in range(self.K):
                    expr1 = LinExpr(0)
                    expr2 = LinExpr(0)
                    expr3 = LinExpr(0)
                    
                    expr1.add(quicksum((epsilon_x_gamma[1][k][u] - epsilon_x_gamma[3][k][u]) for u in range(self.N)))
                    expr2.add(quicksum((epsilon_y_gamma[1][k][u] - epsilon_y_gamma[3][k][u]) for u in range(self.N)))
                    if(self.constr_list[lk][k] != 0):
                        expr3.add(-self.epsilon*self.epsilon_u[self.Q])                
                    self.sm.addConstr(expr1+expr2+expr3 , GRB.LESS_EQUAL, 0, name = 'gamma_k_'+str(k))
                    
                # Gamma_x k n
                for k in range(self.K):
                    for u in range(self.N):
                        expr1 = LinExpr(0)
                        expr1.add(-epsilon_x_gamma[0][k][u]-epsilon_x_gamma[1][k][u]+epsilon_x_gamma[2][k][u]+epsilon_x_gamma[3][k][u]) 
                        l = self.constr_list[lk][k]
                        l = l-1          
                        if(self.H_new[l][u] != 0):
                            if(u != (self.Q-1)):
                                expr1.add(self.H_new[l][u]*self.epsilon_u[u])   
                        self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0, name = 'gamma_x_k_u_'+str(k)+'_'+str(u))

                # Gamma_y k n
                for k in range(self.K):
                    for u in range(self.M):
                        expr1 = LinExpr(0)
                        expr1.add(-epsilon_y_gamma[0][k][u]-epsilon_y_gamma[1][k][u]+epsilon_y_gamma[2][k][u]+epsilon_y_gamma[3][k][u]) 
                        l = self.constr_list[lk][k]
                        l = l-1          
                        if(self.W[l][u] != 0):
                            if(u != (self.Q-1)):
                                expr1.add(self.W[l][u]*self.epsilon_u[self.Q])   
                        self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0, name = 'gamma_y_k_u_'+str(k)+'_'+str(u))

                # alpha_r 
                for r in range(self.R): 
                    expr1 = LinExpr(0)
                    expr1.add(self.b[r]*self.epsilon_u[self.Q])
                    expr1.add(quicksum((self.Ae[r][u]*self.epsilon_u[u]) for u in range(self.Q)))
                    self.sm.addConstr(expr1, GRB.LESS_EQUAL, 0)
                   
                   
        # Add Objective
        expr0 = LinExpr(0)
        if(0 in self.constr_list[lk]):    
            expr0.add(self.epsilon_u[self.Q]*TAU)
            expr0.add(self.epsilon_u[self.Q+self.K+1]*1)
        else:
            expr0.add(self.epsilon_u[self.Q]*(-1))

        
        for l in range(self.L):
            for k in range(self.K):
                for u in range(self.N):                                                                             
                    expr0.add(epsilon_x[1][l][k][u]*(1-X[u])*self.BM)
                    expr0.add(epsilon_x[2][l][k][u]*(X[u])*self.BM)

                for u in range(self.M):
                    expr0.add(epsilon_y[1][l][k][u]*(1 - Y[k][u])*self.BM)
                    expr0.add(epsilon_y[2][l][k][u]*Y[k][u]*self.BM)
                                
        for k in range(self.K):
            for u in range(self.N): 
                expr0.add(epsilon_x_gamma[1][k][u]*(1-X[u])*self.BM)
                expr0.add(epsilon_x_gamma[2][k][u]*X[u]*self.BM)
    
            for u in range(self.M):
                expr0.add(epsilon_y_gamma[1][k][u]*(1 - Y[k][u])*self.BM)
                expr0.add(epsilon_y_gamma[2][k][u]*(Y[k][u])*self.BM)
                         
            if(0 in self.constr_list[lk]):     
                for u in range(self.M):
                    expr0.add(epsilon_y_lambda[1][k][u]*(1 - Y[k][u])*self.BM)
                    expr0.add(epsilon_y_lambda[2][k][u]*(Y[k][u])*self.BM)
         
        self.sm.setObjective(expr0, GRB.MAXIMIZE)     
        self.sm.update()
                 
   
    def create_master_problem(self):   

        print('creating the master problem ...')
        self.tau[0] = self.mm.addVar(lb=-self.N, ub=0, vtype=GRB.CONTINUOUS, name = 'tau')
        
        for u in range(self.N):
            self.x_n[u] = self.mm.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name = 'x_'+str(u))
            self.mm.update()
	    self.x_n[u].start = self.X_init[u]         
        for u in range(self.M):
            for k in range(self.K):
                self.y_k_n[k][u] = self.mm.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name='y_'+str(k)+'_'+str(u))
                #self.y_k_n[k][u].BranchPriority = u

        self.mm.update()

        self.mm.addConstr(quicksum(self.x_n[u] for u in range(self.N)), GRB.EQUAL, self.I, 'constr_0') 
        self.mm.addConstr(self.LB, GRB.LESS_EQUAL, self.tau[0], 'constr_1')  
        
        # --- ADD FAIRNESS CONSTRAINTS TYPE ONE  ---
	
        if(self.fairness_x):
            for c in range(self.categories_num):
                expr = LinExpr(0)
                count = 0
                for u in range(self.N):
                    if(self.node_attribute[u] == (c+1)):
                        count += 1 
                        expr.add(self.x_n[u])

                self.mm.addConstr(expr/count - self.I/self.N, GRB.LESS_EQUAL, self.alpha_x, name = 'fairness_x_1_'+str(c))
                self.mm.addConstr(-self.alpha_x, GRB.LESS_EQUAL, expr/count - self.I/self.N, name = 'fairness_x_2_'+str(c))



        if(self.fairness_y):
            for c in range(self.categories_num):
                for k in range(self.K):

                    expr = LinExpr(0)
                    count = 0
                    for u in range(self.N):
                        if(self.node_attribute[u] == (c+1)):
                            count += 1 
                            expr.add(self.y_k_n[k][u])

                    self.mm.addConstr(self.alpha_y, GRB.LESS_EQUAL, expr/count, name = 'fairness_y_1_'+str(c)+'_'+str(k))
                    #self.mm.addConstr(-self.alpha_y, GRB.LESS_EQUAL, expr/count - self.tau[0]/self.N, name = 'fairness_y_2_'+str(c)+'_'+str(k))

	
        # --- ADD FAIRNESS CONSTRAINTS TYPE TWO  ---
        '''
        if(self.fairness_x):
            for c in range(self.categories_num):
                for c_p in range(self.categories_num):
                    expr = LinExpr(0)
                    count = 0
                    if(c_p != c):
                        expr_p = LinExpr(0)
                        count_p = 0
                        for u in range(self.N):

                            if(self.node_attribute[u] == (c+1)):
                                count += 1 
                                expr.add(self.x_n[u])

                        for u in range(self.N):
                            if(self.node_attribute[u] == (c_p+1)):
                                count_p += 1 
                                expr_p.add(self.x_n[u])

                    	self.mm.addConstr(expr/count - expr_p/count_p, GRB.LESS_EQUAL, self.alpha_x, name = 'fairness_x_1_'+str(c))
                    	self.mm.addConstr(-self.alpha_x, GRB.LESS_EQUAL, expr/count - expr_p/count_p, name = 'fairness_x_2_'+str(c))



        if(self.fairness_y):

            for c in range(self.categories_num):
                for c_p in range(self.categories_num):
                    for k in range(self.K):
                        expr = LinExpr(0)
                        count = 0
                        if(c_p != c):
                            expr_p = LinExpr(0)
                            count_p = 0
                            for u in range(self.N):

                                if(self.node_attribute[u] == (c+1)):
                                    count += 1 
                                    expr.add(self.y_k_n[k][u])

                            for u in range(self.N):
                                if(self.node_attribute[u] == (c_p+1)):
                                    count_p += 1 
                                    expr_p.add(self.y_k_n[k][u])

                            self.mm.addConstr(expr/count - expr_p/count_p, GRB.LESS_EQUAL, self.alpha_y, name = 'fairness_y_1_'+str(c))
                            self.mm.addConstr(-self.alpha_y, GRB.LESS_EQUAL, expr/count - expr_p/count_p, name = 'fairness_y_2_'+str(c))
           
                       
        '''
        '''
        # --------- Add lexicographic constraints --------
        # add an auxiliary variable
        
        self.d_n_k_kp = [[[0 for u in range(self.M)] for k in range(self.K)] for kp in range(self.K)]
        for u in range(self.M):
            for k in range(self.K):
                for kp in range(self.K):
                    self.d_n_k_kp[kp][k][u] = self.mm.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name = 'd_'+str(u)+'_'+str(k)+'_'+str(kp)) # not sure if they are necessary to be bianry
	
        # add necesssary constraints for these variables
        for k in range(self.K):
            for kp in range(self.K):
                for u in range(self.M):
                    self.mm.addConstr(self.y_k_n[k][u] - self.y_k_n[kp][u], GRB.LESS_EQUAL, self.d_n_k_kp[kp][k][u])
                    self.mm.addConstr(self.y_k_n[kp][u] - self.y_k_n[k][u], GRB.LESS_EQUAL, self.d_n_k_kp[kp][k][u])



        for k in range(self.K):
            for kp in range(self.K):
                for u in range(self.M):
                    self.mm.addConstr(self.d_n_k_kp[kp][k][u] , GRB.LESS_EQUAL, self.y_k_n[k][u] + self.y_k_n[kp][u] )
                    self.mm.addConstr(self.d_n_k_kp[kp][k][u] , GRB.LESS_EQUAL, 2 - (self.y_k_n[k][u] + self.y_k_n[kp][u]))

                    #self.mm.addConstr(self.y_k_n[kp][u] - self.y_k_n[k][u], GRB.LESS_EQUAL, self.d_n_k_kp[kp][k][u])

        # add lexicographic constraints
        for k in range(self.K):
            for kp in range(self.K):
                if(k < kp):
                    for u in range(self.M):
                        expr = quicksum(self.d_n_k_kp[k][kp][h] for h in range(u-1))
                        self.mm.addConstr(self.y_k_n[k][u] - expr, GRB.LESS_EQUAL, self.y_k_n[kp][u])

        '''

        # --------- Add lexicographic constraints new constraint --------
        # add an auxiliary variable
        
        self.d_n_k_kp = [[[0 for u in range(self.M)] for k in range(self.K)] for kp in range(self.K)]
        for u in range(self.M):
            for k in range(self.K):
                for kp in range(self.K):
                    self.d_n_k_kp[kp][k][u] = self.mm.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name = 'd_'+str(u)+'_'+str(k)+'_'+str(kp)) # not sure if they are necessary to be bianry
        self.mm.update()
        # add necesssary constraints for these variables
        for k in range(self.K):
            for kp in range(self.K):
                for u in range(self.M):
                    self.mm.addConstr(self.y_k_n[k][u] - self.y_k_n[kp][u], GRB.LESS_EQUAL, self.d_n_k_kp[kp][k][u])
                    self.mm.addConstr(self.y_k_n[kp][u] - self.y_k_n[k][u], GRB.LESS_EQUAL, self.d_n_k_kp[kp][k][u])



        for k in range(self.K):
            for kp in range(self.K):
                for u in range(self.M):
                    self.mm.addConstr(self.d_n_k_kp[kp][k][u] , GRB.LESS_EQUAL, self.y_k_n[k][u] + self.y_k_n[kp][u] )
                    self.mm.addConstr(self.d_n_k_kp[kp][k][u] , GRB.LESS_EQUAL, 2 - (self.y_k_n[k][u] + self.y_k_n[kp][u]))

                    #self.mm.addConstr(self.y_k_n[kp][u] - self.y_k_n[k][u], GRB.LESS_EQUAL, self.d_n_k_kp[kp][k][u])

        # add lexicographic constraints
        for k in range(self.K):
            kp = k+1
            if(kp < self.K):
                for u in range(self.M):
                    expr = quicksum(self.d_n_k_kp[k][kp][h] for h in range(0, u-1))
                    self.mm.addConstr(self.y_k_n[kp][u] - expr, GRB.LESS_EQUAL, self.y_k_n[k][u])


        '''
        # additional cuts ------
        self.zz = [0 for u in range(self.M)]
        for u in range(self.M):
            self.zz[u] = self.mm.addVar(lb = 0, ub = 1, vtype=GRB.BINARY, name = 'zz_'+str(u))
                   
        for u in range(self.M):
            self.mm.addConstr(quicksum(self.x_n[v]*self.A[v][u] for v in range(self.N))-self.J, GRB.LESS_EQUAL, (self.I-self.J)*self.zz[u])
            self.mm.addConstr(quicksum(self.y_k_n[k][u] for k in range(self.K)), GRB.LESS_EQUAL, self.K-self.zz[u])
        '''

        for k in range(self.K):
            self.mm.addConstr(-quicksum(self.y_k_n[k][u] for u in range(self.M)), GRB.LESS_EQUAL, self.tau[0])    
        
        
        for k in range(self.K):
            for u in range(self.M):
                self.mm.addConstr(self.y_k_n[k][u], GRB.LESS_EQUAL, quicksum(self.A[n][u]*self.x_n[n] for n in range(self.N)))


        self.mm.setObjective(self.tau[0], GRB.MINIMIZE)  
        print('CREATED THE MASTER PROBLEM')

    def Master_problem(self, epsilon_u_v, epsilon_x_v, epsilon_y_v, epsilon_x_gamma_v, epsilon_y_gamma_v, epsilon_y_lambda_v, 
                       epsilon_u_c, epsilon_x_c, epsilon_y_c, epsilon_x_gamma_c, epsilon_y_gamma_c, epsilon_y_lambda_c, lk, LB):

        # reset the Master Problem
        # self.mm.reset()

        for pp in range(self.l_p):
            # add more cuts
            expr = LinExpr(0)
            if (0 in self.constr_list[lk]):
                for c in range(len(epsilon_u_v)):
                    u = epsilon_u_v[c]
                    if(u == self.Q):
                        expr.add(epsilon_u_c[c]*self.tau[0])

                        # print('epsilon_u, u: ', epsilon_u_v[c], self.Q, epsilon_u_c[c])
                    if(u == self.Q+self.K+1):
                        expr.add(epsilon_u_c[c])
                        # print('epsilon_u, u: ', epsilon_u_v[c], self.Q+self.K+1, epsilon_u_c[c])

            else:
                for c in range(len(epsilon_u_v)):
                    u = epsilon_u_v[c]
                    if(u == self.Q):
                        expr.add(epsilon_u_c[c]*(-1))

                        # print('epsilon_u, u: ', epsilon_u_v[c], self.Q, epsilon_u_c[c])
                        
            for c in range(len(epsilon_x_v)):
                u = epsilon_x_v[c][0]
                k = epsilon_x_v[c][1]
                #k = self.perms[pp][k]
                if(epsilon_x_v[c][3] == 1):
                    expr.add(epsilon_x_c[c]*(1-self.x_n[u])*self.BM)

                if(epsilon_x_v[c][3] == 2):
                    expr.add(epsilon_x_c[c]*self.x_n[u]*self.BM)
                                 
            for c in range(len(epsilon_y_v)):
                u = epsilon_y_v[c][0]
                k = epsilon_y_v[c][1]
                #k = self.perms[pp][k]

                if(epsilon_y_v[c][3] == 1):
                    expr.add(epsilon_y_c[c]*(1-self.y_k_n[k][u])*self.BM)

                if(epsilon_y_v[c][3] == 2):
                    expr.add(epsilon_y_c[c]*self.y_k_n[k][u]*self.BM)
                                
            for c in range(len(epsilon_x_gamma_v)):
                u = epsilon_x_gamma_v[c][0]
                k = epsilon_x_gamma_v[c][1]
                #k = self.perms[pp][k]

                if(epsilon_x_gamma_v[c][2] == 1):
                    expr.add(epsilon_x_gamma_c[c]*(1-self.x_n[u])*self.BM)

                if(epsilon_x_gamma_v[c][2] == 2):
                    expr.add(epsilon_x_gamma_c[c]*self.x_n[u]*self.BM)
                  
            for c in range(len(epsilon_y_gamma_v)):
                u = epsilon_y_gamma_v[c][0]
                k = epsilon_y_gamma_v[c][1]
                #print(self.perms)
                #k = self.perms[pp][k]

                if(epsilon_y_gamma_v[c][2] == 1):
                    expr.add(epsilon_y_gamma_c[c]*(1-self.y_k_n[k][u])*self.BM)

                if(epsilon_y_gamma_v[c][2] == 2):
                    expr.add(epsilon_y_gamma_c[c]*self.y_k_n[k][u]*self.BM)
                         
            if(0 in self.constr_list[lk]):
                for c in range(len(epsilon_y_lambda_v)):
                    u = epsilon_y_lambda_v[c][0]
                    k = epsilon_y_lambda_v[c][1]
                    #k = self.perms[pp][k]

                    if(epsilon_y_lambda_v[c][2] == 1):
                        expr.add(epsilon_y_lambda_c[c]*(1-self.y_k_n[k][u])*self.BM)

                    if(epsilon_y_lambda_v[c][2] == 2):
                        expr.add(epsilon_y_lambda_c[c]*self.y_k_n[k][u]*self.BM)
                              
            self.mm.addConstr(expr, GRB.LESS_EQUAL, 0)
            #self.mm.addConstr(LB, GRB.LESS_EQUAL, self.tau[0])    

        self.mm.update()


def Benders_decomposition(N, M, I, J, A, K, obj_lb, X_init, Y_init, givenlimit, communities, fairness_x, fairness_y, C, alpha_x, alpha_y):
      
    L  = N              # Number of second stage constraints
    LK = (L+1)**K       # Size of the set L={0, ..., L}^{K}
    Q  = N + 1          # Number of uncertain parameters
    R  = 2*N + 3        # Number of constraints that define the uncertainty set
    epsilon = 1
    elements = [i for i in range(L+1)]      # Number of constraints + 1 -> to form the set L={0, ..., L}^{K}
    Z = [0, 0]


    #constr_list = list(tuple(itertools.product(elements, repeat = K)))   
    constr_list = list(itertools.combinations_with_replacement(elements, K))   



    CMMP = network_monitoring(N, M, I, J, A, K, obj_lb, X_init, Y_init, givenlimit, communities, fairness_x, fairness_y, C, alpha_x, alpha_y)

    reset_num = 1
     
    X = CMMP.X
    Y = CMMP.Y
    TAU = CMMP.TAU
        
    # Initialization
    nCut = 0

    lk = len(constr_list) - 1  # Changed this line it was lk = 0
    
        
    #epsilon_u_c = [u]  
    #epsilon_x_c = [(n,k,l,t)]
    #epsilon_y_c = [(m,k,l,t)]
    #epsilon_x_gamma_c = [(n,k,t)]
    #epsilon_y_gamma_c = [(m,k,t)]
    #epsilon_y_lambda_c = [(m,k,t)]

    s_time = time.time()
    epsilon_u_v = []
    epsilon_x_v = []
    epsilon_y_v = []
    epsilon_x_gamma_v = []
    epsilon_y_gamma_v = []
    epsilon_y_lambda_v = []
       
    epsilon_u_c = []
    epsilon_x_c = []
    epsilon_y_c = []
    epsilon_x_gamma_c = []
    epsilon_y_gamma_c = []
    epsilon_y_lambda_c = []
        
    should_continue = True  
    count = 0  

    # Benders decomposition procedure
    duration = time.time() - s_time
    while should_continue and (duration < givenlimit):

        duration = time.time() - s_time
        print "\n Iteration %s" % (nCut + 1)
        if(nCut > 0.2*(LK)):
            reset_num = 1#N/2
        if(nCut > 0.3*(LK)):
            reset_num = 1#N/4
        if(nCut > 0.6*(LK)):
            reset_num = 1

        epsilon_u_v = []
        epsilon_x_v = []
        epsilon_y_v = []
        epsilon_x_gamma_v = []
        epsilon_y_gamma_v = []
        epsilon_y_lambda_v = []
           
        epsilon_u_c = []
        epsilon_x_c = []
        epsilon_y_c = []
        epsilon_x_gamma_c = []
        epsilon_y_gamma_c = []
        epsilon_y_lambda_c = []
            
        if (lk < 0  and count == 0):
            should_continue = False
            lk = -1
        if(lk<0 and count == 1): 
            lk =  len(constr_list) - 1 
            count = 0
            should_continue = True
        if(lk == -1):
            break    
        cut_added =False


        yes = False

        if(sum(CMMP.constr_list[lk]) > 0):#and (0 not in CMMP.constr_list[lk])):
            for k in range(K):
                if(Y[k][CMMP.constr_list[lk][0]-1] == 1 and (CMMP.constr_list[lk][0]>0)):
                    yes = True

            if(yes):#(Y[0][CMMP.constr_list[lk][0]-1]==1 and (CMMP.constr_list[lk][0]>0)) or (Y[1][CMMP.constr_list[lk][1]-1]==1 and (CMMP.constr_list[lk][1]>0))):# or (Y[2][CMMP.constr_list[lk][2]-1]==1 and (CMMP.constr_list[lk][2]>0))):  
                CMMP.Sub_problem(lk, X, Y, TAU)
                #try:
                print('LK: ', lk)
                print('look: ', CMMP.constr_list[lk])
                CMMP.sm.optimize()
                if CMMP.sm.status == 5:
                    count = 1
                    should_continue = True
                    print "Unbounded"
                    nCut += 1
                    if(0 in CMMP.constr_list[lk]):
                        for u in range(N + 1 + K + 1 + 1): 
                            var = CMMP.sm.getVarByName('epsilon_%s' % (u))
                            if(var.getAttr('UnbdRay') != 0.0):
                                epsilon_u_v.append(u)
                                epsilon_u_c.append(var.getAttr('UnbdRay'))
                                #print('epsilon_u (u, value): ', u, var.getAttr('UnbdRay'))
        
                        for u in range(N): 
                            for k in range(K): 
                                for l in range(L): 
                                    for t in range(4): 
                                        var = CMMP.sm.getVarByName('epsilon_x_%s_%s_%s_%s' % (t, u, k, l))
                                        if(var.getAttr('UnbdRay') != 0):
                                            epsilon_x_v.append([u, k, l, t])
                                            epsilon_x_c.append(var.getAttr('UnbdRay'))
                                            #print('epsilon_x (u, value): ', (u, k, l, t), var.getAttr('UnbdRay'))
        
                        for u in range(M): 
                            for k in range(K): 
                                for l in range(L): 
                                    for t in range(4): 
                                        var = CMMP.sm.getVarByName('epsilon_y_%s_%s_%s_%s' % (t, u, k, l))
                                        if(var.getAttr('UnbdRay') != 0):
                                            epsilon_y_v.append([u, k, l, t])
                                            epsilon_y_c.append(var.getAttr('UnbdRay'))
                                            #print('epsilon_y (u, value): ', (u, k, l, t), var.getAttr('UnbdRay'))
        
                        for u in range(N): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_x_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_x_gamma_v.append([u, k, t])
                                        epsilon_x_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_x_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                             
                        for u in range(M): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_y_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_y_gamma_v.append([u, k, t])
                                        epsilon_y_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_y_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
        
        
                        for u in range(M): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_y_lambda_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_y_lambda_v.append([u, k, t])
                                        epsilon_y_lambda_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_y_lambda (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                                                            
        
                    else:
                        for u in range(N + 1 + 1): 
                            var = CMMP.sm.getVarByName('epsilon_%s' % (u))
                            if(var.getAttr('UnbdRay') != 0.0):
                                epsilon_u_v.append(u)
                                epsilon_u_c.append(var.getAttr('UnbdRay'))
                                #print('epsilon_u (u, value): ', u, var.getAttr('UnbdRay'))
                                                     
                        for u in range(N): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_x_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_x_gamma_v.append([u, k, t])
                                        epsilon_x_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_x_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                                                       
                        for u in range(M): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_y_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_y_gamma_v.append([u, k, t])
                                        epsilon_y_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_y_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                        
                    '''
                    except:
                        print "Sub problem solve error"
                        break
                    '''
                                        
                    #create master problems
                    
                    CMMP.Master_problem(copy.copy(epsilon_u_v), copy.copy(epsilon_x_v), copy.copy(epsilon_y_v), copy.copy(epsilon_x_gamma_v), copy.copy(epsilon_y_gamma_v), copy.copy(epsilon_y_lambda_v), 
                               copy.copy(epsilon_u_c), copy.copy(epsilon_x_c), copy.copy(epsilon_y_c), copy.copy(epsilon_x_gamma_c), copy.copy(epsilon_y_gamma_c), copy.copy(epsilon_y_lambda_c), lk, -sum(Y[0]))           
                    cut_added =True
        else:

            for k in range(K):
                if(Y[k][CMMP.constr_list[lk][0]-1] == 1 and (CMMP.constr_list[lk][0]>0)):
                    yes = True

            if(yes):

                CMMP.Sub_problem(lk, X, Y, TAU)
                #try:
                print('LK: ', lk)
                print('look: ', CMMP.constr_list[lk])
                CMMP.sm.optimize()
                if CMMP.sm.status == 5:
                    count = 1
                    should_continue = True
                    print "Unbounded"
                    nCut += 1
                    if(0 in CMMP.constr_list[lk]):
                        for u in range(N + 1 + K + 1 + 1): 
                            var = CMMP.sm.getVarByName('epsilon_%s' % (u))
                            if(var.getAttr('UnbdRay') != 0.0):
                                epsilon_u_v.append(u)
                                epsilon_u_c.append(var.getAttr('UnbdRay'))
                                #print('epsilon_u (u, value): ', u, var.getAttr('UnbdRay'))
        
                        for u in range(N): 
                            for k in range(K): 
                                for l in range(L): 
                                    for t in range(4): 
                                        var = CMMP.sm.getVarByName('epsilon_x_%s_%s_%s_%s' % (t, u, k, l))
                                        if(var.getAttr('UnbdRay') != 0):
                                            epsilon_x_v.append([u, k, l, t])
                                            epsilon_x_c.append(var.getAttr('UnbdRay'))
                                            #print('epsilon_x (u, value): ', (u, k, l, t), var.getAttr('UnbdRay'))
        
                        for u in range(M): 
                            for k in range(K): 
                                for l in range(L): 
                                    for t in range(4): 
                                        var = CMMP.sm.getVarByName('epsilon_y_%s_%s_%s_%s' % (t, u, k, l))
                                        if(var.getAttr('UnbdRay') != 0):
                                            epsilon_y_v.append([u, k, l, t])
                                            epsilon_y_c.append(var.getAttr('UnbdRay'))
                                            #print('epsilon_y (u, value): ', (u, k, l, t), var.getAttr('UnbdRay'))
        
                        for u in range(N): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_x_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_x_gamma_v.append([u, k, t])
                                        epsilon_x_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_x_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                             
                        for u in range(M): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_y_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_y_gamma_v.append([u, k, t])
                                        epsilon_y_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_y_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
        
        
                        for u in range(M): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_y_lambda_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_y_lambda_v.append([u, k, t])
                                        epsilon_y_lambda_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_y_lambda (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                                                            
        
                    else:
                        for u in range(N + 1 + 1): 
                            var = CMMP.sm.getVarByName('epsilon_%s' % (u))
                            if(var.getAttr('UnbdRay') != 0.0):
                                epsilon_u_v.append(u)
                                epsilon_u_c.append(var.getAttr('UnbdRay'))
                                #print('epsilon_u (u, value): ', u, var.getAttr('UnbdRay'))
                                                     
                        for u in range(N): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_x_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_x_gamma_v.append([u, k, t])
                                        epsilon_x_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_x_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                                                       
                        for u in range(M): 
                            for k in range(K): 
                                for t in range(4): 
                                    var = CMMP.sm.getVarByName('epsilon_y_gamma_%s_%s_%s' % (t, u, k))
                                    if(var.getAttr('UnbdRay') != 0):
                                        epsilon_y_gamma_v.append([u, k, t])
                                        epsilon_y_gamma_c.append(var.getAttr('UnbdRay'))
                                        #print('epsilon_y_gamma (u, value): ', (u, k, t), var.getAttr('UnbdRay'))
                                        
                   
                    # Update the master problem
                    CMMP.Master_problem(copy.copy(epsilon_u_v), copy.copy(epsilon_x_v), copy.copy(epsilon_y_v), copy.copy(epsilon_x_gamma_v), copy.copy(epsilon_y_gamma_v), copy.copy(epsilon_y_lambda_v), 
                               copy.copy(epsilon_u_c), copy.copy(epsilon_x_c), copy.copy(epsilon_y_c), copy.copy(epsilon_x_gamma_c), copy.copy(epsilon_y_gamma_c), copy.copy(epsilon_y_lambda_c), lk, -sum(Y[0]))           
            
                    cut_added =True

        if((((nCut+1) / reset_num - (nCut+1) // reset_num) == 0) and cut_added) :
            #print((nCut+1) / reset_num , (nCut+1) // reset_num)
            CMMP.mm.optimize()        
            c = CMMP.mm.getVarByName('tau')
            TAU = round(c.x)   
            for u in range(M):
                for k in range(K):
                    c = CMMP.mm.getVarByName('y_%s_%s' %(k,u))
                    Y[k][u] = round(c.x)    
            for u in range(N):
                c = CMMP.mm.getVarByName('x_%s' %(u))
                X[u] = round(c.x)
            
            Z[0] = CMMP.mm.NumVars
            Z[1] =  CMMP.mm.NumConstrs

            print 'Solved the MASTER PROBLEM AND THIS IS THE SOLUTION: ', X, sum(Y[0]), TAU
        lk = lk - 1    
    print('Done')    
    return X, Y, TAU, Z








    


    
    
