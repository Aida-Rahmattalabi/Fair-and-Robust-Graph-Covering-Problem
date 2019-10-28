#!/usr/bin/python2.7
from __future__ import division

import random
import copy
import time
from guppy import hpy

random.seed(time)  

class baseline():
    def __init__(self, N, M, A, I, J):
        self.N = N
        self.M = M
        self.A = A
        self.I = I
        self.J = J

    def random_selection(self):
        ' selects nodes uniformly at random '
        # initialize 
        X = [0 for i in range(self.N)]

        # random selection 
        while(sum(X) < self.I):
            n = random.randint(0, self.N-1)
            X[n] = 1
        
        return X

    def DC_selection(self):
        ' selects the highest degree nodes '
        
        # initialize
        X = [0 for i in range(self.N)]
        # find the out-degrees for the monitors
        deg = [0 for u in range(N)]
        nodes = {}
        for u in range(self.N):
            deg[u] = sum(self.A[u]) 
            nodes.update({u:deg[u]})

        nodes_sorted = {}
            
        sorted_nodes = sorted(nodes.iteritems(), key=lambda (k,v): (v,k), reverse = True)
        
        # choose the highest-degree nodes
        for u in range(self.I):
            node_to_add = sorted_nodes[u][0]
            X[node_to_add] = 1

        return X

    def greedy_non_robust(self):
        ' greedily chooses the seed set, without considering the uncertainty'

        # initialize 
        S0 = [0 for u in range(self.N)]
        S1 = self.greedy_subroutine(self.I, S0)

        return S1

    def greedy_tzoumas(self): 
        ' this baseline is based on Tzoumas (Jadbababei) paper '
        
        # initialize
        S0 = [0 for u in range(self.N)]
        S1 = [0 for u in range(self.N)]
        
        # find the out-degrees for the monitors
        deg = [0 for u in range(self.N)]
        nodes = {}
        for u in range(self.N):
            deg[u] = sum(self.A[u]) 
            nodes.update({u:deg[u]})

        nodes_sorted = {}
            
        sorted_nodes = sorted(nodes.iteritems(), key=lambda (k,v): (v,k), reverse = True)
        
        # choose the first (robust) set
        for u in range(self.J):
            node_to_add = sorted_nodes[u][0]
            S0[node_to_add] = 1

        # choose the second set
        S1 = self.greedy_subroutine(self.I-self.J, S0)
        S_final = [x + y for x, y in zip(S0, S1)]    

        return S_final


    def baseline_Bogunovic(self):
        ' this baseline is based on Bogunovic paper'
        
        # initialize 
        eta = 1
        
        S0 = [0 for u in range(N)]
        S1 = [0 for u in range(N)]
        
        for i in range(int(math.ceil(math.log(J, 2)))+1): #3
            for j in range(1,int(math.ceil(J/pow(2, i)))+1): 
                res = greedy_subroutine(pow(2,i), S0)
                for u in range(N):
                    if(res[u]):
                        if(sum(S0) < I):
                            S0[u] = 1
                        else:
                            return S0
                
        S1 = greedy_subroutine(I-sum(S0), S0)    
        S_final = [x + y for x, y in zip(S0, S1)]    
     
        return S_final 


    def greedy_subroutine(self, size, ground_set):
        ' greedily chooses seed nodes of size  "size" '
        
        # initialize 
        num_added = 0
        nodes_add = []
        cur_set = [0 for u in range(self.N)]
        
        while(num_added < size):
            nodes_add = []
            cur_coverage = self.coverage(cur_set)
            max_benefit = -1
            
            for u in range(self.N):
                if(ground_set[u] == 0 and cur_set[u]==0):
                    set_temp = copy.deepcopy(cur_set)
                    set_temp[u] = 1
                    benefit = self.coverage(set_temp) - cur_coverage
                    if(benefit >= max_benefit):
                        max_benefit = benefit
            
            for u in range(self.N):
                if(ground_set[u] == 0 and cur_set[u]==0):
                    set_temp = copy.deepcopy(cur_set)
                    set_temp[u] = 1
                    benefit = self.coverage(set_temp) - cur_coverage
                    if(benefit >= max_benefit):
                        nodes_add.append(u)

            random.seed(time)  
            if(len(nodes_add) == 1) :   
                node_to_add = nodes_add[0]
            else: 
                node_to_add = nodes_add[random.randint(0, len(nodes_add)-1)]   
            cur_set[node_to_add] = 1 
            num_added += 1
            
        return cur_set
        
    def coverage(self, set): 
        ' evaluates the coverage of a seed set "set" '
        
        # initialize 
        num_covered = 0

        # counts the number of nodes that have a neighbor in the seed set "set"  
        for u in range(self.M):
            cu = 0        
            for v in range(self.N):
                if(self.A[v][u] and set[v]):
                    cu = 1
            num_covered += cu
            
        return num_covered


