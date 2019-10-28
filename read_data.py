#!/usr/bin/python2.7

class read():
    def __init__(self, address, N, M, subgraph):
        self.address = address
        self.N = N
        self.M = M
        self.subgraph = subgraph # choose a subgraph?

        self.ID = {}
    def read_network_data_youthnet(self, net_name):
        ''' ------ returns a 1) dictionary, mapping the PID to a number, and the 2) adjacency matrix for the surveyed people ------ '''
        ''' ----------------------------------------------------'''
        net = open(self.address + net_name)

        self.ID = {}
        count = 0

        first_line = net.readline()

        while True:
            # read line
            line = net.readline()
            # check if line is not empty
            if line: 
                line = line.split(',')
                a = line[0]  # person
                b = line[1]  # friend
                if(a not in self.ID.keys()):
                    self.ID.update({a: count})
                    count += 1 
                if(b not in self.ID.keys()):
                    self.ID.update({b: count})
                    count += 1                     
            if not line:
                break
        net.close()
        ''' -----------------------------------------------------'''
        ''' ------------ create the adjacency matrix ------------'''
        A = [[0 for u in range(len(self.ID))] for v in range(len(self.ID))]
        # read the network (again)
        net = open(self.address + net_name)
        first_line = net.readline()
        count = 0
        while (True):
            # read line
            line = net.readline()
            # check if line is not empty
            if line: 
                line = line.split(',')
                a = line[0]  # person
                b = line[1]  # friend
                count +=1 
                if(b in self.ID.keys()):
                    A[self.ID[a]][self.ID[b]] = 1
                    A[self.ID[b]][self.ID[a]] = 1
            if not line:
                break
        net.close()

        if(self.subgraph):
            AA = [A[i][0:self.M] for i in range(self.N)]
        else: 
            AA = A

        print (len(AA))    
        return AA
        
        ''' ----------------------------------------------------'''
    def read_affliation_data_youthnet(self, aff_name, attribute):

        C = [-1 for i in range(len(self.ID))]

        # returns a list of characteristics
        aff = open(self.address + aff_name)
        first_line = aff.readline()
        
        first_line = first_line.split(',')
        
        count = 0

        while True:
            # read line
            line = aff.readline()
            # check if line is not empty
            if line: 
                line = line.split(',')
                a = line[0]
                if(a in self.ID.keys()):                         
                    if(attribute == 'Gender'):
                        if(line[1]):
                            C[self.ID[line[0]]] = int(line[1])
                        else:
                            C[self.ID[line[0]]] = -1
                    if(attribute == 'Sexori'):
                        if(line[2]):
                            C[self.ID[line[0]]] = int(line[2])
                        else:
                            C[self.ID[line[0]]] = -1
                    if(attribute == 'Race'):
                        if(line[3] and line[3] != '\r\n'):
                            C[self.ID[line[0]]] = int(line[3])
                        else:
                            C[self.ID[line[0]]] = -1
                    if(attribute == 'Suicide_idea'):
                        if(line[4]):
                            C[self.ID[line[0]]] = int(line[4])
                        else:
                            C[self.ID[line[0]]] = -1
                    if(attribute == 'Suicide_plan'):
                        if(line[5]):
                            C[self.ID[line[0]]] = int(line[5])
                        else:
                            C[self.ID[line[0]]] = -1
                    if(attribute == 'Suicide_attempt'):
                        if(line[6]):
                            C[self.ID[line[0]]] = int(line[6])
                        else:
                            C[self.ID[line[0]]] = -1

            if not line:
                break
        aff.close()

        categories_num = max(C)
        if(-1 in C):
            categories_num += 1

        for u in range(len(C)):
            if(C[u] == -1):
                C[u] = categories_num

        



        if(self.subgraph):
            C = [C[i] for i in range(self.N)]

        return C



    def read_nodes_attribute_real(self, D):

        node_attribute = [0 for u in range(self.N)]

        for c in range(len(D)):
            for l in range(len(D[c])):
                u = D[c][l]
                node_attribute[D[c][u]] = c + 1

        return node_attribute



    def define_categories_real(self, C):   #lump groups 
        
        # ----------------- create categories -------------------
        categories_num = max(C)

        categories = [0 for c in range(categories_num)]
        for c in range(categories_num):
            for u in range(self.N):
                if(C[u] == (c+1)):
                    categories[c] += 1


        to_be_merged = []

        for c in range(categories_num):
            if(categories[c] <0.1*self.N):
                to_be_merged.append(c+1)


        node_attribute = [0 for u in range(self.N)]

        count = 1
        old_new = {}
        code = {1:'Other'}
        
        for c in range(len(to_be_merged)):
            old_new.update({(to_be_merged[c]):1})



        for u in range(self.N):
            if(C[u] in old_new.keys()):
                node_attribute[u] = old_new[C[u]]


            if(C[u] not in old_new.keys()):
                old_new.update({C[u]:count+1})
                node_attribute[u] = old_new[C[u]]
                if(C[u] == 1):
                    code.update({count+1: 'American Indian'})
                if(C[u] == 2):
                    code.update({count+1: 'Asian'})
                if(C[u] == 3):
                    code.update({count+1: 'Black/African American'})
                if(C[u] == 4):
                    code.update({count+1: 'Native Hawaiian'})
                if(C[u] == 5):
                    code.update({count+1: 'White'})
                if(C[u] == 6):
                    code.update({count+1: 'Hispanic'})
                if(C[u] == 7):
                    code.update({count+1: 'Mixed'})
                if(C[u] == -1):
                    code.update({count+1: 'Unknown'})


                count += 1


        categories_num = max(node_attribute)
        D = [[]]
        for c in range(categories_num):
            for u in range(self.N):
                if(node_attribute[u] == (c+1)):
                    D[(len(D)-1)].append(u)
            if(c != categories_num-1):
                D.append([])


        #-------------------------------------------------------        
        return node_attribute, D           
                






    def read_community_data(self, address, filename, graphType):
    
    	communities = []
    	graph_file = open(address+"/"+filename,'r')
    	count = 0
    	for row in graph_file:
        	row = row.split(',')
        	communities.append([])

       		#print row
        	if(graphType == 'WS'):
            		pass
        	else:
            		for u in range(len(row)):
                		if(row[u] != '\n'):
                    		  communities[count].append(int(row[u]))
                
            		#print communities
            		count += 1
        
        return communities

	
    def read_SBM_data(self, address, filename, graphType):
    
    	A = [[0 for u in range(self.N)] for v in range(self.N)]
    
    	graph_file = open(address+"/"+filename,'r')
    
    	for row in graph_file:
        	row = row.split(',')
        	if(graphType == 'WS'):
            		if(int(row[0])<= self.N and int(row[1])<=self.N):
                		A[int(row[0])-1][int(row[1])-1] = 1
        	else:
            		if(int(row[0])<= self.N and int(row[1])<=self.N):
                		A[int(row[0])][int(row[1])] = 1
                
    	return A


    def read_nodes_attribute_SBM(self, D):

        node_attribute = [0 for u in range(self.N)]

        for c in range(len(D)):
            for l in range(len(D[c])):
                u = D[c][l]
                node_attribute[D[c][u]] = c + 1

        return node_attribute




