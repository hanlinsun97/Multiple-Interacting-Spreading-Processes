import numpy as np
import networkx as nx
from itertools import product

def framework(G, T):
    for node in G.nodes():
        G.nodes[node]["message_A_star"] = np.zeros([T])
        G.nodes[node]["message_B_star"] = np.zeros([T])
        G.nodes[node]["message_AB_star"] = np.zeros([T])

        G.nodes[node]["A"] = np.zeros([T])
        G.nodes[node]["B"] = np.zeros([T])
        G.nodes[node]["S"] = np.zeros([T])
    for edge in G.edges():
        G.edges[edge]["message_A_star"] = np.zeros([T])
        G.edges[edge]["message_B_star"] = np.zeros([T])
        G.edges[edge]["message_AB_star"] = np.zeros([T])
    return G

def initialization(G, Apos, Bpos):
    for node in G.nodes():
        if node in Apos:
            G.nodes[node]["message_A_star"][0] = 0
            G.nodes[node]["message_B_star"][0] = 1
            G.nodes[node]["message_AB_star"][0] = 0

            G.nodes[node]["A"][0] = 1
            G.nodes[node]["B"][0] = 0
            G.nodes[node]["S"][0] = 0
        elif node in Bpos:
            G.nodes[node]["message_A_star"][0] = 1
            G.nodes[node]["message_B_star"][0] = 0
            G.nodes[node]["message_AB_star"][0] = 0

            G.nodes[node]["A"][0] = 0
            G.nodes[node]["B"][0] = 1
            G.nodes[node]["S"][0] = 0
        else:
            G.nodes[node]["message_A_star"][0] = 0
            G.nodes[node]["message_B_star"][0] = 0
            G.nodes[node]["message_AB_star"][0] = 1

            G.nodes[node]["A"][0] = 0
            G.nodes[node]["B"][0] = 0
            G.nodes[node]["S"][0] = 1
    
    for edge in G.edges():
        i = edge[0]
        if i in Apos:
            G.edges[edge]["message_A_star"][0] = 0
            G.edges[edge]["message_B_star"][0] = 1
            G.edges[edge]["message_AB_star"][0] = 0
        elif i in Bpos:
            G.edges[edge]["message_A_star"][0] = 1
            G.edges[edge]["message_B_star"][0] = 0
            G.edges[edge]["message_AB_star"][0] = 0
        else:
            G.edges[edge]["message_A_star"][0] = 0
            G.edges[edge]["message_B_star"][0] = 0
            G.edges[edge]["message_AB_star"][0] = 1
            
    return G


def ind(tau, t):
    if tau <= t:
        indd = 1
    else:
        indd = 0
    return indd

def q_SA(G, t, setA, setB):
    prod_A = 1
    prod_B = 1
    for j in setA:
        prod_A = prod_A * (1 - alpha_A * ind(G.nodes[j]["tau_A"], t))
    for l in setB:
        prod_B = prod_B * (1 - alpha_B * ind(G.nodes[l]["tau_B"], t))
    numerator = (1 - prod_A) * prod_B
    denominator = prod_A + prod_B - prod_A * prod_B
    qsa = numerator / denominator
    return qsa

def q_SB(G, t, setA, setB):
    prod_A = 1
    prod_B = 1
    for j in setA:
        prod_A = prod_A * (1 - alpha_A * ind(G.nodes[j]["tau_A"], t))
    for l in setB:
        prod_B = prod_B * (1 - alpha_B * ind(G.nodes[l]["tau_B"], t))
    numerator = prod_A * (1 - prod_B)
    denominator = prod_A + prod_B - prod_A * prod_B
    qsb = numerator / denominator
    return qsb

def q_SS(G, t, setA, setB):
    prod_A = 1
    prod_B = 1
    for j in setA:
        prod_A = prod_A * (1 - alpha_A * ind(G.nodes[j]["tau_A"], t))
    for l in setB:
        prod_B = prod_B * (1 - alpha_B * ind(G.nodes[l]["tau_B"], t))
    numerator = prod_A * prod_B
    denominator = prod_A + prod_B - prod_A * prod_B
    qss = numerator / denominator
    return qss

def SetFromConfig(ele, neighbor):  

    #  Config should look like [("A,A"), ("A,B")]. Choose one element such as ("A","A")
    #  Number of states should equal to the number of neighbors.

    setA = [];
    setB = [];
    setS = [];
    for index in range(len(neighbor)):
        if ele[index] == "A":
            setA.append(neighbor[index])
        elif ele[index] == "B":
            setB.append(neighbor[index])  ## TODO: Use a loop. 
        elif ele[index] == "S":
            setS.append(neighbor[index])
    return setA, setB, setS

def message_node(G, GG, node, t):
    neighbor = list(nx.all_neighbors(GG, node))
    neighbor_num = len(neighbor)
    SetConfig = list(product(["A", "B", "S"], repeat=neighbor_num))

    for ele in SetConfig:  # Choose \Theta_A \Theta_B and \Theta_S, the first sum
        setA, setB, setS = SetFromConfig(ele, neighbor)    
        

        # Roll a configuration for {\tau_k^A} for k \in \Theta_A. 
        tau_A_config = list(product(range(t), repeat=len(setA)))
        tau_B_config = list(product(range(t), repeat=len(setB)))
        for eleA in tau_A_config:
            for i in range(len(setA)):
                G.nodes[setA[i]]["tau_A"] = eleA[i]
            for eleB in tau_B_config:
                for j in range(len(setB)):
                    G.nodes[setB[j]]["tau_B"] = eleB[j]
                # continue. 

                # B_star:
                prod_qss_Bstar = 1
                for t_prime in range(t-1):
                    prod_qss_Bstar = prod_qss_Bstar * q_SS(G, t_prime, setA, setB)
                qsa_Bstar = q_SA(G, t-1, setA, setB)
                prod_B_star = prod_qss_Bstar * qsa_Bstar

                # A_star:
                prod_qss_Astar = 1
                for t_prime in range(t-1):
                    prod_qss_Astar = prod_qss_Astar * q_SS(G, t_prime, setA, setB)
                qsb_Astar = q_SB(G, t-1, setA, setB)
                prod_A_star = prod_qss_Astar * qsb_Astar

                # AB_star:
                prod_qss_ABstar = 1
                for t_prime in range(t):
                    prod_qss_ABstar = prod_qss_ABstar * q_SS(G, t_prime, setA, setB)
                prod_AB_star = prod_qss_ABstar

                prod_message_B_star = 1
                prod_message_A_star = 1
                prod_message_AB_star = 1

                for k in setA:
                    
                    prod_message_B_star = prod_message_B_star * \
                        G.edges[k, node]["message_B_star"][G.nodes[k]["tau_A"]]  # Confusing a bit here. tau_A is stored in the network (parameter of a node)
                for l in setB:
                    prod_message_A_star = prod_message_A_star * \
                        G.edges[l, node]["message_A_star"][G.nodes[l]["tau_B"]]
                for n in setS:
                    prod_message_AB_star = prod_message_AB_star * \
                        G.edges[n, node]["message_AB_star"][t-1]
                prod_thesame = prod_message_B_star * prod_message_A_star * prod_message_AB_star
                
                G.nodes[node]["message_B_star"][t] += prod_B_star * prod_thesame
                G.nodes[node]["message_A_star"][t] += prod_A_star * prod_thesame
                G.nodes[node]["message_AB_star"][t] += prod_AB_star * prod_thesame

        G.nodes[node]["message_B_star"][t] *= G.nodes[node]["S"][0] 
        G.nodes[node]["message_A_star"][t] *= G.nodes[node]["S"][0]
        G.nodes[node]["message_AB_star"][t] *= G.nodes[node]["S"][0]
    return G


def message_edge(G, GG, edge, t):
    neighbor_old = list(nx.all_neighbors(GG, edge[0]))
    neighbor = []
    for ele in neighbor_old:
        if ele != edge[1]:
            neighbor.append(ele)
    neighbor_num = len(neighbor)
    SetConfig = list(product(["A", "B", "S"], repeat=neighbor_num))

    for ele in SetConfig:  # Choose \Theta_A \Theta_B and \Theta_S, the first sum
        setA, setB, setS = SetFromConfig(ele, neighbor)

        # Roll a configuration for {\tau_k^A} for k \in \Theta_A.
        tau_A_config = list(product(range(t), repeat=len(setA)))
        tau_B_config = list(product(range(t), repeat=len(setB)))
      
        for eleA in tau_A_config:
            for i in range(len(setA)):
                G.nodes[setA[i]]["tau_A"] = eleA[i]
            for eleB in tau_B_config:
                for j in range(len(setB)):
                    G.nodes[setB[j]]["tau_B"] = eleB[j]
                # continue.

                # B_star:
                prod_qss_Bstar = 1
                for t_prime in range(t-1):
                    prod_qss_Bstar = prod_qss_Bstar * \
                        q_SS(G, t_prime, setA, setB)
                qsa_Bstar = q_SA(G, t-1, setA, setB)
                prod_B_star = prod_qss_Bstar * qsa_Bstar
            
                # A_star:
                prod_qss_Astar = 1
                for t_prime in range(t-1):
                    prod_qss_Astar = prod_qss_Astar * \
                        q_SS(G, t_prime, setA, setB)
                qsb_Astar = q_SB(G, t-1, setA, setB)
                prod_A_star = prod_qss_Astar * qsb_Astar

                # AB_star:
                prod_qss_ABstar = 1
                for t_prime in range(t):
                    prod_qss_ABstar = prod_qss_ABstar * \
                        q_SS(G, t_prime, setA, setB)
                prod_AB_star = prod_qss_ABstar

                # THE LAST THREE PRODUCTS. 

                prod_message_B_star = 1
                prod_message_A_star = 1
                prod_message_AB_star = 1

                for k in setA:
                    prod_message_B_star = prod_message_B_star * \
                        G.edges[k, edge[0]]["message_B_star"][G.nodes[k]["tau_A"]]  # Confusing a bit here. tau_A is stored in the network (parameter of a node)
                    
                for l in setB:
                    prod_message_A_star = prod_message_A_star * \
                        G.edges[l, edge[0]]["message_A_star"][G.nodes[l]["tau_B"]]
                    
                for n in setS:
                    prod_message_AB_star = prod_message_AB_star * \
                        G.edges[n, edge[0]]["message_AB_star"][t-1]
                prod_thesame = prod_message_B_star * prod_message_A_star * prod_message_AB_star
                
                
                G.edges[edge]["message_B_star"][t] += prod_B_star * prod_thesame
                G.edges[edge]["message_A_star"][t] += prod_A_star * prod_thesame
                G.edges[edge]["message_AB_star"][t] += prod_AB_star * prod_thesame

        
        G.edges[edge]["message_B_star"][t] *= G.nodes[edge[0]]["S"][0]
        G.edges[edge]["message_A_star"][t] *= G.nodes[edge[0]]["S"][0]
        G.edges[edge]["message_AB_star"][t] *= G.nodes[edge[0]]["S"][0]
    return G


def marginal(G, node, t):
    G.nodes[node]["A"][t] = G.nodes[node]["A"][t-1] + G.nodes[node]["message_B_star"][t]
    G.nodes[node]["B"][t] = G.nodes[node]["B"][t-1] + G.nodes[node]["message_A_star"][t]
    G.nodes[node]["S"][t] = G.nodes[node]["message_AB_star"][t]
    return G

######################################################################################################
########################################### Start the main ###########################################
######################################################################################################


GG = nx.random_tree(5, seed=1)
G = GG.to_directed()

alpha_A = 0.4
alpha_B = 0.6
Apos = [2]
Bpos = [3]
target = 0
T = 6
G = framework(G, T)
G = initialization(G, Apos, Bpos)


for t in range(1, T):
    for edge in G.edges():
        G = message_edge(G, GG, edge, t)

    for node in G.nodes():
        G = message_node(G, GG, node, t)
    
    for node in G.nodes():
        G = marginal(G, node, t)

print("S", G.nodes[target]["S"])
print("A", G.nodes[target]["A"])
print("B", G.nodes[target]["B"])

