
import numpy as np
import networkx as nx

p_A = 0.4
p_B = 0.6
TIME = 4
node_num = 5

A_pos = [2]
B_pos = [3]
target = 0
ITERATION_STEP = 10000000

S_list = []
A_list = []
B_list = []

def Choosing(G, P_A, P_B, P_AB):
    p = np.random.rand()
    if p < P_A:
        G.nodes[i]["state_new"] = "A"
    elif P_A <= p < P_A + P_B:
        G.nodes[i]["state_new"] = "B"
    elif P_A + P_B <= p < P_A + P_B + P_AB:
         Choosing(G, P_A, P_B, P_AB)
    else:
        G.nodes[i]["state_new"] = "S"
    return G

for iter in range(1):
    S = []
    A = []
    B = []
    AB = []
    for Chosen_time in range(TIME):
        P_s = []
        GG = nx.random_tree(node_num, seed=1)
        G = GG.to_directed()

        for iteration_step in range(ITERATION_STEP):
            for i in range(node_num):
                G.nodes[i]["state"] = "S"
                G.nodes[i]["state_new"] = "S"
            # Initial Condition
            for Apos in A_pos:
                G.nodes[Apos]["state"] = "A"
                G.nodes[Apos]["state_new"] = "A"
            for Bpos in B_pos:
                G.nodes[Bpos]["state"] = "B"
                G.nodes[Bpos]["state_new"] = "B"

            for time in range(TIME):
                for i in range(node_num):
                    if G.nodes[i]["state"] == "S":
                        Neigh = []
                        for k in nx.all_neighbors(GG, i):
                            Neigh.append(G.nodes[k]["state"])
                        Num_A = Neigh.count("A")
                        Num_B = Neigh.count("B")
                        P_A = (1 - np.power((1 - p_A), Num_A)) * \
                            np.power((1 - p_B), Num_B)
                        P_B = (1 - np.power((1 - p_B), Num_B)) * \
                            np.power((1 - p_A), Num_A)
                        P_AB = (1 - np.power((1 - p_A), Num_A)) * \
                            (1 - np.power((1 - p_B), Num_B))
                        P_Ar = P_A / (1 - P_AB)
                        P_Br = P_B / (1 - P_AB)
                        
                        G = Choosing(G, P_A, P_B, P_AB)
                for i in range(node_num):
                    G.nodes[i]["state"] = G.nodes[i]["state_new"]
        
                for i in range(node_num):
                    if i == target and time == Chosen_time:
                        P_s.append(G.nodes[i]["state"])
        

        P1 = P_s.count("S") / len(P_s)
        P2 = P_s.count("A") / len(P_s)
        P3 = P_s.count("B") / len(P_s)

        S.append(P1)
        A.append(P2)
        B.append(P3)

    S_list.append(S)
    A_list.append(A)
    B_list.append(B)

print("S", S_list)
print("A", A_list)
print("B", B_list)
