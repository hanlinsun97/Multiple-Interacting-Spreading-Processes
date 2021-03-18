import networkx as nx
import numpy as np


def initialization(G, T):
    for edge in G.edges():
        G.edges[edge]["theta_AB_simple"] = np.full_like(np.empty([T+1, T+1]),np.nan)
        G.edges[edge]["theta_AB"] = np.full_like(np.empty([T+1]),np.nan)
        G.edges[edge]["theta_BA"] = np.full_like(np.empty([T+1]),np.nan)
        G.edges[edge]["theta_ABAB"] = np.full_like(np.empty([T+1, T+1, T+1, T+1]),np.nan)
        G.edges[edge]["theta_BAAB"] = np.full_like(np.empty([T+1, T+1, T+1, T+1]),np.nan)
        G.edges[edge]["theta_A"] = np.full_like(np.empty([T+1, T+1]),np.nan)
        G.edges[edge]["theta_B"] = np.full_like(np.empty([T+1, T+1]),np.nan)

        G.edges[edge]["phi_AB_simple"] = np.full_like(np.empty([T+1, T+1]),np.nan)
        G.edges[edge]["phi_AB"] = np.full_like(np.empty([T+1]),np.nan)
        G.edges[edge]["phi_BA"] = np.full_like(np.empty([T+1]),np.nan)
        G.edges[edge]["phi_ABAB"] = np.full_like(np.empty([T+1, T+1, T+1, T+1]),np.nan)
        G.edges[edge]["phi_BAAB"] = np.full_like(np.empty([T+1, T+1, T+1, T+1]),np.nan)
        G.edges[edge]["phi_A"] = np.full_like(np.empty([T+1,T+1]),np.nan)
        G.edges[edge]["phi_B"] = np.full_like(np.empty([T+1,T+1]),np.nan)

        G.edges[edge]["message"] = np.full_like(np.empty([T+1, T+1, T+1]),np.nan)
        G.edges[edge]["message_star_A"] = np.full_like(np.empty([T+1, T+1]),np.nan)
        G.edges[edge]["message_star_B"] = np.full_like(np.empty([T+1, T+1]),np.nan)

        G.edges[edge]["mu_A"] = np.full_like(np.empty([T+1]),np.nan)
        G.edges[edge]["mu_B"] = np.full_like(np.empty([T+1]),np.nan)
    
    for node in G.nodes():
        G.nodes[node]["message"]=np.full_like(np.empty([T+1, T+1, T+1]),np.nan)
        G.nodes[node]["message_star_A"] = np.full_like(np.empty([T+1, T+1]),np.nan)
        G.nodes[node]["message_star_B"] = np.full_like(np.empty([T+1, T+1]),np.nan)

        G.nodes[node]["mu_A"] = np.full_like(np.empty([T+1]),np.nan)
        G.nodes[node]["mu_B"] = np.full_like(np.empty([T+1]),np.nan)
        G.nodes[node]["A"] = np.full_like(np.empty([T+1]),np.nan)
        G.nodes[node]["B"] = np.full_like(np.empty([T+1]),np.nan)
        G.nodes[node]["S"]=np.full_like(np.empty([T+1]),np.nan)
        G.nodes[node]["AB"] = np.full_like(np.empty([T+1]),np.nan)
        G.nodes[node]["A_star"] = np.full_like(np.empty([T+1]), np.nan)
        G.nodes[node]["B_star"]=np.full_like(np.empty([T+1]),np.nan)

    return G

def first_values(G, alpha_A, alpha_B, alpha_AB, alpha_BA, A_pos, B_pos, AB_pos):
    for node in G.nodes():
        if node in A_pos:
            G.nodes[node]["S"][0] = 0
            G.nodes[node]["A_star"][0] = 1
            G.nodes[node]["A"][0] = 1
            G.nodes[node]["B_star"][0] = 0
            G.nodes[node]["B"][0] = 0
            G.nodes[node]["AB"][0] = 0
        elif node in B_pos:
            G.nodes[node]["S"][0] = 0
            G.nodes[node]["A_star"][0] = 0
            G.nodes[node]["A"][0] = 0
            G.nodes[node]["B_star"][0] = 1
            G.nodes[node]["B"][0] = 1
            G.nodes[node]["AB"][0] = 0
        elif node in AB_pos:
            G.nodes[node]["S"][0] = 0
            G.nodes[node]["A_star"][0] = 0
            G.nodes[node]["A"][0] = 1
            G.nodes[node]["B_star"][0] = 0
            G.nodes[node]["B"][0] = 1
            G.nodes[node]["AB"][0] = 1
        else:
            G.nodes[node]["S"][0] = 1
            G.nodes[node]["A_star"][0] = 0
            G.nodes[node]["A"][0] = 0
            G.nodes[node]["B_star"][0] = 0
            G.nodes[node]["B"][0] = 0
            G.nodes[node]["AB"][0] = 0

    for node in G.nodes():
        G.nodes[node]["message"][0,0,0] = G.nodes[node]["AB"][0]
        G.nodes[node]["message_star_B"][0,0] = G.nodes[node]["A_star"][0]
        G.nodes[node]["message_star_A"][0,0] = G.nodes[node]["B_star"][0]
        G.nodes[node]["mu_A"][0] = G.nodes[node]["A"][0]
        G.nodes[node]["mu_B"][0] = G.nodes[node]["B"][0]

    for edge in G.edges():
        i = edge[0]
        j = edge[1]
        G.edges[edge]["message"][0,0,0] = G.nodes[i]["AB"][0]
        G.edges[edge]["message_star_B"][0,0] = G.nodes[i]["A_star"][0]
        G.edges[edge]["message_star_A"][0,0] = G.nodes[i]["B_star"][0]

        G.edges[edge]["theta_AB_simple"][0,0] = 1
        G.edges[edge]["theta_AB"][0] = 1
        G.edges[edge]["theta_BA"][0] = 1
        
        G.edges[edge]["phi_A"][0,0] = G.nodes[i]["A"][0]
        G.edges[edge]["phi_B"][0,0] = G.nodes[i]["B"][0]
        G.edges[edge]["phi_AB_simple"][0,0] = G.nodes[i]["AB"][0]
        G.edges[edge]["phi_AB"][0] = G.nodes[i]["B"][0]
        G.edges[edge]["phi_BA"][0] = G.nodes[i]["A"][0]
    return G

def marginals(G, t):
    for node in G.nodes():
        prod_S = 1
        G.nodes[node]["A"][t] = 0
        G.nodes[node]["B"][t] = 0
        for t_prime in range(t+1):
            G.nodes[node]["A"][t] += G.nodes[node]["mu_A"][t_prime]
            G.nodes[node]["B"][t] += G.nodes[node]["mu_B"][t_prime]
        for k in nx.all_neighbors(G.to_undirected(), node):
            prod_S = prod_S * G.edges[k, node]["theta_AB_simple"][t,t]
        
        G.nodes[node]["S"][t] = G.nodes[node]["S"][0] * prod_S
        G.nodes[node]["AB"][t] = G.nodes[node]["A"][t] + G.nodes[node]["B"][t] + G.nodes[node]["S"][t] - 1
        G.nodes[node]["A_star"][t] = G.nodes[node]["A"][t] - G.nodes[node]["AB"][t]
        G.nodes[node]["B_star"][t] = G.nodes[node]["B"][t] - G.nodes[node]["AB"][t]
    return G
    
def message_edge(G, edge, tau_A, tau_B,t):
    i = edge[0]
    j = edge[1]
    if tau_A == 0 and tau_B !=0:
        prod_1 = 1
        prod_2 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            if k != j:
                prod_1 = prod_1 * G.edges[k, i]["theta_AB"][tau_B-1]
                prod_2 = prod_2 * G.edges[k, i]["theta_AB"][tau_B]
        G.edges[edge]["message"][tau_A, tau_B, t] = G.nodes[i]["A_star"][0] * (prod_1 - prod_2)
    elif tau_B == 0 and tau_A != 0:
        prod_1 = 1
        prod_2 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            if k != j:
                prod_1 = prod_1 * G.edges[k, i]["theta_BA"][tau_A-1]
                prod_2 = prod_2 * G.edges[k, i]["theta_BA"][tau_A]
        G.edges[edge]["message"][tau_A, tau_B, t] = G.nodes[i]["B_star"][0] * (prod_1 - prod_2)
    elif tau_A > 0 and tau_B > 0 and tau_B > tau_A:
        prod_1 = 1
        prod_2 = 1
        prod_3 = 1
        prod_4 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            if k != j:
                prod_1 = prod_1 * G.edges[k, i]["theta_ABAB"][tau_A-1, tau_A, tau_B-1, tau_A]
                prod_2 = prod_2 * G.edges[k, i]["theta_ABAB"][tau_A-1, tau_A, tau_B, tau_A]
                prod_3 = prod_3 * G.edges[k, i]["theta_ABAB"][tau_A, tau_A, tau_B-1, tau_A]
                prod_4 = prod_4 * G.edges[k, i]["theta_ABAB"][tau_A, tau_A, tau_B, tau_A]
        G.edges[edge]["message"][tau_A, tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2 - prod_3 + prod_4)
    elif tau_A > 0 and tau_B > 0 and tau_B < tau_A:
        prod_1 = 1
        prod_2 = 1
        prod_3 = 1
        prod_4 = 1
        
        for k in nx.all_neighbors(G.to_undirected(), i):
            if k != j:
                prod_1 = prod_1 * G.edges[k, i]["theta_BAAB"][tau_B-1, tau_B, tau_A-1, tau_B]
                prod_2 = prod_2 * G.edges[k, i]["theta_BAAB"][tau_B-1, tau_B, tau_A, tau_B]
                prod_3 = prod_3 * G.edges[k, i]["theta_BAAB"][tau_B, tau_B, tau_A-1, tau_B]
                prod_4 = prod_4 * G.edges[k, i]["theta_BAAB"][tau_B, tau_B, tau_A, tau_B]
        G.edges[edge]["message"][tau_A, tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2 - prod_3 + prod_4)   
    elif tau_A > 0 and tau_B > 0 and tau_A == tau_B:
        prod_1 = 1
        prod_2 = 1
        prod_3 = 1
        prod_4 = 1
        
        for k in nx.all_neighbors(G.to_undirected(), i):
            if k != j:
                prod_1 = prod_1 * G.edges[k, i]["theta_AB_simple"][tau_A-1, tau_A-1]
                prod_2 = prod_2 * G.edges[k, i]["theta_AB_simple"][tau_A-1, tau_A]
                prod_3 = prod_3 * G.edges[k, i]["theta_AB_simple"][tau_A, tau_A-1]
                prod_4 = prod_4 * G.edges[k, i]["theta_AB_simple"][tau_A, tau_A]

        G.edges[edge]["message"][tau_A, tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2 - prod_3 + prod_4)
    return G

def message_edge_star(G, edge, tau_A, tau_B,t):
    i = edge[0]
    j = edge[1]
    if tau_B == "*" and tau_A != "*":
        if tau_A == 0:
            prod = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                if k != j:
                    prod = prod * G.edges[k, i]["theta_AB"][t]
            G.edges[edge]["message_star_B"][tau_A, t] = G.nodes[i]["A_star"][0] * prod
        elif tau_A > 0:
            prod_1 = 1
            prod_2 = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                if k != j:
                    prod_1 = prod_1 * G.edges[k, i]["theta_ABAB"][tau_A-1, tau_A, t, tau_A]
                    prod_2 = prod_2 * G.edges[k, i]["theta_ABAB"][tau_A, tau_A, t, tau_A]
            G.edges[edge]["message_star_B"][tau_A, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2)
       
    elif tau_A == "*" and tau_B != "*":
        if tau_B == 0:
            prod = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                if k != j:
                    prod = prod * G.edges[k, i]["theta_BA"][t]
            G.edges[edge]["message_star_A"][tau_B, t] = G.nodes[i]["B_star"][0] * prod
        elif tau_B > 0:
            prod_1 = 1
            prod_2 = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                if k != j:
                    prod_1 = prod_1 * G.edges[k, i]["theta_BAAB"][tau_B-1, tau_B, t, tau_B]
                    prod_2 = prod_2 * G.edges[k, i]["theta_BAAB"][tau_B, tau_B, t, tau_B]
            G.edges[edge]["message_star_A"][tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2)
    elif tau_A == "*" and tau_B == "*":
        prod = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            if k != j:
                prod = prod * G.edges[k, i]["theta_AB_simple"][t,t]
        G.edges[edge]["message_star_A_B"][t] = G.nodes[i]["S"][0] * prod
    return G

def message_node(G, node, tau_A, tau_B, t):
    
    i = node
    if tau_A == 0 and tau_B !=0:
        prod_1 = 1
        prod_2 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            prod_1 = prod_1 * G.edges[k, i]["theta_AB"][tau_B-1]
            prod_2 = prod_2 * G.edges[k, i]["theta_AB"][tau_B]
        G.nodes[node]["message"][tau_A, tau_B, t] = G.nodes[i]["A_star"][0] * (prod_1 - prod_2)
    elif tau_B == 0 and tau_A != 0:
        prod_1 = 1
        prod_2 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            prod_1 = prod_1 * G.edges[k, i]["theta_BA"][tau_A-1]
            prod_2 = prod_2 * G.edges[k, i]["theta_BA"][tau_A]
        G.nodes[node]["message"][tau_A, tau_B, t] = G.nodes[i]["B_star"][0] * (prod_1 - prod_2)
        
    elif tau_A > 0 and tau_B > 0 and tau_B > tau_A:
        prod_1 = 1
        prod_2 = 1
        prod_3 = 1
        prod_4 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            prod_1 = prod_1 * G.edges[k, i]["theta_ABAB"][tau_A-1, tau_A, tau_B-1, tau_A]
            prod_2 = prod_2 * G.edges[k, i]["theta_ABAB"][tau_A-1, tau_A, tau_B, tau_A]
            prod_3 = prod_3 * G.edges[k, i]["theta_ABAB"][tau_A, tau_A, tau_B-1, tau_A]
            prod_4 = prod_4 * G.edges[k, i]["theta_ABAB"][tau_A, tau_A, tau_B, tau_A]
        G.nodes[node]["message"][tau_A, tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2 - prod_3 + prod_4)
    elif tau_A > 0 and tau_B > 0 and tau_B < tau_A:
        prod_1 = 1
        prod_2 = 1
        prod_3 = 1
        prod_4 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            prod_1 = prod_1 * G.edges[k, i]["theta_BAAB"][tau_B-1, tau_B, tau_A-1, tau_B]
            prod_2 = prod_2 * G.edges[k, i]["theta_BAAB"][tau_B-1, tau_B, tau_A, tau_B]
            prod_3 = prod_3 * G.edges[k, i]["theta_BAAB"][tau_B, tau_B, tau_A-1, tau_B]
            prod_4 = prod_4 * G.edges[k, i]["theta_BAAB"][tau_B, tau_B, tau_A, tau_B]
        G.nodes[node]["message"][tau_A, tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2 - prod_3 + prod_4)
    elif tau_A > 0 and tau_B > 0 and tau_A == tau_B:
        prod_1 = 1
        prod_2 = 1
        prod_3 = 1
        prod_4 = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            prod_1 = prod_1 * G.edges[k, i]["theta_AB_simple"][tau_A-1, tau_A-1]
            prod_2 = prod_2 * G.edges[k, i]["theta_AB_simple"][tau_A-1, tau_A]
            prod_3 = prod_3 * G.edges[k, i]["theta_AB_simple"][tau_A, tau_A - 1]
            prod_4 = prod_4 * G.edges[k, i]["theta_AB_simple"][tau_A, tau_A]
        G.nodes[node]["message"][tau_A, tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2 - prod_3 + prod_4)
    return G

def message_node_star(G, node, tau_A, tau_B, t):
    i = node
    if tau_B == "*" and tau_A != "*":
        if tau_A == 0:
            prod = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                prod = prod * G.edges[k, i]["theta_AB"][t]
            G.nodes[node]["message_star_B"][tau_A, t] = G.nodes[i]["A_star"][0] * prod
        elif tau_A > 0:
            prod_1 = 1
            prod_2 = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                prod_1 = prod_1 * G.edges[k, i]["theta_ABAB"][tau_A-1, tau_A, t, tau_A]
                prod_2 = prod_2 * G.edges[k, i]["theta_ABAB"][tau_A, tau_A, t, tau_A]
            G.nodes[node]["message_star_B"][tau_A, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2)
    elif tau_A == "*" and tau_B != "*":
        if tau_B == 0:
            prod = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                prod = prod * G.edges[k, i]["theta_BA"][t]
            G.nodes[node]["message_star_A"][tau_B, t] = G.nodes[i]["B_star"][0] * prod
        elif tau_B > 0:
            prod_1 = 1
            prod_2 = 1
            for k in nx.all_neighbors(G.to_undirected(), i):
                prod_1 = prod_1 * G.edges[k, i]["theta_BAAB"][tau_B-1, tau_B, t, tau_B]
                prod_2 = prod_2 * G.edges[k, i]["theta_BAAB"][tau_B, tau_B, t, tau_B]
            G.nodes[node]["message_star_A"][tau_B, t] = G.nodes[i]["S"][0] * (prod_1 - prod_2)
    elif tau_A == "*" and tau_B == "*":
        prod = 1
        for k in nx.all_neighbors(G.to_undirected(), i):
            prod = prod * G.edges[k, i]["theta_AB_simple"][t,t]
        G.nodes[node]["message_star_A_B"][t] = G.nodes[i]["S"][0] * prod
    return G

def mu_edge(G, edge, t):
    sum_1 = 0
    for t_prime in range(t+1):
        sum_1 += G.edges[edge]["message"][t, t_prime, t]
    G.edges[edge]["mu_A"][t] = sum_1 + G.edges[edge]["message_star_B"][t, t]

    # G.edges[edge]["theta_A"][0,t] = G.edges[edge]["mu_A"][t]

    sum_2 = 0
    for t_prime in range(t+1):
        sum_2 += G.edges[edge]["message"][t_prime, t, t]
    G.edges[edge]["mu_B"][t] = sum_2 + G.edges[edge]["message_star_A"][t, t]

    # G.edges[edge]["theta_B"][t,0] = G.edges[edge]["mu_B"][t]

    return G

def mu_node(G, node, t):
    sum_1 = 0
    for t_prime in range(t+1):
        sum_1 += G.nodes[node]["message"][t, t_prime, t]
    G.nodes[node]["mu_A"][t] = sum_1 + G.nodes[node]["message_star_B"][t, t]
    
    sum_2 = 0
    for t_prime in range(t+1):
        sum_2 += G.nodes[node]["message"][t_prime, t, t]
    G.nodes[node]["mu_B"][t] = sum_2 + G.nodes[node]["message_star_A"][t, t]
    return G

def theta_AB_simple(G, edge, t):
    G.edges[edge]["theta_AB_simple"][t+1,t+1] = G.edges[edge]["theta_AB_simple"][t,t]- alpha_A * G.edges[edge]["phi_A"][t,t] - alpha_B * G.edges[edge]["phi_B"][t,t] + alpha_A * alpha_B * G.edges[edge]["phi_AB_simple"][t,t]
    G.edges[edge]["theta_AB_simple"][t+1, t] = G.edges[edge]["theta_AB_simple"][t,t] - alpha_A * G.edges[edge]["phi_A"][t,t]
    G.edges[edge]["theta_AB_simple"][t, t+1] = G.edges[edge]["theta_AB_simple"][t,t] - alpha_B * G.edges[edge]["phi_B"][t,t]
   
    for t_AB in range(T+1):
        for tau_A in range(T+1):
            if tau_A >= t_AB:
                G.edges[edge]["theta_ABAB"][t, t+1, t_AB, tau_A] = G.edges[edge]["theta_AB_simple"][t, t+1]
                G.edges[edge]["theta_ABAB"][t+1, t, t_AB, tau_A] = G.edges[edge]["theta_AB_simple"][t+1, t]
                G.edges[edge]["theta_ABAB"][t+1, t+1, t_AB, tau_A]= G.edges[edge]["theta_AB_simple"][t+1, t+1]
    for t_AB in range(T+1):
        for tau_B in range(T+1):
            if tau_B >= t_AB:
                G.edges[edge]["theta_BAAB"][t+1, t, t_AB, tau_B]= G.edges[edge]["theta_AB_simple"][t, t+1]  #BUG 
                G.edges[edge]["theta_BAAB"][t, t+1, t_AB, tau_B]= G.edges[edge]["theta_AB_simple"][t+1, t]
                G.edges[edge]["theta_BAAB"][t+1, t+1, t_AB, tau_B]= G.edges[edge]["theta_AB_simple"][t+1, t+1]
    return G

def ind(tau, t_prime):
    if t_prime >= tau:
        indd = 1
    else:
        indd = 0
    return indd

def phi_AB_simple(G, edge, t):
    sum_1 = 0
    sum_2 = 0
    
    for tau_kA in range(t+1):
        prod_1 = 1
        for t_prime in range(t):
            prod_1 = prod_1 * (1 - alpha_A * ind(tau_kA, t_prime)) 
        sum_1 = sum_1 + prod_1 * G.edges[edge]["message"][tau_kA, t, t]
    for tau_kB in range(t+1):
        prod_2 = 1
        for t_prime in range(t):
            prod_2 = prod_2 * (1 - alpha_B * ind(tau_kB, t_prime)) 
        sum_2 = sum_2 + prod_2 * G.edges[edge]["message"][t, tau_kB, t]
    G.edges[edge]["phi_AB_simple"][t,t] = (1 - alpha_A) * (1 - alpha_B)* G.edges[edge]["phi_AB_simple"][t-1,t-1] - G.edges[edge]["message"][t,t,t] + sum_1 + sum_2
    
    return G

def phi(G, edge, t):
    
    G.edges[edge]["phi_A"][t,t] = (1 - alpha_A) * G.edges[edge]["phi_A"][t-1,t-1] - alpha_B * (1 - alpha_A) * G.edges[edge]["phi_AB_simple"][t-1,t-1] + G.edges[edge]["theta_B"][t,t]
    G.edges[edge]["phi_B"][t,t] = (1 - alpha_B) * G.edges[edge]["phi_B"][t-1,t-1] - alpha_A * (1 - alpha_B) * G.edges[edge]["phi_AB_simple"][t-1,t-1] + G.edges[edge]["theta_A"][t,t]
    
    G.edges[edge]["phi_A"][t,t-1] = (1 - alpha_A) * G.edges[edge]["phi_A"][t-1, t-1] + G.edges[edge]["theta_B"][t, t-1]
    G.edges[edge]["phi_B"][t-1,t] = (1 - alpha_B) * G.edges[edge]["phi_B"][t-1, t-1] + G.edges[edge]["theta_A"][t-1, t]
    
    
    for t_AB in range(T+1):
        for tau_A in range(T+1):
            if tau_A >= t_AB and t_AB == t:
                G.edges[edge]["phi_ABAB"][t, t, t_AB, tau_A] = G.edges[edge]["phi_B"][t, t]
                G.edges[edge]["phi_ABAB"][t-1, t, t_AB, tau_A] = G.edges[edge]["phi_B"][t-1, t]
        for tau_B in range(T+1):
            if tau_B >= t_AB and t_AB == t:
                G.edges[edge]["phi_BAAB"][t,t,t_AB, tau_B] = G.edges[edge]["phi_A"][t,t]
                G.edges[edge]["phi_BAAB"][t-1,t,t_AB, tau_B] = G.edges[edge]["phi_A"][t,t-1]
    return G

def theta(G, edge, t):

    G.edges[edge]["theta_A"][0, t] = G.edges[edge]["mu_B"][t]
    G.edges[edge]["theta_B"][t, 0] = G.edges[edge]["mu_A"][t]
    for tt in range(1, t+1):
        sum_A = 0
        for tau_A in range(tt):
            prod_A = 1
            for t_prime in range(tt-1):
                value_A = ind(tau_A, t_prime)
                prod_A = prod_A * (1 - alpha_A * value_A)
            sum_A += prod_A * G.edges[edge]["message"][tau_A, t,t]
        
        G.edges[edge]["theta_A"][tt, t] = G.edges[edge]["theta_A"][tt-1, t] - alpha_A * sum_A

        sum_B = 0
        for tau_B in range(tt):
            prod_B = 1
            for t_prime in range(tt-1):
                value_B = ind(tau_B, t_prime)
                prod_B = prod_B * (1 - alpha_B * value_B)
            sum_B += prod_B * G.edges[edge]["message"][t, tau_B,t]
        G.edges[edge]["theta_B"][t, tt]= G.edges[edge]["theta_B"][t, tt-1] - alpha_B * sum_B
        # if t == 1 and edge==(2,0):
        #     print("A",G.edges[edge]["theta_B"][1,1])
            
    return G

def theta_AB(G, edge, t):
    G.edges[edge]["theta_AB"][t+1] = G.edges[edge]["theta_AB"][t] - alpha_BA * G.edges[edge]["phi_AB"][t]
    G.edges[edge]["theta_BA"][t+1] = G.edges[edge]["theta_BA"][t] - alpha_AB * G.edges[edge]["phi_BA"][t]
    return G

def phi_AB(G, edge, t):
    G.edges[edge]["phi_AB"][t] = (1 - alpha_BA) * G.edges[edge]["phi_AB"][t-1] + G.edges[edge]["mu_B"][t]
    G.edges[edge]["phi_BA"][t] = (1 - alpha_AB) * G.edges[edge]["phi_BA"][t-1] + G.edges[edge]["mu_A"][t]
    return G

def theta_ABAB(G, edge, t_A, t_B, t_AB, tau_A):
    if tau_A <= t_AB:
        indd = 1
    else:
        indd = 0
    G.edges[edge]["theta_ABAB"][t_A, t_B, t_AB+1, tau_A] = G.edges[edge]["theta_ABAB"][t_A, t_B, t_AB, tau_A] - indd * alpha_BA * G.edges[edge]["phi_ABAB"][t_A, t_B, t_AB, tau_A]
    return G

def theta_BAAB(G, edge, t_B, t_A, t_AB, tau_B):
    if tau_B <= t_AB:
        indd = 1
    else:
        indd = 0  # BUG
    G.edges[edge]["theta_BAAB"][t_B, t_A, t_AB+1, tau_B] = G.edges[edge]["theta_BAAB"][t_B, t_A, t_AB, tau_B] - indd * alpha_AB * G.edges[edge]["phi_BAAB"][t_B, t_A, t_AB, tau_B]
    return G

def phi_ABAB(G, edge, t_A, t_B, t_AB, tau_A):
    sum_1 = 0
    for tau_kA in range(t_AB+1):
        prod_1 = 1
        for t_prime in range(t_A):
            prod_1 = prod_1 * (1 - alpha_A * ind(tau_kA, t_prime)) 
        sum_1 = sum_1 + prod_1 * G.edges[edge]["message"][tau_kA, t_AB, t_AB]
    G.edges[edge]["phi_ABAB"][t_A, t_B, t_AB, tau_A] = (1 - alpha_BA) * G.edges[edge]["phi_ABAB"][t_A, t_B, t_AB-1, tau_A] + sum_1 + G.edges[edge]["message_star_A"][t_AB, t_AB]
    return G

def phi_BAAB(G, edge, t_B, t_A, t_AB, tau_B):
    sum_2 = 0
    for tau_kB in range(t_AB+1):
        prod_2 = 1
        for t_prime in range(t_B):
            prod_2 = prod_2 * (1 - alpha_B * ind(tau_kB, t_prime)) 
        sum_2 = sum_2 + prod_2 * G.edges[edge]["message"][t_AB, tau_kB, t_AB]
    G.edges[edge]["phi_BAAB"][t_B, t_A, t_AB, tau_B] = (1 - alpha_AB) * G.edges[edge]["phi_BAAB"][t_B, t_A, t_AB-1, tau_B] + sum_2 + G.edges[edge]["message_star_B"][t_AB, t_AB]
    return G


T = 7
alpha_A = 0.6
alpha_B = 0.5
alpha_AB = 0.2
alpha_BA = 0.3
A_pos = [2]
B_pos = [0]
AB_pos = []
target = 3

G = nx.random_tree(6, seed=2) 
# G = nx.star_graph(2)

G = G.to_directed()

G = initialization(G, T)

for t in range(T):
    if t == 0:
        G = first_values(G, alpha_A, alpha_B, alpha_AB, alpha_BA, A_pos, B_pos, AB_pos)
    elif t > 0:
        for edge in G.edges():
            G = theta_AB_simple(G, edge, t-1)
        for edge in G.edges():
            G = theta_AB(G, edge, t-1)
        for edge in G.edges():
            G = message_edge_star(G, edge, t, "*", t)
            G = message_edge_star(G, edge, "*", t, t)
        for node in G.nodes():
            G = message_node_star(G, node, t, "*", t)
            G = message_node_star(G, node, "*", t, t)
        for t_prime in range(1, t):
            for edge in G.edges():  
                G = theta_ABAB(G, edge, t_prime, t_prime, t-1, t_prime)
                G = theta_ABAB(G, edge, t_prime - 1, t_prime, t-1, t_prime)
                G = theta_BAAB(G, edge, t_prime, t_prime, t-1, t_prime)
                G = theta_BAAB(G, edge, t_prime - 1, t_prime, t-1, t_prime)
        for t_prime in range(t+1):
            for edge in G.edges():
                G = message_edge(G, edge, t, t_prime, t)
                G = message_edge(G, edge, t_prime, t, t)
            for node in G.nodes():
                G = message_node(G, node, t, t_prime, t)
                G = message_node(G, node, t_prime, t, t)
        for edge in G.edges():
            G = mu_edge(G, edge, t)
        for node in G.nodes():
            G = mu_node(G, node, t)
        for edge in G.edges():
            G = theta(G, edge, t)
        for t_prime in range(1, t):
            for edge in G.edges():
                G = phi_ABAB(G, edge, t_prime, t_prime, t, t_prime)
                G = phi_ABAB(G, edge, t_prime-1, t_prime, t, t_prime)
                G = phi_BAAB(G, edge, t_prime, t_prime, t, t_prime)
                G = phi_BAAB(G, edge, t_prime-1, t_prime, t, t_prime)
        for edge in G.edges():
            G = phi(G, edge, t)
            G = phi_AB_simple(G, edge, t)
            G = phi_AB(G, edge, t)
        G = marginals(G, t)


print("S", G.nodes[target]["S"])
print("A", G.nodes[target]["A_star"])
print("B", G.nodes[target]["B_star"])
print("AB", G.nodes[target]["AB"])
    
