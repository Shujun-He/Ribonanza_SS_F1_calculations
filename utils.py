import numpy as np

def detect_crossed_pairs(bp_list):
    """
    Detect crossed base pairs in a list of base pairs in RNA secondary structure.

    Args:
    bp_list (list of tuples): List of base pairs, where each tuple (i, j) represents a base pair.
    
    Returns:
    list of tuples: List of crossed base pairs.
    """
    crossed_pairs_set = set()
    crossed_pairs = []
    # Iterate through each pair of base pairs
    for i in range(len(bp_list)):
        for j in range(i+1, len(bp_list)):
            bp1 = bp_list[i]
            bp2 = bp_list[j]

            # Check if they are crossed
            if (bp1[0] < bp2[0] < bp1[1] < bp2[1]) or (bp2[0] < bp1[0] < bp2[1] < bp1[1]):
                crossed_pairs.append(bp1)
                crossed_pairs.append(bp2)
                crossed_pairs_set.add(bp1[0])
                crossed_pairs_set.add(bp1[1])
                crossed_pairs_set.add(bp2[0])
                crossed_pairs_set.add(bp2[1])
    return crossed_pairs, crossed_pairs_set

def dotbrackte2bp(structure):
    stack={'(':[],
           '[':[],
           '<':[],
           '{':[]}
    pop={')':'(',
         ']':'[',
         '>':"<",
         '}':'{'}       
    bp_list=[]
    matrix=np.zeros((len(structure),len(structure)))
    for i,s in enumerate(structure):
        if s in stack:
            stack[s].append((i,s))
        elif s in pop:
            forward_bracket=stack[pop[s]].pop()
            #bp_list.append(str(forward_bracket[0])+'-'+str(i))
            #bp_list.append([forward_bracket[0],i])
            bp_list.append([forward_bracket[0],i])

    return bp_list  


def calculate_f1_score_with_pseudoknots(true_pairs, predicted_pairs):
    true_pairs=[f"{i}-{j}" for i,j in true_pairs]
    predicted_pairs=[f"{i}-{j}" for i,j in predicted_pairs]
    
    true_pairs=set(true_pairs)
    predicted_pairs=set(predicted_pairs)

    # Calculate TP, FP, and FN
    TP = len(true_pairs.intersection(predicted_pairs))
    FP = len(predicted_pairs)-TP
    FN = len(true_pairs)-TP

    # Calculate Precision, Recall, and F1 Score
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1_score