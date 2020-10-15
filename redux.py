import os
import math
import copy

paper = {"Q": {0,1,2,3}
    , "Sigma": range(0,2)
    ,"delta": [[0,1],[2,1],[0,3],[2,3]]
    , "q0": 0
    ,"Gamma": {0,1}
    ,   "tau": [0,0,1,1]}
    
paper2 = {"Q": {0,1,2,3,4}
    , "Sigma": range(0,2)
    ,"delta": [[4,1],[2,1],[0,3],[2,3],[0,1]]
    , "q0": 0
    ,"Gamma": {0,1}
    ,   "tau": [0,0,1,1,0]}
    
thue = {"Q": {0,1}
    , "Sigma": range(0,2)
    , "delta": [[0,1],[1,0]]
    ,"Gamma": {0,1}
    , "tau": [0,1]}
    
odd0 = {"Q": {0,1}
    , "Sigma": range(0,2)
    , "delta": [[1,0],[0,1]]
    ,"Gamma": {0,1}
    , "tau": [0,1]}


def kDFAO(Gamma, k):
    return {"Q": set()
        , "Sigma" :range(0,k)
        ,"delta" :[]
        ,"q0" : 0
        ,"Gamma": Gamma
        , "tau" : []}

def minRegkDFAO(kDFAO, k):
    C = [kDFAO["delta"][kDFAO["q0"]][sigma] for sigma in range(1,k)]
    F = copy.deepcopy(C)
    while len(F) != 0 :
        q = F.pop()
        for sigma in range(0,k):
            p = kDFAO["delta"][q][sigma]
            if not p in C:
                C.append(p)
                F.append(p)
    return C
    
def hopcroftRegkDFAO(kDFAO, k):
    P = set()
    for gamma in kDFAO["Gamma"]:
        P.add(frozenset({q for q in kDFAO["Q"] if kDFAO["tau"][q] == gamma}))
    W = copy.deepcopy(P)
    while len(W) != 0 :
        A = W.pop();
        for sigma in range(0,k):
            X = {q for q in kDFAO["Q"] if kDFAO["delta"][q][sigma] in A}
            P2 = set()
            for Y in P:
                intersection = frozenset(X.intersection(Y))
                difference = frozenset(Y.difference(X))
                if len(intersection)!=0 and len(difference)!=0:
                    P2.add(intersection)
                    P2.add(difference)
                    if Y in W:
                        W.remove(Y)
                        W.add(intersection)
                        W.add(difference)
                    else:
                        W.add(intersection)
                        W.add(difference)
                else:
                    P2.add(Y)
            P = P2
    return P
    
def findInEquiv(state, equiv):
    counter = 0
    for eqClass in equiv:
        if state in eqClass:
            return counter
        counter += 1
    
def equivToDFAO(equiv,original,k):
    equiv = list(equiv) #convert equiv from type set to type list
    new = kDFAO({0,1},k);
    counter = 0;
    for eqClass in equiv:
        if original["q0"] in eqClass:
            new["q0"] = counter
        new["Q"].add(counter)
        representative = list(eqClass)[0]
        new["delta"].append([findInEquiv(original["delta"][representative][sigma],equiv) for sigma in range(0,k)])
        new["tau"].append(original["tau"][representative])
        counter += 1;
    return new 
        
def reduxRegkDFAO(kDFAO, k):
    reachable = minRegkDFAO(kDFAO,k)
    if kDFAO['q0'] not in reachable :
        print(len(reachable)+1)
        kDFAO["Q"] = set(reachable)
        equiv = hopcroftRegkDFAO(kDFAO,k)
        new = equivToDFAO(equiv,kDFAO,k)
        q0s = list()
        initials = [q for q in new["Q"] if kDFAO['tau'][kDFAO['q0'] == new['tau'][q]]]
        for sigma in range(1,k):
            counter = 0
            q0 = None
            for eqClass in equiv:
                if kDFAO["delta"][kDFAO["q0"]][sigma] in  eqClass:
                    q0 = counter
                    break
                counter += 1
            q0s.append(q0)
            initials = [q for q in initials if new["delta"][q][sigma] == q0]
        if len(initials) == 0:
            new["q0"] = len(new["Q"])
            new["Q"].add(new["q0"])
            q0s.insert(0,0)
            new["delta"].append(q0s)
            new["tau"].append(kDFAO['tau'][kDFAO['q0']])
        else:
            print(len(reachable))
            new["q0"] = list(initials)[0]
        return new
    else:
        reachable.append(0)
        kDFAO["Q"] = reachable
        equiv = hopcroftRegkDFAO(kDFAO,k)
        kDFAO = equivToDFAO(equiv,kDFAO,k)        
        return kDFAO
        
        
def fuse(b, n, a, k):
    new = kDFAO(a["Gamma"] | b["Gamma"],k)
    p = int(math.log(n+1,k))
    #P index = |b|*|a|*i +|a|*b + a 
    B = len(b["Q"])
    A = len(a["Q"])
    new["Q"] = set(range(0, A + p*A*B))
    new["q0"] = 0
    for i in range(0, p-1):
        for qb in b["Q"]:
            for qa in a["Q"]:
                new["delta"].append(
                    [A*B*(i+1)+A*(b["delta"][qb][sigma]) + a["delta"][qa][sigma]
                        for sigma in range(0,k)
                    ])
                new["tau"].append(b["tau"][qb])
    for qb in b["Q"]:
        for qa in a["Q"]:
            new["delta"].append(
                [B*A*p+a["delta"][qa][sigma]
                    for sigma in range(0,k)
                ])
            new["tau"].append(b["tau"][qb])
    for qa in a["Q"]:
        new["delta"].append(
            [A*B*p+a["delta"][qa][sigma] for sigma in range(0,k)]
            )
        new["tau"].append(a["tau"][qa])    
    return new
    
def convertDFAOToTabular(DFAO, k, start, end):
    alignments = "|"
    qs = "q "
    deltas = [f"$\\delta(q,{sigma})$" for sigma in range(0,k)]
    tau = "$\\tau(q)$"
    for q in range(start, end):
        alignments = "|c"+alignments
        qs += f"& {q}"
        for sigma in range(0,k):
            deltas[sigma] += f"& {DFAO['delta'][q][sigma]}"
        tau += f"& {DFAO['tau'][q]}"
    alignments = " c |" + alignments
    qs += "\\\\"
    tau += "\n"
    output = "\\begin{tabular}{"+alignments+"}\n"+qs
    for delta in deltas:
        output += delta + "\\\\ \n"
    output += tau + "\\end{tabular}"
    return output

def printDFAO(DFAO,k, maxTableSize = 13):
    size =  len(DFAO['Q'])
    start = 0
    output = ""
    while start+maxTableSize < size-1:
        output += convertDFAOToTabular(DFAO, k, start, start+maxTableSize) + "\n \n"
        start += maxTableSize
    output += convertDFAOToTabular(DFAO, k, start,size)
    file = open("texOutput.txt",'w')
    file.write(output)
    file.close()
for p in range(1,21):
    print(str(p) +":")
    fused = fuse(odd0,2**(p) - 1,thue,2)
    print(len(fused['Q']))
    reducedDFAO = reduxRegkDFAO(fused, 2)
    print(len(reducedDFAO['Q']))