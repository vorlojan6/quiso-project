
CC_3 = CyclotomicField(3) # number field
w = CC_3.gen() # third root of unity

def v(input_tup): # vector constructor for shorter notation
    return vector(CC_3, input_tup)

phi = [[None for _ in range(4)] for _ in range(10)] # init list of vectors defining the Penrose Dodecahedron  

# Define Penrose dodecahedron vectors (NOTE: ~w = w.conjugate())
##############################
phi[0][0] = v((1,  0,  0,  0))
phi[0][1] = v((0,  1,  0,  0))
phi[0][2] = v((0,  0,  1,  0))
phi[0][3] = v((0,  0,  0,  1))
phi[1][0] = v((0,  1, -1,  1))
phi[1][1] = v((1,  0, -1, -1))
phi[1][2] = v((1, -1,  0,  1))
phi[1][3] = v((1,  1,  1,  0))
phi[2][0] = v((0,  1, -w, ~w))
phi[2][1] = v((1,  0, -w,-~w))
phi[2][2] = v((1, -w,  0, ~w))
phi[2][3] = v((1,  w, ~w,  0))
phi[3][0] = v((0,  1,-~w,  w))
phi[3][1] = v((1,  0,-~w, -w))
phi[3][2] = v((1,-~w,  0,  w))
phi[3][3] = v((1, ~w,  w,  0))
phi[4][0] = v((0,  1, -w,  1))
phi[4][1] = v((1,  0, -1, -w))
phi[4][2] = v((1,-~w,  0, ~w))
phi[4][3] = v((1,  w,  1,  0))
phi[5][0] = v((0,  1,-~w, ~w))
phi[5][1] = v((1,  0, -w, -1))
phi[5][2] = v((1, -1,  0,  w))
phi[5][3] = v((1, ~w, ~w,  0))
phi[6][0] = v((0,  1, -1,  w))
phi[6][1] = v((1,  0,-~w,-~w))
phi[6][2] = v((1, -w,  0,  1))
phi[6][3] = v((1,  1,  w,  0))
phi[7][0] = v((0,  1,-~w,  1))
phi[7][1] = v((1,  0, -1,-~w))
phi[7][2] = v((1, -w,  0,  w))
phi[7][3] = v((1, ~w,  1,  0))
phi[8][0] = v((0,  1, -1, ~w))
phi[8][1] = v((1,  0, -w, -w))
phi[8][2] = v((1,-~w,  0,  1))
phi[8][3] = v((1,  1, ~w,  0))
phi[9][0] = v((0,  1, -w,  w))
phi[9][1] = v((1,  0,-~w, -1))
phi[9][2] = v((1, -1,  0, ~w))
phi[9][3] = v((1,  w,  w,  0))
##############################

phii = [phi[i][j] for i in range(10) for j in range(4)] # flatten the list (for easier iterating)
assert phii[1] == phi[0][1]

def is_basis(vecs) -> bool:
    return Matrix(vecs).rank() == 4

def dot(x, y): 
    return x.hermitian_inner_product(y)

def is_orthonormal_basis(vecs) -> bool:
    return is_basis(vecs) and all([dot(vecs[ii], vecs[jj]) == 0 for ii in range(4) for jj in range(ii + 1, 4)])

# print("Checking if each row forms a basis / orthonormal basis:")
# for i in range(10):
#     print(i, is_basis(phi[i]), is_orthonormal_basis(phi[i]))

n_Q, n_Qstar = 40, 40

#---------------------------------------------------------------------------------------#
# Analyse the orthogonality graph (Kochen-Specker-Diagram) of the Penrose Dodecahedron. #
#   - There is one node representing each vector                                        #
#   - Two nodes share an edge iff their corresponding vectors are orthogonal            #
#---------------------------------------------------------------------------------------#


# Compute neighbours list + parameter k_Q
####################################
neighbours_Q = [[] for _ in range(40)]
for i in range(40):
    for j in range(40):
        if i == j: continue
        if dot(phii[i], phii[j]) == 0:
            neighbours_Q[i].append(j)

assert len(set([len(ngh) for ngh in neighbours_Q])) == 1, "The graph Q is not strongly regular"
k_Q = list(set([len(ngh) for ngh in neighbours_Q]))[0]
####################################

# Compute parameter lambda_Q
####################################
lambda_list_Q = {}
for i in range(40):
    for j in neighbours_Q[i]:
        lambda_list_Q[(i, j)] = []
        for l in range(40):
            if l in (i, j): continue # dont count the neighbouring nodes themselves
            if l in neighbours_Q[i] and l in neighbours_Q[j]: lambda_list_Q[(i, j)].append(l)
assert len(set([len(val) for val in lambda_list_Q.values()])) == 1, "The graph Q is not strongly regular" # assert that lambda is the same for every node-pair
lambda_Q = list(set([len(val) for val in lambda_list_Q.values()]))[0]
####################################

# Compute parameter mu_Q
####################################
mu_list_Q = {}
for i in range(40):
    for j in [jj for jj in range(i + 1, 40) if jj not in neighbours_Q[i]]:
        mu_list_Q[(i, j)] = []
        for l in range(40):
            if l in neighbours_Q[i] and l in neighbours_Q[j]: mu_list_Q[(i, j)].append(l)

# print([len(val) for val in mu_list.values()])
assert len(set([len(val) for val in mu_list_Q.values()])) == 1, "The graph Q is not strongly regular" # assert that lambda is the same for every node-pair
mu_Q = list(set([len(val) for val in mu_list_Q.values()]))[0]
####################################
print(f"The graph Q  is strongly regular with parameters {n_Q, k_Q, lambda_Q, mu_Q}")
















#---------------------------------------------------------------------------------------#
# Analyse the graph Q* of tetrads related to the Penrose Dodecahedron                   #
#   - There is one node representing each orthogonal tetrad                             #
#   - Two nodes share an edge iff their corresponding tetrads share a vector            #
#---------------------------------------------------------------------------------------#

# Determine the nodes of Q*
####################################
tetrads = []
for i1 in range(40):
    for i2 in range(i1 + 1, 40):
        if dot(phii[i1], phii[i2]) != 0: continue
        for i3 in range(i2 + 1, 40):
            if (dot(phii[i1], phii[i3]) != 0) or (dot(phii[i2], phii[i3]) != 0): continue
            for i4 in range(i3 + 1, 40):
                if (dot(phii[i1], phii[i4]) != 0) or (dot(phii[i2], phii[i4]) != 0) or (dot(phii[i3], phii[i4]) != 0): continue
                tetrads.append((i1, i2, i3, i4))
# print(tetrads)
####################################


# Compute parameter k_Qstar
####################################
neighbours_Qstar = [[] for _ in range(40)]
for i in range(40):
    for j in range(40):
        if i == j: continue
        if len(set(tetrads[i] + tetrads[j])) < 8: # if the tetrads dont share any vector, their set-length should be 4 + 4 = 8
            neighbours_Qstar[i].append(j)
assert len(set([len(ngh) for ngh in neighbours_Qstar])) == 1, "The graph Q* is not strongly regular"
k_Qstar = list(set([len(ngh) for ngh in neighbours_Qstar]))[0]
####################################


# Compute parameter lambda_Qstar
####################################
lambda_list_Qstar = {}
for i in range(40):
    for j in neighbours_Qstar[i]:
        lambda_list_Qstar[(i, j)] = []
        for l in range(40):
            if l in (i, j): continue # dont count the neighbouring nodes themselves
            if l in neighbours_Qstar[i] and l in neighbours_Qstar[j]: lambda_list_Qstar[(i, j)].append(l)
assert len(set([len(val) for val in lambda_list_Qstar.values()])) == 1, "The graph Q* is not strongly regular" # assert that lambda is the same for every node-pair
lambda_Qstar = list(set([len(val) for val in lambda_list_Qstar.values()]))[0]
####################################

# Compute parameter mu_Qstar
####################################
mu_list_Qstar = {}
for i in range(40):
    for j in [jj for jj in range(i + 1, 40) if jj not in neighbours_Qstar[i]]:
        mu_list_Qstar[(i, j)] = []
        for l in range(40):
            if l in neighbours_Qstar[i] and l in neighbours_Qstar[j]: mu_list_Qstar[(i, j)].append(l)

assert len(set([len(val) for val in mu_list_Qstar.values()])) == 1, "The graph Q* is not strongly regular" # assert that lambda is the same for every node-pair
mu_Qstar = list(set([len(val) for val in mu_list_Qstar.values()]))[0]
####################################
print(f"The graph Q* is strongly regular with parameters {n_Qstar, k_Qstar, lambda_Qstar, mu_Qstar}")