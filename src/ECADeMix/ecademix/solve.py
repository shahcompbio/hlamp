import numpy as np
from sklearn.decomposition import NMF
import gurobipy as gp
from gurobipy import GRB
from datetime import datetime
import scipy.special

# set up and run ILP
def run_NMF(data, n_amplicons, seed=0):
    """
    Run NMF on the given data and return the W and H matrices.
    """
    nmf = NMF(n_components=n_amplicons, init='random', random_state=seed)
    W = nmf.fit_transform(data)
    H = nmf.components_
    return W, H
    
def construct_solution(W, H, max_amplicon_cn, verbose=True):
    """
    Ensure that the NMF solution is a valid solution to the ILP.
    """
    loading_start = W.copy()
    amplicon_start = np.round(H)

    # transfer weight from amplicons to loadings to ensure max_amplicon_cn is observed
    for p in range(amplicon_start.shape[0]):
        if np.max(amplicon_start[p]) > max_amplicon_cn:
            factor = np.max(amplicon_start[p]) // max_amplicon_cn + 1
            if verbose:
                print(f'transfering factor of {factor} from amplicon {p} to loadings for amplicon {p}')
            loading_start[:, p] = W[:, p] * factor
            amplicon_start[p] = np.round(H[p] / factor)

    # make sure each row and col of loadings is nonzero
    for i in range(loading_start.shape[0]):
        if np.sum(loading_start[i]) == 0:
            if verbose:
                print(f"supplementing W[{i}, {np.argmax(W[i])}] = 1 (zero row)")
            loading_start[i, np.argmax(W[i])] = 1
    for p in range(loading_start.shape[1]):
        if np.sum(loading_start[:, p]) == 0:
            loading_start[np.argmax(W[:, p]), p] = 1
            if verbose:
                print(f"supplementing W[{np.argmax(W[:, p])}, {p}] = 1 (zero col)")
    amplicon_start = amplicon_start.astype(int)
    return loading_start, amplicon_start

def compute_main_objective(data, loadings, amplicons):
    n_amplicons = amplicons.shape[0]
    return sum(
            np.abs(data[i, j] - sum(loadings[i, p] * amplicons[p, j] 
            for p in range(n_amplicons)))
            for i in range(data.shape[0]) 
            for j in range(data.shape[1])
            )


def get_model(data, n_amplicons, 
              min_loading, max_loading, max_amplicon_cn, 
              loading_start=None, amplicon_start=None,
              loadings=None, amplicons=None, 
              lambda_x=1, lambda_y=1, lambda_x_l1=0,
              n_threads=1, time_per_run=60, cell_weights=False, lambda_dot=0):
    
    if loadings is None and amplicons is None:
        raise ValueError("Loadings or amplicons must be specified as input.")
    elif loadings is None:
        infer_loadings = True
        assert loading_start is not None, "Initial [loading_start] must be provided for the [loadings] to be inferred."
    else:
        infer_loadings = False
        assert amplicon_start is not None, "Initial [amplicon_start] values must be provided for the [amplicons] to be inferred."

    data_matrix = data.copy()
    
    m, n = data_matrix.shape  # cells, bins
    model = gp.Model()

    ## Variables
    # X = amplicons (n_amplicons x bin)
    # Y = amplicon loadings (cell x n_amplicons)
    if infer_loadings:
        X = amplicons.copy()
        Y = model.addVars(m, n_amplicons, vtype=GRB.CONTINUOUS, lb=min_loading, ub=max_loading, name="Y")  
        for i in range(m):
            for p in range(n_amplicons):
                Y[i, p].Start = loading_start[i, p]
                Y[i, p].VarHintVal = loading_start[i, p]
    else:
        X = model.addVars(n_amplicons, n, lb=0, ub=max_amplicon_cn, vtype=GRB.INTEGER, name="X")      
        for p in range(n_amplicons):
            for j in range(n):
                X[p, j].Start = amplicon_start[p, j]
                X[p, j].VarHintVal = amplicon_start[p, j]
        Y = loadings.copy()

    Xpos = model.addVars(n_amplicons, n, vtype=GRB.BINARY, lb=min_loading, ub=max_loading, name="Y")  
    Ypos = model.addVars(m, n_amplicons, vtype=GRB.BINARY, lb=min_loading, ub=max_loading, name="Y")  

    # indicator variable showing when amplicons X is positive
    for p in range(n_amplicons):
        for j in range(n):
            model.addConstr(Xpos[p, j] * max_amplicon_cn >= X[p, j], f'Xpos_def_{p}_{j}')

    # indicator variable showing when loadings Y is positive
    for i in range(m):
        for p in range(n_amplicons):
            model.addConstr(Ypos[i, p] * max_loading >= Y[i, p], f'Ypos_def_{i}_{p}')
    
    # auxiliary variables representing objective terms
    obj_abs = model.addVars(m, n, vtype=GRB.CONTINUOUS, lb=0, name="obj_abs")   # abs value obj_abs
    
    # Constraints for absolute value (obj_abs)
    for i in range(m):
        for j in range(n):
            expr = data_matrix[i, j] - gp.quicksum(Y[i, p] * X[p, j]
                                                          for p in range(n_amplicons))
            model.addConstr(obj_abs[i, j] >= expr, name=f"obj_abs_pos_{i}_{j}")
            model.addConstr(obj_abs[i, j] >= -expr, name=f"obj_abs_neg_{i}_{j}")
    

    dot_term = 0
    n_pairs = max(1, scipy.special.comb(n_amplicons, 2))
    if lambda_dot > 0:
        for p1 in range(n_amplicons):
            for p2 in range(p1 + 1, n_amplicons):
                for j in range(n):
                    dot_term += X[p1, j] * X[p2, j]

    if cell_weights:
        # cells are weighted so that they all contribute similarly to the objective
        weights = 1 / np.sum(data, axis=1)
        weights /= np.mean(weights)    
    else:
        weights = np.ones(m)
    model.setObjective((gp.quicksum(obj_abs[i, j] * weights[i] for i in range(m) for j in range(n))
                       + gp.quicksum(lambda_x*Xpos[p, j] for p in range(n_amplicons) for j in range(n)) # L0 regularization on amplicons X
                       + gp.quicksum(lambda_x_l1*X[p, j] for p in range(n_amplicons) for j in range(n)) # L1 regularization on amplicons X
                       + gp.quicksum(lambda_y*Ypos[i, p] for i in range(m) for p in range(n_amplicons)) # L0 regularization on loadings Y
                       + (lambda_dot / n_pairs) * dot_term), # regularization on dot product of amplicons X
        GRB.MINIMIZE)
    model.setParam('OutputFlag', 0)
    model.setParam('Threads', n_threads)
    model.setParam('TimeLimit', time_per_run)
    return model, {'X':X, 'Y':Y}


def optimize_model(data, loading_start, amplicon_start, n_amplicons, min_loading, max_loading, max_amplicon_cn,
                           n_threads=7, time_per_run=120, lambda_x=0, lambda_y=0, lambda_x_l1=0, verbose=True,
                           cell_weights=False, lambda_dot=0):    
    converged = False   
    loadings = loading_start.copy()
    amplicons = amplicon_start.copy()
    prev_amplicons = np.zeros(shape=amplicons.shape)
    m, n = data.shape
    
    constants = {'data':data,
                 'n_amplicons':n_amplicons, 
                 'max_amplicon_cn':max_amplicon_cn,
                 'min_loading':min_loading,
                 'max_loading':max_loading,
                 'lambda_x':lambda_x,
                 'lambda_y':lambda_y,
                 'lambda_x_l1':lambda_x_l1,
                'n_threads':n_threads,
                'time_per_run':time_per_run,
                'cell_weights':cell_weights,
                'lambda_dot':lambda_dot}
    
    n_iter = 0
    objectives = []
    if verbose:
        print(datetime.now(), 'start')
    while not converged:
        n_iter += 1
        
        # fix amplicons and infer loadings
        model, variables = get_model(amplicons=amplicons, loading_start=loadings, **constants)
        model.optimize()
        loadings = np.zeros(shape=(m,n_amplicons))
        for i in range(m):
            for p in range(n_amplicons):
                loadings[i, p] = variables['Y'][i, p].X
        objectives.append(model.ObjVal)
        if verbose:
            print(datetime.now(), f'loadings step {n_iter}: model obj={model.ObjVal}, main obj={compute_main_objective(data, loadings, amplicons)}')

        # fix loadings and infer amplicons
        model, variables = get_model(loadings=loadings, amplicon_start=amplicons, **constants)
        model.optimize()
        amplicons = np.zeros(shape=(n_amplicons,n), dtype=int)
        for p in range(n_amplicons):
            for j in range(n):
                amplicons[p, j] = variables['X'][p, j].X
        objectives.append(model.ObjVal)
        if verbose:
            print(datetime.now(), f'amplicons step {n_iter}: model obj={model.ObjVal}, main obj={compute_main_objective(data, loadings, amplicons)}')

        if np.array_equal(prev_amplicons, amplicons):
            converged=True
        else:
            prev_amplicons = amplicons

    return {'loadings':loadings, 'amplicons':amplicons, 'iterations':n_iter, 'objective_value':model.ObjVal, 'objectives':objectives}

