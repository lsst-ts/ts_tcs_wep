"""
    Pure Python/Numpy implementation of the Nelder-Mead algorithm.
    Reference: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method

    Download from https://github.com/fchollet/nelder-mead.
    Modify this script to fit the interger step for need.
"""

import copy
import numpy as np


def feval(func, vars=()):
    """
    
    Evaluate the function.
    
    Arguments:
        func {[object]} -- Function to evaluate.
    
    Keyword Arguments:
        vars {tuple} -- Inputs of function (default: {()}).
    
    Returns:
        [double] -- Output of evaluated function.
    """
    return eval("func")(*vars)


def nelderMeadModify(func, x_start, args=(), step=0.1, no_improve_thr=10e-6, no_improv_break=10, max_iter=0,
                     alpha=1., gamma=2., rho=-0.5, sigma=0.5):
    """
    
    Optimization of the Nelder-Mead algorithm.
    
    Arguments:
        func {[object]} -- Function to optimize, must return a scalar score and operate over a 
                           numpy array of the same dimensions as x_start.
        x_start {[numpy array]} -- Initial position.
    
    Keyword Arguments:
        args {tuple} -- [description] (default: {()})
        step {float} -- Look-around radius in initial step (default: {0.1}).
        no_improve_thr {float} -- Break after no_improv_break iterations with an improvement lower 
                                  than no_improv_thr (default: {10e-6}).
        no_improv_break {int} -- Break after no_improv_break iterations with an improvement lower 
                                 than no_improv_thr (default: {10})
        max_iter {int} -- Always break after this number of iterations. 
                          Set it to 0 to loop indefinitely (default: {0}).
        alpha {float} -- Reflection parameter of the algorithm (default: {1.}).
        gamma {float} -- Expansion parameter of the algorithm (default: {2.}).
        rho {float} -- Contraction parameter of the algorithm (default: {-0.5}).
        sigma {float} -- Reduction parameter of the algorithm (default: {0.5}).

    Returns:
        [tuple] -- Best parameter array and best score for the evaluated function.
    """

    # init
    dim = len(x_start)

    vars = (x_start,) + args
    prev_best = feval(func, vars)

    no_improv = 0
    res = [[x_start, prev_best]]

    for i in range(dim):
        x = copy.copy(x_start)
        x[i] = x[i] + step

        vars = (x,) + args
        score = feval(func, vars)  

        res.append([x, score])

    # simplex iter
    iters = 0
    while 1:
        # order
        res.sort(key=lambda x: x[1])
        best = res[0][1]

        # break after max_iter
        if max_iter and iters >= max_iter:
            return res[0]
        iters += 1

        # break after no_improv_break iterations with no improvement
        # print '...best so far:', best

        if best < prev_best - no_improve_thr:
            no_improv = 0
            prev_best = best
        else:
            no_improv += 1

        if no_improv >= no_improv_break:
            return res[0]

        # centroid
        x0 = [0.] * dim
        for tup in res[:-1]:
            for i, c in enumerate(tup[0]):
                x0[i] += c / (len(res)-1)

        # reflection
        xr = x0 + alpha*(x0 - res[-1][0])

        vars = (xr,) + args
        rscore = feval(func, vars)

        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]
            res.append([xr, rscore])
            continue

        # expansion
        if rscore < res[0][1]:
            xe = x0 + gamma*(x0 - res[-1][0])

            vars = (xe,) + args
            escore = feval(func, vars)

            if escore < rscore:
                del res[-1]
                res.append([xe, escore])
                continue
            else:
                del res[-1]
                res.append([xr, rscore])
                continue

        # contraction
        xc = x0 + rho*(x0 - res[-1][0])

        vars = (xc,) + args
        cscore = feval(func, vars)

        if cscore < res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # reduction
        x1 = res[0][0]
        nres = []
        for tup in res:
            redx = x1 + sigma*(tup[0] - x1)

            vars = (redx,) + args
            score = feval(func, vars)

            nres.append([redx, score])
        res = nres


if __name__ == "__main__":
    pass
