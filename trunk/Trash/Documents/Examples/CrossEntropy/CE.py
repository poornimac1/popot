from sampling import *
import math

def algorithm_2_1(v0, S, sampling_func, pdf, params, verbose=0):
    # Algorithm 2.1 of "A Tutorial on the Cross Entropy Method"
    # v0 : initial parameter vector
    # S  : function to compute the performance
    # sampling_func : function to draw a sample from the pdf
    # params = [N, rho, gamma]
    t = 1
    v = v0
    N = params[0]
    rho = params[1]
    gamma = params[2]
    loop = True
    while(loop):
        ########
        # Step 2
        # Let's draw N samples from the pdf
        X = [sampling_func(v) for i in range(N)]
        # Compute their performances
        perfs = [(id, S(Xi)) for id,Xi in enumerate(X)]
        # Sort the performances
        sorted_perf = sorted(perfs, key=lambda elem: elem[1])
        # Compute the (1-rho) sample quantile of the performances
        gamma_t = sorted_perf[int(math.floor((1-rho)*N))][1]
        if(gamma_t > gamma):
            gamma_t = gamma
        ####
        # Step 3 : Compute the new parameter vector
        # The computation of the sum can be made faster
        # as we sum only on a part at the end of the sorted perf.
        vt = [0.0 for i in range(len(v0))]
        norm_factor = 0.0
        for i in range(N):
            if(perfs[i][1] >= gamma_t):
                w = pdf(X[i], v0) / pdf(X[i], v)
                norm_factor += w
                for j in range(len(v0)):
                    vt[j] = vt[j] + w * X[i][j]
        vt = [x/norm_factor for x in vt]
        if(gamma_t == gamma):
            loop = False
        else:
            if(verbose):
                print "t = %i ; gamma_t = %e ; v = [%s]" % (t, gamma_t,' '.join('{0:.3f};'.format(k) for k in vt))
            #print gamma_t
            #print vt
            #print ""
            t = t + 1
        v = vt
    if(verbose): print "Finished at time %i" % t
    return [v, gamma_t]


 

def algorithm_2_2(p0, S, sampling_func, params, verbose=0):
    # p0 : initial parameter vector
    # S : performance function
    # sampling func
    
    # params = [ N, rho]
    # where :
    # N : number of samples
    # rho : use the (1 - rho) sample quantile

    def is_degenerate(p):
        res = True
        for pi in p:
            res = res & ((pi == 0) or (pi == 1))   
        return res

    t = 1
    p = list(p0)
    N = params[0]
    rho = params[1]
    loop = True
    while(loop):
        ########
        # Step 2
        # Let's draw N samples from the pmf
        X = [sampling_func(p) for i in range(N)]
        # Compute their performances
        perfs = [(id, S(Xi)) for id,Xi in enumerate(X)]
        # Sort the performances
        sorted_perf = sorted(perfs, key=lambda elem:elem[1])
        # Compute the (1-rho) sample quantile of the performances
        gamma_t = sorted_perf[int(math.floor((1-rho)*N))][1]
        ########
        # Step 3
        #
        # Compute the new parameter vector
        pt = [0.0 for i in range(len(p0))]
        for i in range(int(math.floor((1-rho)*N)),N):
            pt = [sum(x) for x in zip(pt, X[sorted_perf[i][0]])]
        pt = [x/(rho*N) for x in pt]
        t = t + 1
        p = pt
        if(verbose): print gamma_t, pt
        if(is_degenerate(pt)):
            loop = False
    return p


def optimization(v0, f, sampling_func, stochastic_pb_solver, params, should_stop,  alpha=1.0, d = 5, verbose=0):
    # p0 : initial estimate
    # f : function to maximize
    # sampling_func : function to draw samples from f(.,v)
    # params = [N, rho]
    # where :
    # N : number of samples
    # rho : used to select the (1-rho) sample quantile
    # alpha : for smoothed updating
    # d : history over which we define the stopping criteria
    v = list(v0)
    t = 1
    N = params[0]
    rho = params[1]
    loop = True
    history_gamma = [0.0 for i in range(d)]
    while(loop):
        ########
        # Step 2
        # Generate N samples
        X = [sampling_func(v) for i in range(N)]
        # Compute their performances
        perfs = [f(Xi) for Xi in X]
        # Sort the performances
        sorted_perf = sorted(perfs)
        # Compute the sample (1-rho) quantile
        gamma_t = sorted_perf[int(math.floor((1-rho)*N))]
        # Use these samples to solve the stochastic
        # program :
        # max_v D(v) = max_v 1/N sum_1^N I{f(Xi) >= gamma_t} ln(f(Xi,v))
        # returning v*
        vstar = stochastic_pb_solver(X, perfs, gamma_t, v0, v)
        vt = [alpha * vstar_i + (1.0 - alpha) * v_i for (vstar_i,v_i) in zip(vstar, v)]
        history_gamma = [history_gamma[i] for i in range(1,d-1)] + [gamma_t]
        loop = not should_stop(history_gamma,t)
        t = t + 1
        v = vt
        print gamma_t
    return v
