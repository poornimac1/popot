import numpy as np
import matplotlib.pyplot as plt

# For properly displaying the LaTeX legends
from matplotlib import rc
rc('text', usetex=True)


def monte_carlo(Nmc, f, sampling_func):
    means = [f(sampling_func()) for i in range(Nmc)]
    mu = sum(means,0.0) / Nmc
    variances = map((lambda x: (x-mu)**2),means)
    var = sum(variances, 0.0) / Nmc
    return [mu, var]

def generate_sample_IT(cdf_1):
    # Given a probability density function f, with an associated cumulative
    # distribution function cdf with a known reciprocal cdf_1
    # we can generate a sample according to f by
    # - generating a sample y following a uniform distribution in [0,1]
    # - returns cdf_1(y)
    return cdf_1(np.random.random())

def generate_sample_AR(pdf_f, pdf_g, draw_sample_g, c):
    # Generate a sample according to the accept/reject algorithm
    # pdf_f is the probability density function from which we want to draw a sample
    # pdf_g is the instrumental probability density function
    # draw_sample_g is the function thanks to which we can draw a sample from pdf_g
    # c is the constant such that :
    # Forall y : f(y) / (c g(y)) <= 1
    found = False
    while(not(found)):
        y = draw_sample_g()
        u = np.random.random()
        if(u <= pdf_f(y) / (c * pdf_g(y))):
            found = True
    return y


def cdf(x_min, x_max, N, samples):
    # [x_min, x_max] : bounds for the cumulative distribution
    # N : number of subsamples
    # samples : '''list''' of samples from which we want to compute the CDF
    val = np.zeros(N,2)
    for i in range(N):
        val[i,0] = x_min + i * (x_max - x_min)/(N-1)
        

##### Exponential distribution

def exp_pdf(x, params):
    # Exponential probability density function
    # Params : [Lambda]
    lbda = params[0]
    return lbda * np.exp(-lbda * x)

def exp_cdf(x, params):
    # Cumulative distribution
    # Params : [Lambda]
    lbda = params[0]
    return 1.0 - np.exp(-lbda * x)   

def exp_cdf_1(x, params):
    # Reciprocal Cumulative Distribution
    # Params : [Lambda]
    lbda = params[0]
    if(not lbda):
        raise Exception('Infinite mean !')
    return -np.log(1.0 - x) / lbda 
    
##### Bernouilli distribution

def bernouilli_pmf(x, params):
    # Bernouilly probability mass function
    # params = [p]
    # P(X=0) = 1 - p
    # P(X=1) = p
    if(x==0):
        return 1 - params[0]
    elif(x==1):
        return params[0]
    else:
        raise Exception("Bernouilli distribution works only on binary random variables")

def bernouilli_cdf(x, params):
    # Cumulative Distribution
    if( x < 0):
        return 0.0
    elif(x < 1):
        return 1.0 - params[0]
    else:
        return 1.0



if __name__ == "__main__":
    # Let's try our sampling methods
    plt.figure()

    ###########################
    # Exponential distributions
    plt.subplot(121)
    plt.title("Exponential distributions")
    x_min = 0
    x_max = 5
    for lbda in [0.5, 1.0, 1.5]:
        Nsamples = 1000
        X = []
        for i in range(Nsamples):
            X.append(generate_sample_IT(lambda x:exp_cdf_1(x,[lbda])))
        # Compute the Cumulative Distribution from the samples
        cdf = [[],[]]
        for i in range(300):
            xval = x_min + i * (x_max - x_min)/299.0
            cdf[0].append(xval)
            cdf[1].append(np.where(np.array(X) <= xval)[0].size/float(Nsamples))
        plt.plot(cdf[0], cdf[1])
        plt.xlim([x_min, x_max])
    # We do the same with the accept/reject
    # Say, we want to sample according to exp(lbda_0)
    # and we don't know its cdf_1
    # We take exp(lbda_1) from which we know the cdf_1 
    # with c = lbda_0 / lbda_1
    for lbda in [0.5, 1.0, 1.5]:
        Nsamples = 1000
        X = []
        # Define our Instrumental distribution
        # Take lbda1 = lbda/2.0 so that lbda1 <= lbda
        lbda_1 = lbda/2.0
        pdf_g = lambda x:exp_pdf(x, [lbda_1])
        draw_sample_g = lambda :generate_sample_IT(lambda x:exp_cdf_1(x,[lbda_1]))
        pdf_f = lambda x:exp_pdf(x, [lbda])
        # c equals to f(0)/g(0)=lbda / lbda1 (since lbda1 <= lbda)
        c = lbda / lbda_1
        for i in range(Nsamples):
            X.append(generate_sample_AR(pdf_f, pdf_g, draw_sample_g, c))
        # Compute the Cumulative Distribution from the samples
        cdf = [[],[]]
        for i in range(300):
            xval = x_min + i * (x_max - x_min)/299.0
            cdf[0].append(xval)
            cdf[1].append(np.where(np.array(X) <= xval)[0].size/float(Nsamples))
        plt.plot(cdf[0], cdf[1])
        plt.xlim([x_min, x_max])
    plt.legend([r'$IT - \lambda = 0.5$',r'$IT - \lambda = 1.0$',r'$IT - \lambda = 1.5$',
                r'$AR - \lambda = 0.5$',r'$AR - \lambda = 1.0$',r'$AR - \lambda = 1.5$'])


    plt.show()
