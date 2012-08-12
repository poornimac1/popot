import sampling
import CE


def f(X):
    # Let's maximize
    return -((1.0 - X[0])**2 + 100.0 *(X[1] - X[0]**2)**2)

def my_sampling_func(v):
    # The weight of the connections are independently drawn from
    # an exponential distribution with means in v
    # i.e. the parameter of the exponential distributions are 1/v_i
    return [sampling.generate_sample_IT(lambda x: sampling.exp_cdf_1(x,[1.0/vi])) for vi in v]

def my_pdf(x, u):
    # The random variables being independent, the pdf of X
    # is the product of the individual pdfs
    return reduce(lambda x,y: x*y, [sampling.exp_pdf(xi, [1.0/ui]) for xi,ui in zip(x,u)])

def my_stochastic_solver(samples, performances, gamma, u, vt_1):
  norm_fac = 0.0
  vt = [0.0 for i in range(len(vt_1))]
  for (Xi, Sxi) in zip(samples,performances):
      norm_fac += (Sxi >= gamma) * my_pdf(Xi, u)/my_pdf(Xi,vt_1)
      vt = [vtj + Xij * (Sxi >= gamma) * my_pdf(Xi, u)/my_pdf(Xi,vt_1) for (vtj,Xij) in zip(vt,Xi)]
  vt = [vtj/norm_fac for vtj in vt]
  return vt

def should_stop(history, t):
    return (t >= 10)

N = 100
rho = 0.1
alpha = 0.1
v0 = [1.0, 1.00]
params = [N, rho, alpha]
p = CE.optimization(v0, f, my_sampling_func, my_stochastic_solver, params, should_stop)

print p
