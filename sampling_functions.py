# Imports

import numpy as np
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.special import gamma as gamma_func
from scipy.special import kv
from scipy.linalg import expm
import time as t

# 16th October Session

# reply to Joe Johnson on email

# We need to simulate the gamma process
'''
We do a rejection sampling approach

'''


class GammaDistr:
    def __init__(self, alpha, beta):
        self.alpha = alpha
        self.beta = beta
        self.distribution = 'Gamma'

        self.T = None
        self.start_time = None
        self.sample_size = None
        self.rate = None

    def set_process_conditions(self, t0, T, sample_size):
        self.start_time = t0
        self.T = T
        self.sample_size = sample_size
        self.rate = 1./(T-t0)


class NormalDistr:
    def __init__(self, mean, std, secondary_distr):
        self.mean = mean
        self.std = std
        self.distribution = 'Normal'

        self.secondary_distr = secondary_distr

    def NormalGammaPDF(self, x):
        # plot each term to find issues

        gamma_distr = self.secondary_distr


        # alpha = gamma_distr.alpha
        # beta = gamma_distr.beta
        #
        # t1 = (beta ** alpha)*np.sqrt(1)/(gamma_func(alpha)*np.sqrt(2*np.pi))
        # t2 = np.array(times*(alpha-1/2))
        # t3 = np.exp(-beta*times)
        # t4 = np.exp(-0.5*times*(x - self.mean)**2)

        # return t1*t2*t3*t4

        #
        #
        #
        t1 = 2*np.exp(self.mean*x/self.std**2)
        #
        # print(gamma_distr)
        #

        delta = 2*self.std**2/gamma_distr.beta + self.mean**2
        tau = 1/gamma_distr.beta - 0.5



        t2 = np.power(gamma_distr.beta, 1/gamma_distr.beta)*np.sqrt(2*np.pi*(self.std**2))*gamma_func(1/gamma_distr.beta)
        t3 = (np.abs(x)/(self.std**2)) * np.sqrt(delta)
        t4 = (1/gamma_distr.beta) - 0.5
        t5 = (1./self.std**2) * np.sqrt(self.mean**2 + (2*(self.std**2)/gamma_distr.beta)*np.abs(x))

        # plt.plot(np.arange(0, len(t5)), t5)
        # plt.show()

        # print(t2)
        # print(t4)
        # print(kv(t4, t5))

        return (t1 / t2) * (x**2/delta)**(tau/2) * kv(tau,t3)

        # return (t1/t2)*np.power(t3, t4) * kv(t4, t5)

debug = False

class DistributionSimulator:

    def __init__(self, DistributionObject): #*OtherObject
        self.distribution = DistributionObject.distribution
        self.distr_obj = DistributionObject
        # self.secondary_distr = OtherObject

        self.sorted_process_set = None
        self.time_arr = None
        self.process_path = None
        self.NG_jump_series = None

    def plot_simulation_distribution(self, process_endpoints):

        plt.figure(1)
        fig, ax = plt.subplots()
        xx, bins, p = ax.hist(process_endpoints, bins=200, density=True)

        if self.distribution == 'Gamma':  # this was for future if we plot a diff distrs end. classes now handle this - This is to overlay the PDF
            T = self.distr_obj.T
            x = np.linspace(0, 10, 1000)
            # x = self.

            # TODO: had self.distr_obj as gamma obj before
            shape = (self.distr_obj.alpha ** 2) * T / self.distr_obj.beta
            rate = self.distr_obj.beta / self.distr_obj.alpha
            y = stats.gamma.pdf(x, a=shape, loc=0, scale=rate) #TODO: put loc argument = 0?

        if self.distribution == 'Normal':
             x = np.linspace(-10, 10, 5000)
             #x = self.NG_jump_series

             y = self.distr_obj.NormalGammaPDF(x)

            # x, y = 0, 0

        ax.plot(x, y)
        # for item in p:
        #     item.set_height(item.get_height()/sum(xx))
        fig.suptitle('Endpoints for {} Process Simulation'.format(self.distribution if self.distribution != 'Normal' else 'Normal Gamma'))
        fig.show()
        return fig, ax

    def tractable_inverse_gamma_tail(self, alpha, beta, x):
        return 1 / ((alpha / beta) * (np.exp(beta * x / alpha ** 2) - 1))

    def acceptance_probability(self, alpha, beta, x):

        return (1 + alpha * x / beta) * np.exp(-alpha * x / beta)

    def process_simulation(self, DistributionObject, *prev_sim_data):

        if self.distribution == 'Gamma':
            # start_time = t.time()
            T = DistributionObject.T
            t0 = DistributionObject.start_time
            sample_size = DistributionObject.sample_size
            alpha = DistributionObject.alpha
            beta = DistributionObject.beta

            exp_rvs = np.random.exponential(scale=DistributionObject.rate, size=sample_size)
            poisson_epochs = np.cumsum(exp_rvs)

            # generate_jump_sizes = np.vectorize(self.tractable_inverse_gamma_tail)
            # generate_acceptance_probabilities = np.vectorize(self.acceptance_probability)

            jump_sizes = self.tractable_inverse_gamma_tail(alpha, beta, poisson_epochs) #NOTE - used to be caling the np vectorized functions
            acceptance_probabilities = self.acceptance_probability(alpha, beta, jump_sizes)


            # Now need to use this acceptance probability array to accept or reject the samples
            # samples = [0]

            # rnd = np.random.choice([0,1], size = acceptance_probabilities.shape, p = acceptance_probabilities)

            rnd_gen = np.random.default_rng()
            unif_r = rnd_gen.random(len(acceptance_probabilities))

            # accept samples who's acceptance probs are higher than this randomly generated list of numbers

            samples = np.where(acceptance_probabilities > unif_r, jump_sizes, 99999)
            samples = samples[samples < 99999]
            if t0 == 0: # if time range is 0-> 1, we start at initial point (0,0). Else, we do not insert a 0
                samples = [0] + list(samples)
            else:
                samples = list(samples)


            # for i in range(0, len(jump_sizes)): # TODO: SPEED UP!
            #     accept = random.choices([1, 0], weights=[acceptance_probabilities[i], 1 - acceptance_probabilities[i]],
            #                             k=1)
            #     if accept[0] == 1:
            #         samples.append(jump_sizes[i])
            if len(samples) == 0:
                jump_times = []
            else:
                if t0 == 0:
                    jump_times = [t0] + list(np.random.uniform(t0, T, len(samples)-1)) # jump times are uniform spread, property of exponen distr
                else:
                    jump_times = list(np.random.uniform(t0, T, len(samples)-1))

            # if len(samples) != len(jump_times):
            #     print('Length Issues')
            jumps_and_times = zip(samples, jump_times)

            self.sorted_process_set = sorted(jumps_and_times, key=lambda x: x[1])
            self.process_path = np.cumsum([jump_time_set[0] for jump_time_set in self.sorted_process_set])
            self.time_arr = [jump_time_set[1] for jump_time_set in self.sorted_process_set]


            # TODO: Ask about the case when time just starts
            # print('Time Elapsed for a Single Gamma Proces Path = {}'.format(t.time() - start_time))
            return self.process_path, self.time_arr, self.sorted_process_set


        if self.distribution == 'Normal': # we want to simulate a N-G process. Run a gamma first, then put it through a normal
            # DO NORMAL DISTR FUNCTIONS
            mean = DistributionObject.mean
            std = DistributionObject.std

            # Create a GammaDistr Object (with necessary params). Use this to make DistrSim Object. Run the Gamma Sim, then return the sorted process set
            hidden_gamma_distr = DistributionObject.secondary_distr
            hidden_gamma_sim = DistributionSimulator(hidden_gamma_distr)
            hidden_gamma_sim.process_simulation(hidden_gamma_distr)
            gamma_process_set = hidden_gamma_sim.sorted_process_set

            # make the gamma process set the sorted process set of the hidden_gamma_distr

            # process_set = prev_sim_data # this is the sorted process set of the gamma which we just ran


            normal_gamma_jump_series = []
            jump_time = list(zip(*gamma_process_set))
            # for jump_time in gamma_process_set:  # process_set is a set of (jump, time). Has length Ns
            normal_gamma_jump_series.append(np.random.normal(loc=mean * np.array(jump_time[0]), scale=(std) * np.sqrt(np.array(jump_time[0]))))

            self.process_path = np.cumsum(normal_gamma_jump_series)
            self.time_arr = [tuple[1] for tuple in gamma_process_set]

            self.NG_jump_series = normal_gamma_jump_series

            return self.process_path, self.time_arr


    def process_endpoint_sampler(self, iterations, DistributionObject, **kwargs): #TODO:

        '''
        Issue: So, below there is a class instance called gamma_sim and another called gamma_obj
                Right now these are global variables as I have defined these outside the function
                But, I want to pass these class instances to this function process_endpoint_sampler and so
                I tried to use **kwargs (and *args) but when I do so, I cannot figure out how to assign each of the variables in
                the kwargs to gamma_sim and gamma_obj within this function
        '''

        # This is not right, but for example - if kwargs is specified then I want this to happen:
        # gamma_sim = kwargs[0]
        # gamma_obj = kwargs[1]

        process_endpoints = []

        # start = t.time()
        for i in range(0, iterations):
            # if DistributionObject.distribution == 'Normal':
            # #     gamma_sim.process_simulation(gamma_obj) #we generate a new gamma sim to sample from each iteration
            #prev_sim_data = gamma_sim.process_simulation(gamma_obj).sorted_process_set  # TODO: test this. If not work, try and re-do the whole gamma sim each iteration

            #gamma_sim.process_simulation(gamma_obj) #TODO: Right Now This Shadows Name of Outer Scope - so

            self.process_simulation(DistributionObject) #, gamma_sim.sorted_process_set) # this gamma sim is from outer scope

            process_endpoints.append(self.process_path[-1])
        # print('Elapsed time for Endpoint Sampling: {}'.format(t.time() - start))
        return process_endpoints


class StateSpaceSimulation:
    """A Class Enabling Definition of State Space Model to Forward Simulate a Process With, And
    Provides Functionality to Simulate & Observe Such Process
    """

    # calculate jumps in interval
    # calculate mean vec and cov mat
    # sample for n (stoc int === noise)
    # incriment state vector by X_tb = exp(A(tb - ta) Xta + n

    def __init__(self, DistrSimObject):
        """Pass DistrSim Object, who's attributre distr_obj is our Distribution ie. here, it will be the Gamm distr"""
        # Have to Pass GammaSim object with a simulation already ran
        self.sorted_obs_times = DistrSimObject.time_arr #takes the initial observation times from the initial distrsimObj, which is a simulation based on our initial gamma object

        self.distr_obj = DistrSimObject.distr_obj
        self.distr_sim_obj = DistrSimObject

        self.X = np.zeros((2,1))
        self.t_0 = 0

        self.A = np.zeros((2,2))
        self.h = np.zeros((2))

        self.NG_mean = 0
        self.NG_var = 0

    def set_NG_conditions(self, mean, var):
        self.NG_mean = mean
        self.NG_var = var

    def define_A_matrix(self, flatterned_A):
        if len(flatterned_A) == 1:
            self.A = flatterned_A[0]
        else:
            self.A[0][0] = flatterned_A[0]
            self.A[0][1] = flatterned_A[1]
            self.A[1][0] = flatterned_A[2]
            self.A[1][1] = flatterned_A[3]

    def define_model_h(self, flatterned_h):
        self.h[0] = flatterned_h[0]
        self.h[1] = flatterned_h[1]


    def generate_jumps_between_time_intervals(self, t_a, t_b):
        gamma_obj_small_interval = GammaDistr(self.distr_obj.alpha, self.distr_obj.beta) #           nCheck time intervals of gamma obj are done correctly
        #give our gamma_object over this small time interval the alpha and beta of the distr_obj; which is the distribution object behind the gamma_sim passed when creating the first instance of the class
        gamma_obj_small_interval.set_process_conditions(t_a, t_b, self.distr_obj.sample_size)

        ### DistrSimObj = self.distr_sim_obj

        gamma_sim_small_interval = DistributionSimulator(gamma_obj_small_interval)
        path, time, jump_time_set = gamma_sim_small_interval.process_simulation(gamma_obj_small_interval)

        # DistrSimObj.process_simulation(gamma_obj_small_interval) #, gamma_obj_small_interval) # this overwrites the existing process set for the DistrSimObj w a new set of jumps for new gamma_obj_small_interval
        ### self.distr_sim_obj = DistrSimObj # update the DistrSimulationObject attribute of the state space simulator

        self.distr_sim_obj = gamma_sim_small_interval #Update the distr_sim_obj attribute with the recent run for small interval.

        ### process_path, sorted_times, sorted_set = self.distr_sim_obj.process_simulation(self.distr_obj)

        # each time generate_jumps_btw_time_intervals is called in each iter of run_state_space_sim,
        #we need to re-instate our gamma_objject (over the time interval minT, ta, to maxT tb)
        #which will allow us to THEN call aglgo 2 ie run a gamma process sim over ta to tb
        #we then return these jumps

    def calculate_jumps_mean_and_cov(self, start_time, end_time, Mat):
        """Calculate the mean vector given a collection of gamma jumps, gamma jump times, and a time interval"""

        #sum over time interval, matrix exponential (end time - jump time) x associated_gamma_jump
        mean_vec = np.array([ [float(0)],
                              [float(0)] ])
        cov_mat = np.array([ [float(0),float(0)],
                             [float(0),float(0)] ])

        # print()
        # print(self.distr_sim_obj.sorted_process_set)

        for jumpsize, jumptime in self.distr_sim_obj.sorted_process_set: # we get the sorted process set for the gamma sim on small interval
            ft_Vi = self.calc_matrix_exponent(end_time, jumptime, Mat)
            # print('ft_Vi:')
            # print(ft_Vi)

            mean_vec += (ft_Vi * jumpsize)
            cov_mat += (ft_Vi * ft_Vi.T * jumpsize)
        return mean_vec, cov_mat

    def sample_stoc_int(self, mean, cov):
        m = mean * self.NG_mean
        cov = cov * self.NG_var

        # cholesky decomposition for sampling of noise
        try:
            cov_chol = np.linalg.cholesky(cov)
            n = cov_chol @ np.column_stack([np.random.normal(size=2)]) + np.column_stack([m])
        except np.linalg.LinAlgError:
            # truncate innovation to zero if the increment is too small for Cholesky decomposition
            n = np.zeros(2)
        # print(np.shape(n))
        return n

    def update_state_vec(self, n, start_time, end_time, Mat):
        A = self.A
        if Mat:
            multiplier = expm(A * (end_time - start_time))
        else:
            multiplier = np.exp(A * (end_time - start_time))

        self.X = multiplier*self.X + n # TODO: Ask


    def calc_matrix_exponent(self, t, jump_time, Mat): # This is the ft(Vi) CALCULATION - returns the ft(Vi)
        """Calculate the matrix exponent given an A matrix, time interval end, and times within such time interval"""
        # we need to calc the matrix exponent for every summation in our time interval, as we have
        #exp(A(t - Vi)) where Vi are the jump times (of the gamma sim), and t is the END TIME of our interval.
        #-> this finds later use in our M
        # print('A:')
        # print(self.A)
        #
        # print('h:')
        # print(self.h)

        A = self.A
        h = self.h
        return expm(A*(t - jump_time))@np.reshape(np.array(h), (2,1)) if Mat else np.reshape(np.array(h), (2,1))*np.exp(A*(t-jump_time)) # TODO: have ft_VI something else if we do not use A matrix


    def run_state_space_simulation(self, Mat):
        # debug = True
        x0_evolution = [0]
        x1_evolution = [0]
        for i in range(0,len(self.sorted_obs_times)-1):
            start_time = self.sorted_obs_times[i]
            end_time = self.sorted_obs_times[i+1]

            self.generate_jumps_between_time_intervals(start_time, end_time)
            mean_vec, cov_mat = self.calculate_jumps_mean_and_cov(start_time, end_time, Mat)
            n = self.sample_stoc_int(mean_vec, cov_mat)
            self.update_state_vec(n, start_time, end_time, Mat)

            x0_evolution.append(self.X[0][0])
            x1_evolution.append(self.X[1][0])

        fig1, ax1 = plt.subplots()
        ax1.scatter(self.sorted_obs_times, x0_evolution, color='r', s=7, zorder=2)
        ax1.plot(self.sorted_obs_times, x0_evolution, zorder=1)
        fig1.suptitle('Evolution of State X0')
        plt.xlabel('Time')
        plt.ylabel('X0')
        plt.show()

        fig2, ax2 = plt.subplots()
        ax2.plot(self.sorted_obs_times, x1_evolution)
        fig2.suptitle('Evolution of State X1')
        plt.xlabel('Time')
        plt.ylabel('X0')
        plt.show()



        # plt.plot(self.sorted_obs_times, x0_evolution)
        # plt.title('Evolution of state vector x0')
        # plt.show()

        return fig1, fig2



def plotter(x, y, title, x_lbl, y_lbl):
    plt.figure(99)
    plt.step(x,y)
    plt.title(title)
    plt.xlabel(x_lbl)
    plt.ylabel(y_lbl)
    plt.show()
    return plt


"""
Now I want to pass our state vector Xt via some system.
dXt = A Xt dt + h dZt
We can define A and h to be state transition matrices 

We get Xt = exp(At) Xo + stoc.int(f_t)

2.2.9 states the stoc.int collapses to a sum, where the int is the effect of passing each of the jumps dZi
through the system at time Vi for Vi <= t

to Forward simulate (for the Variance Gamma):
sample from the integral  - because I(ft) is a guass rv; sum of gauss rvs = gauss rv
mean vector is given by: 3.1.2
Cov is given by : 3.1.3

The key eqns 



"""


if __name__ == '__main__':

    '''
    Generate Gamma Distr Object & Set The Conditions For The Process
    '''
    gamma_obj = GammaDistr(alpha=1, beta=1)
    gamma_obj.set_process_conditions(t0=0, T=1, sample_size=300) # sample_size refers to number of random times we generate

    '''
    Create a Simulation Object With the Parameters From the Existing Gamma Distribution Object
    Run The Simulation on This Object & Gather the Path, Times and set of Jumps and Times together
    '''
    gamma_sim = DistributionSimulator(gamma_obj) # create our simulation
    path, time, jump_time_set = gamma_sim.process_simulation(gamma_obj) # run the simulation. after this comment, the gamma_sim object has the jump_time_sets (sorted process set) as an attribute
    plt = plotter(time, path, 'Gamma Process Simulation TEST', 'Time', 'Value')

    '''
    Now Call .process_endpoint_sampler to Run The Sim 10,000 Times For Our Gamma Object (which defines the parameters of the gamma distr)
    & Plot Histogram of Endpoints. NOTE - in streamlit, we will want the sliders to correspond to changing the params of gamma_obj, and re-making the gamma obj (and gamma sim?)
    '''
    fig, ax = gamma_sim.plot_simulation_distribution(gamma_sim.process_endpoint_sampler(10000, gamma_obj))


    '''
    Define a Normal Distr Object To Be Used to Create a Normal Gamma Simulation Object
    '''
    normal_obj = NormalDistr(mean=0, std=np.sqrt(1), secondary_distr=gamma_obj) # We MUST define the secondary distribution of this NG sim ie. define the gamma distr to use
    #re-write a new gamma_obj to use?
    #gamma_obj_2 = GammaDistr(alpha=0.5, beta=0.5)
    normal_gamma_sim = DistributionSimulator(normal_obj)

    '''
    Run The Process Simulation, Using The Normal Obj & The Sorted Process Set of The Desired Gamma Simulation Object to Sample From
    '''
    path, time = normal_gamma_sim.process_simulation(normal_obj) #, gamma_sim.sorted_process_set)
    plotter(time, path, 'Normal Gamma Process Sim TEST', 'Time', 'Value')

    '''
    Now Call .process_endpoint_sampler With The Desired Normal Object, Gamma Object & Gamma Simulation Object to Sample 
    Process Endpoints & Plot
    '''
    fig, ax = normal_gamma_sim.plot_simulation_distribution(normal_gamma_sim.process_endpoint_sampler(10000, normal_obj)) #, G_obj=gamma_obj, G_sim_obj=gamma_sim))


    '''
    State Space Simulations
    
    -> Define an A matrix: [ [a1, a2], [a3, a4] ]
    -> Define a H matrix: [h1, h2]
    
    '''

    state_simulation = StateSpaceSimulation(gamma_sim) # Create the State Space class with our Gamma Sim object; gamma_sim = DistributionSimulator(gamma_obj)

    state_simulation.define_model_h([0,1])
    theta = -1
    # state_simulation.define_A_matrix([0,1,0,theta])
    state_simulation.define_A_matrix([-0.01])
    state_simulation.set_NG_conditions(normal_obj.mean, (normal_obj.std)**2)

    state_simulation.run_state_space_simulation(Mat=False)











# path, times, jumps = gamma_process_sim(alpha=1, beta=0.05, sample_size=50, T=1)
# plotter(times, path, 'Gamma Process Simulation', 'Time', 'Value')
#
# normal_gamma_path, times = generate_normal_gamma_process(0, 0.5, jumps, times)
# plotter(times, normal_gamma_path, 'Normal Gamma Simulation', 'Time', 'Value')

# TODO: Ask Simon about if we start at (0,0)

# TODO: Make classes, make functions dynamic; make the actions distribution independent eg. sampling endpoints should be the same function for any distr



# OLD CODE:


# #DONE - MOVED
# def tractable_inverse_gamma_tail(alpha, beta, x):
#     return 1/((alpha/beta) * (np.exp(beta*x/alpha**2) - 1))
#
# #DONE - MOVED
# def acceptance_probability(alpha, beta, x):
#     return (1 + alpha*x/beta)*np.exp(-alpha*x/beta)
#
# #DONE - MOVED
# def gamma_process_sim(alpha, beta, sample_size, T):
#
#     exp_rvs = np.random.exponential(scale=T, size=sample_size)
#     poisson_epochs = np.cumsum(exp_rvs)
#
#     generate_jump_sizes = np.vectorize(tractable_inverse_gamma_tail)
#     generate_acceptance_probabilities = np.vectorize(acceptance_probability)
#
#     jump_sizes = generate_jump_sizes(alpha, beta, poisson_epochs)
#     acceptance_probabilities = generate_acceptance_probabilities(alpha, beta, jump_sizes)
#
#     # Now need to use this acceptance probability array to accept or reject the samples
#     samples = list()
#     for i in range(0, len(jump_sizes)):
#         accept = random.choices([1,0], weights=[acceptance_probabilities[i], 1-acceptance_probabilities[i]], k=1)
#         if accept[0] == 1:
#             samples.append(jump_sizes[i])
#
#     jump_times = np.random.uniform(0,T,len(samples))
#
#     jumps_and_times = zip(samples, jump_times)
#
#
#     sorted_gamma_path = sorted(jumps_and_times, key = lambda x: x[1])
#
#     process_path = np.cumsum([jump_time_set[0] for jump_time_set in sorted_gamma_path])
#     time_arr = [jump_time_set[1] for jump_time_set in sorted_gamma_path]
#
#      # TODO: Ask about the case when time just starts
#
#     return process_path, time_arr, sorted_gamma_path
#
# # Make this a function to sample any process's endpoints
# def sample_gamma_process(alpha, beta, sample_size, iterations):
#     gamma_endpoints = []
#     for i in range(0, iterations):
#         process_path, time_arr, sorted_process_set = gamma_process_sim(alpha, beta, sample_size, T=1)
#         gamma_endpoints.append(process_path[-1])
#
#     return gamma_endpoints
#
#
# def plot_process_endpoints_samples(process_endpoints, gamma_distribution):
#     if gamma_distribution: # this was for future if we plot a diff distrs end. classes now handle this - This is to overlay the PDF
#         x = np.linspace(0, 20, 1000)
#         shape = (alpha**2) * t / beta
#         rate = alpha / beta
#         y = stats.gamma.pdf(x, a=shape, scale=rate)
#
#     fig, ax = plt.subplots()
#     ax.hist(process_endpoints, bins=300, density=True)
#     ax.plot(x, y)
#
#     return fig, ax
#
# #DONE - MOVED
# def generate_normal_gamma_process(mean, std, gamma_jump_series, times):
#     normal_gamma_jump_series = []
#     for jump in gamma_jump_series: # jump is a set of (jump, time)
#         normal_gamma_jump_series.append(np.random.normal(loc=mean*jump[0], scale=(std**2)*jump[0]))
#     return np.cumsum(normal_gamma_jump_series), times
#
# #







