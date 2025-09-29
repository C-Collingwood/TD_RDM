# Two-Drift_Toolbox
Draft toolbox for the two-drift race diffusion model, in MATLAB R2024b

The two-drift race diffusion model is an RL-EAM which works in a multi-alternate choice situation,
with an internal, time-dependent switch in drift-rate.

Included in this toolbox are four models:
Habit1_Race - *DEFAULT* The early/fast drift rate for this model is determined by H, an 
             measure of action expectation (S-R habits), which learns through
             an APE.
             The late/slow drift rate is determined by the sum of H and a 
             Q-value, from classic RL, which learns using an RPE.
Habit2_Race - As with Habit1_Race, but the temperature term for H adapts at t_2
RL2_Race   - This model uses two Q-values with different learning rates. The early
             drift rate applies Q1, and the later drift rate applies Q1 + Q2. No limitations
             are placed on the relative learning rates.
RL_Race    - Most closely resembles the RL-RDM proposed by Miletic et al. 2021, 
             using a single constant drift rate, based on a Q-value.

## How to use:
The end-user function is *fitTwoDriftRM.m*, which returns best fitting parameters, the
overall negative loglikelihood of the provided data, the loglikelihood
of each trial and an object containing all the model data.
An additional function *Simulate* is provided to create surrogate data.


## Getting Started
We recommend beginning with *Example_Script.m*, which demonstrates how fitTwoDriftRM.m can be used on the data of one participant provided by Hardwick et al. (2019).
This script walks through the required data-cleaning, calculates the t_1 heuristic and creates the required input data structures (and optional)
in order to run *fitTwoDriftRM*. 
This is followed by an example of how to use the *Simulate* function on the Habit1_Race algorithm to create data for 100 free-RT trials followed by 50 time-controlled.

### Input: Required
choice: A vector of choices made across the trials (e.g. \[1 2 3 1 4 3 4])
stimulus: A vector of the stimuli seen in each trial (e.g., \[1 2 4 1 3 3 4])
rt: A vector of reaction times in each trial  (in seconds)
outcome: A vector of outcomes (0/1) of each trial (e.g., \[1 1 0 1 0 1 1])
map: Matrix of stimulus mapping, Rows = pairings, Columns = \[stimulus, responseA, responseB] 
     (e.g., \[1,1,1;
             2,2,2;
             3,3,4;
             4,4,3])
time_cont: logical vector of flags for time_controlled (=1, default) or free trials (=0) (e.g., \[0 0 0 0 0 1 1])
non_decision: non-decision time, t_1, in seconds (= 0.1s, default)


### Optional Input (Name-Value arguments):
model: Habit1_Race (default), Habit2_Race, RL2_race, RL_race (STRING)
remap_trial: index of trial where map switches from A to B, (=1, default)
include_fit: logical vector for trials to be included in fitting
             procedure. (=ones(1,number trials), default).
weight: Weighting towards forced trials, 0-1, (=0.5,default)
new_cond: Trial number where a new condition starts, (return to map A, all 
          variables reset to initial values), (= [], default).
mult_start: How many start-points to run fitting from, using the MultiStart MATLAB package. 
            If this is a non-zero value, the function runs in parallel (=0, default (non-parallel, 1 start point)). 


### Output
PARAM: The best-fitting parameters in a vector. Order is model-dependent:
    1. "Habit1_Race"
        ..* \[aq, ah, bq, bh, t2, theta] 
    2. "Habit2_Race
        ..* \[aq, ah, bq, bh1, bh2, t2, theta] 
    3. "RL_Race"
        ..* \[aq, bq, theta]
    4. "RL2_Race"
        ..* \[aq2, aq1, bq2, bq1, t2, theta]

NLL: The total negative loglikelihood of the best-fitting parameters.
LL: A vector containing the loglikelihood of each trial that was included in fitting. Excluded trials are returned as NaN.
MODEL: An object containing all data regarding the best-fitting parameters, including:
    * MODEL.type: String of model name (e.g., "Habit1_Race")
    * MODEL.free_par: Cell array of parameter names (e.g., {'aq','ah','bq','bh','t2','theta'})
    * MODEL.par: Table of best-fitting parameter values, including t1 and s)
    * MODEL.values: Structure containing all Q/H and dq/dh values for each trial, (e.g., MODEL.values.Q)
    * MODEL.weight: The weighting value, w_c, used to balance time-controlled and free-RT LL values when calculating NLL.
    * MODEL.map: The mapping structure of stimulus to correct response.
    * MODEL.par_bound: The upper and lower boundaries on each free parameter.
    * MODEL.data: A table containing all the data used during the course of model learning and fitting, including:
        ..* .data.Choice
        ..* .data.Stimulus
        ..* .data.Reaction_Time
        ..* .data.Outcome
        ..* .data.Time-Controlled
        ..* .data.Remap_trial (0 = Map A, 1 = Map B)


## Additional Functions
This toolbox contains one other end-user function, *Simulate.m*, which can be used to create surrogate data.
Many of the input arguments (required and optional) are shared with fitTwoDriftRM.

### Input: Required
model: Habit1_Race (default), Habit2_Race, RL2_race, RL_race (STRING).
trials: total number of trials (free-RT and time-controlled).
parameter: a vector of parameters in the order \[alpha(s), beta(s), t_2, theta], (as with the PARAM output).
map: Matrix of stimulus mapping, Rows = pairings, Columns = \[stimulus, responseA, responseB] 
     (e.g., \[1,1,1;
             2,2,2;
             3,3,4;
             4,4,3])
non_decision: non-decision time, t_1, in seconds (= 0.1s, default)

### Optional Input (Name-Value arguments):
stimulus: A vector of the stimuli seen in each trial (e.g., \[1 2 4 1 3 3 4], default = random)
time_cont: logical vector of flags for time_controlled (=1, default) or free trials (=0) (e.g., \[0 0 0 0 0 1 1])
forced_rt: a vector (1,trial number) with the RT (seconds) at the elements corresponding to time-controlled, default = uniform\[0,2].
remap_trial: index of trial where map switches from A to B, (=1, default)
new_cond: Trial number where a new condition starts, (return to map A, all 
          variables reset to initial values), (= \[], default).
error: Value for output on error trials (= 0, default).
reward: Value for output on correct trials (=1, default).


### Output
MODEL object as above.


## Required MATLAB Packages
Optimization Toolbox
Statistics Toolbox
MATLAB R2024b