# Two-Drift_Toolbox
Draft toolbox for the two drift race model

The two-drift race model is an RL-EAM which works in a multi-alternate choice situation,
with an internal, time-dependent switch in drift rate.

Included in this toolbox are three models:
Habit_race - (DEFAULT) The early/fast drift rate for this model is determined by H, an 
             measure of action expectation (S-R habits), which learns through
             an APE.
             The late/slow drift rate is determined by the sum of H and a 
             Q-value, from classic RL, which learns using an RPE.
RL2_race   - This model uses two Q-values with different learning rates. The early
             drift rate applies Q1, and the later drift rate applies Q1 + Q2. No limitations
             are placed on the relative learning rates.
RL_race    - Most closely resembles the RL-RDM proposed by Miletic et al. 2021, 
             using a single constant drift rate, based on a Q-value.

How to use:

User function is fitTwoDriftRM.m, which returns best fitting parameters, the loglikelihood
used to calculate them and an object containing all the model data.

Input: Required
choice: A vector of choices made across the trials (e.g. [1 2 3 1 4 3 4])
stimulus: A vector of the stimuli seen in each trial (e.g., [1 2 4 1 3 3 4])
rt: A vector of reaction times in each trial  (in seconds)
outcome: A vector of outcomes (0/1) of each trial (e.g., [1 1 0 1 0 1 1])
map: Matrix of stimulus mapping, Rows = pairings, Columns = [stimulus, responseA, responseB] 
     (e.g., [1,1,1;
             2,2,2;
             3,3,4;
             4,4,3])
force: logical vector of flags for forced (=1) or free trials (=0)


Optional Input (Name-Value arguments):
model: Habit_race (default), RL2_race, RL_race (STRING)

remap_trial:index of trial where map switches from A to B, (=1, default)
non_decision: non-decision time, default uses heuristic described by
              approx_th.m
include_fit: logical vector for trials to be included in fitting
             procedure. (=ones(1,number trials), default).
weight: Weighting towards forced trials, (0-1, =0.5,default)
