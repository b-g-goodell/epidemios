"""Outbreak module
========================

Purpose
-------

This module is a toy implementation of a Gillespie algorithm realization of a stochastic ordinary differential equation
model of a disease outbreak in the standard SIR (susceptible, infectious, resistant) model.

Installation
------------

Matplotlib...

Hosting, Issues, Complaints, Bugs
---------------------------------

This code is being hosted and maintained on `my repo <https://github.com/b-g-goodell/epidemios>`_; please take out an
issue there and @b-g-goodell to get my attention if there is a problem.

Examples and Sample Code
------------------------

Lorem ipsum

License
-------

Copyright (c) 2018-2020, Brandon Goodell. The license for this code is under `GNU GPL v3.0
<https://github.com/b-g-goodell/trace/blob/master/LICENSE>`_.

"""
from random import random
from math import log, sqrt


def sample_index(inp):
    """ Sample an index from the list of non-negative values in inp from a pmf proportional to inp

    Example: if inp = [0.1, 0.2, 0.7], then we produce 2 with probability 70%, we produce 1 with probability 20%..."""
    assert isinstance(inp, list)
    assert len(inp) >= 1
    for _ in inp:
        assert _ >= 0.0
    target = random()*sum(inp)
    cum = inp[0]
    i = 0
    while cum < target and i + 1 < len(inp):
        i += 1
        cum += inp[i]
    if cum >= target:
        return i
    return None


def sample_delay(inp):
    """ Generate an exponential random variable independent of previous generated variables"""
    assert inp > 0.0
    return -1.0/inp * log(random())


class Outbreak(object):
    """ Outbreak class tracks parameters, state, runtime, and has a lookup of different possible events"""
    def __init__(self, inp):
        self.parameters = inp['parameters']
        self.runtime = inp['runtime']
        self.state = inp['initial state']
        self.events = inp['events']
        for f in self.events:
            assert callable(f[0])
            assert isinstance(f[1], list)

    def run(self):
        """ Execute a simulation"""
        result = dict()
        t = 0
        while t < self.runtime:
            # Compute the PMF from self.events
            pmf = [f[0]((self.state, self.parameters)) for f in self.events]
            # Extract the event consequences from self.events
            tmp = [f[1] for f in self.events]

            # Check the PMF is non-negative and has a positive sum
            assert all(_ >= 0.0 for _ in pmf)
            assert sum(pmf) > 0.0

            # Compute the delay parameter for generating wait time between events
            delay_parameter = float(sum(pmf))
            # Generate wait time
            event_delay = sample_delay(delay_parameter)
            # Update time
            t += event_delay
            # Generate an event identity
            event_identity = sample_index(pmf)
            # Get the event differential and apply it
            event_differential = tmp[event_identity]
            self.state = [sum(_) for _ in zip(event_differential, self.state)]
            # Record results
            result[t] = [self.state, event_delay, event_identity, event_differential]

        return result


def standard_SIR(fn, gestation_period=280.0, proportion_population_preggers=0.25525, human_lifespan=30842.5, sqrt_median_time_between_trips_outside_the_house=sqrt(2.0), social_distancing_modifier=0.0, probability_of_transmission_per_contact=0.01, median_time_before_recovery=28.0, proportion_of_people_with_capacity_to_recover_at_all=1.0, infection_fatality_period=14.0):
    """ The parameter values included as defaults are either magic or random numbers, or made up, or whatever; pass in your own parameters to make your own simulations!"""

    # We show how to compute parameters from source-able numbers like human gestation period.
    parameters = []  # parameters are per-capita

    birth_rate = proportion_population_preggers / gestation_period
    parameters += [birth_rate]

    death_rate = 1.0/human_lifespan  # 1.0/median human lifespan in days
    parameters += [death_rate]

    rate_of_interaction = 1.0/sqrt_median_time_between_trips_outside_the_house  # 1.0/median wait time between trips outside the house per person
    parameters += [rate_of_interaction * (1.0 - social_distancing_modifier) * probability_of_transmission_per_contact]

    recovery_rate = 1.0/median_time_before_recovery  # recovery time
    parameters += [recovery_rate * proportion_of_people_with_capacity_to_recover_at_all]

    ifr = 1.0/infection_fatality_period
    parameters += [ifr]

    runtime = 200.0
    state = [90, 10, 0]

    # This is a complete list of all possible events in a standard SIR model... sort of.
    # For example, we are allowing infectious people to give birth to infectious people (lol)
    # Note that if some event occurs with a very low probability, it can generally be discarded from the simulation
    # entirely. For example, in a population with an extremely long gestation period, there's no reason to even
    # include birth rates. Similarly, if the lifespan is very high compared to how rapidly an infection moves through
    # the population, natural death rates can be ignored too
    sus_births_sus = [lambda x: x[1][0]*x[0][0], [1, 0, 0]]
    sus_dies_naturally = [lambda x: x[1][1]*x[0][0], [-1, 0, 0]]
    sus_gets_infected = [lambda x: x[1][2]*x[0][0]*x[0][1], [-1, 1, 0]]
    res_births_sus = [lambda x: x[1][0]*x[0][2], [1, 0, 0]]
    inf_recovers = [lambda x: x[1][3]*x[0][1], [0, -1, 1]]
    inf_dies_naturally = [lambda x: x[1][1]*x[0][1], [0, -1, 0]]
    inf_dies_from_disease = [lambda x: x[1][4]*x[0][1], [0, -1, 0]]
    inf_births_inf_lol = [lambda x: x[1][0]*x[0][1], [0, 1, 0]]
    res_dies_naturally = [lambda x: x[1][1]*x[0][2], [0, 0, -1]]

    events = [sus_births_sus, sus_dies_naturally, sus_gets_infected, res_births_sus, inf_recovers, inf_dies_naturally]
    events += [inf_dies_from_disease, inf_births_inf_lol, res_dies_naturally]

    # Create outbreak and record results
    inp = dict()
    inp['parameters'] = parameters
    inp['runtime'] = runtime
    inp['initial state'] = state
    inp['events'] = events
    ollie = Outbreak(inp)
    result = ollie.run()
    with open(fn, "w+") as wf:
        wf.write("t,S,I,R\n")
    timeline = sorted(list(result.keys()))
    for next_t in timeline:
        next_state = result[next_t][0]
        with open(fn, "a+") as wf:
            line = str(next_t) + ","
            line += str(next_state[0]) + ","
            line += str(next_state[1]) + ","
            line += str(next_state[2]) + "\n"
            wf.write(line)

standard_SIR("simulation_results_standard_SIR.csv")
