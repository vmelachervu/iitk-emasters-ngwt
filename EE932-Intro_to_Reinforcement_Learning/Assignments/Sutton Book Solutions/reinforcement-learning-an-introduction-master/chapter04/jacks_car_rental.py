# copyright huseyinabanox@gmail.com

import time
import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns

LOC2_EXPECTED_RETURN = 2
LOC2_EXPECTED_REQUEST = 4
LOC1_EXPECTED_RETURN = 3
LOC1_EXPECTED_REQUEST = 3
MAX_EXPECTED_UPDATE = 11
FREE_CARS_MOVE_TO_LOC2 = 0

CAR_MOVE_COST = 2
CAR_PARK_COST = 0
RENT_PRICE = 10
MAX_CARS = 20
MAX_CARS_MOVED = 5
GAMMA = 0.9
PMF_MATRIX_MAP = {expected: [poisson.pmf(occurred, expected) for occurred in range(MAX_CARS + 1)]
                  for expected in [2, 3, 4]}
USE_EXPECTED_RETURNS = True


class State:
    def __init__(self, nLoc1, nLoc2):
        self.nLoc1 = nLoc1
        self.nLoc2 = nLoc2


class ProbableUpdate:
    def __init__(self, req1, preq1, req2, preq2, ret1, pret1, ret2, pret2):
        self.req1 = req1
        self.preq1 = preq1

        self.req2 = req2
        self.preq2 = preq2

        self.ret1 = ret1
        self.pret1 = pret1

        self.ret2 = ret2
        self.pret2 = pret2


def step(cars_moved, s, V):
    """

    :param cars_moved: negative actions are towards location 1
    :param s: state
    :param V: value matrix
    :return: action return
    """

    nLoc1 = s.nLoc1
    nLoc2 = s.nLoc2

    cars_moved_to_charge = cars_moved
    if cars_moved > 0:  # to location 2
        cars_moved_to_charge = max(cars_moved_to_charge - FREE_CARS_MOVE_TO_LOC2, 0)
        nLoc1 -= cars_moved
    elif cars_moved < 0:
        nLoc2 -= abs(cars_moved)

    action_return = -abs(cars_moved_to_charge) * CAR_MOVE_COST
    action_return -= ((nLoc1 > 10) + (nLoc2 > 10)) * CAR_PARK_COST

    if cars_moved > 0:  # to location 2
        nLoc2 = min(nLoc2 + cars_moved, MAX_CARS)
    elif cars_moved < 0:
        nLoc1 = min(nLoc1 + abs(cars_moved), MAX_CARS)

    for u in states_updates():
        req1, req2, ret1, ret2 = min(u.req1, nLoc1), min(u.req2, nLoc2), u.ret1, u.ret2
        prob = u.preq1 * u.preq2 * u.pret1 * u.pret2

        spNLoc1 = min(nLoc1 - req1 + ret1, MAX_CARS)
        spNLoc2 = min(nLoc2 - req2 + ret2, MAX_CARS)

        reward = (req1 + req2) * RENT_PRICE + GAMMA * V[spNLoc1, spNLoc2]
        action_return += prob * reward

    return action_return


def states():
    for nLoc1 in range(MAX_CARS + 1):
        for nLoc2 in range(MAX_CARS + 1):
            yield State(nLoc1, nLoc2)


def states_updates():
    for req1 in range(MAX_EXPECTED_UPDATE):
        for req2 in range(MAX_EXPECTED_UPDATE):
            if USE_EXPECTED_RETURNS:
                yield ProbableUpdate(req1, PMF_MATRIX_MAP.get(LOC1_EXPECTED_REQUEST)[req1],
                                     req2, PMF_MATRIX_MAP.get(LOC2_EXPECTED_REQUEST)[req2],
                                     LOC1_EXPECTED_RETURN, 1.0,
                                     LOC2_EXPECTED_RETURN, 1.0
                                     )
            else:
                for ret1 in range(MAX_EXPECTED_UPDATE):
                    for ret2 in range(MAX_EXPECTED_UPDATE):
                        yield ProbableUpdate(req1, PMF_MATRIX_MAP.get(LOC1_EXPECTED_REQUEST)[req1],
                                             req2, PMF_MATRIX_MAP.get(LOC2_EXPECTED_REQUEST)[req2],
                                             ret1, PMF_MATRIX_MAP.get(LOC1_EXPECTED_RETURN)[ret1],
                                             ret2, PMF_MATRIX_MAP.get(LOC2_EXPECTED_RETURN)[ret2]
                                             )


def policy_evaluation(V, pi):
    print("policy_evaluation started")

    theta = 0.0001
    delta = 1
    iteration = 0

    while delta > theta:
        delta = 0
        iteration += 1
        istart = time.time()

        for s in states():
            v = V[s.nLoc1, s.nLoc2]

            a = pi[s.nLoc1, s.nLoc2]
            action_return = step(a, s, V)
            V[s.nLoc1, s.nLoc2] = action_return

            delta = max(abs(action_return - v), delta)

        # if iteration % 10 == 0:
        print("policy_evaluation iteration {}, max delta='{}' in {} seconds"
              .format(iteration, delta, time.time() - istart))
    print("policy_evaluation completed")


def policy_improvement(V, pi):
    start = time.time()
    print("policy_improvement started")

    policy_stable = True

    for s in states():
        old_action = pi[s.nLoc1, s.nLoc2]
        actions = np.arange(-MAX_CARS_MOVED, MAX_CARS_MOVED + 1)  # negative actions goes to location 1
        action_returns = []

        for a in actions:
            if (0 <= a <= s.nLoc1) or (-s.nLoc2 <= a <= 0):
                action_return = step(a, s, V)
                action_returns.append(action_return)
            else:
                action_returns.append(-1e+10)

        new_action = actions[np.argmax(action_returns)]
        pi[s.nLoc1, s.nLoc2] = new_action

        if new_action != old_action: policy_stable = False

    print("policy_improvement completed in {} seconds".format(time.time() - start))

    return policy_stable


def policy_iteration(file_name='figure_4_2.png'):
    V = np.zeros((MAX_CARS + 1, MAX_CARS + 1))
    pi = np.zeros((MAX_CARS + 1, MAX_CARS + 1), dtype=int)

    _, axes = plt.subplots(2, 3, figsize=(40, 20))
    plt.subplots_adjust(wspace=0.1, hspace=0.2)
    axes = axes.flatten()
    plot_policy(axes, 0, pi)

    policy_stable = False
    iteration = 0
    while not policy_stable:
        iteration += 1

        policy_evaluation(V, pi)
        policy_stable = policy_improvement(V, pi)

        if iteration < 5:
            plot_policy(axes, iteration, pi)

        print("Policy iteration {} completed".format(iteration))

    plot_value_function(V, axes)

    plt.savefig('../images/' + file_name)
    plt.close()


def plot_value_function(V, axes):
    fig = sns.heatmap(np.flipud(V), cmap="YlGnBu", ax=axes[-1])
    fig.set_ylabel('# cars at first location', fontsize=30)
    fig.set_yticks(list(reversed(range(MAX_CARS + 1))))
    fig.set_xlabel('# cars at second location', fontsize=30)
    fig.set_title('optimal value', fontsize=30)


def plot_policy(axes, iteration, pi):
    fig = sns.heatmap(np.flipud(pi), cmap="PuBuGn", ax=axes[iteration])
    fig.set_ylabel('# cars at first location', fontsize=30)
    fig.set_yticks(list(reversed(range(MAX_CARS + 1))))
    fig.set_xlabel('# cars at second location', fontsize=30)
    fig.set_title('policy {}'.format(iteration), fontsize=30)


def solve_jacks_car_rental():
    policy_iteration()


def solve_jacks_car_rental_e_4_7():
    global FREE_CARS_MOVE_TO_LOC2
    FREE_CARS_MOVE_TO_LOC2 = 1
    global CAR_PARK_COST
    CAR_PARK_COST = 4
    policy_iteration(file_name='figure_4_2_e_4_7.png')


if __name__ == '__main__':
    # solve_jacks_car_rental()
    solve_jacks_car_rental_e_4_7()
