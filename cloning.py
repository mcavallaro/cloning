import numpy as np
from random import choice, sample
# heapq provides the Heap Queue structure which allows to manipulate efficiently sorted list.
# Here we use a list pof copies of the system, sorted according to their next time of evolution.
import heapq 

# procedure to remove a element of index i from a heap, keeping the heap structure intact
def heapq_remove(heap, index):
    # Move slot to be removed to top of heap
    while index > 0:
        up = int((index + 1) / 2 - 1)
        heap[index] = heap[up]
        index = up
    # Remove top of heap and restore heap property
    heapq.heappop(heap)


def escape(state):
    if state == 1:
        return 5
    else:
        return 1


def flip(state):
    if state == 1:
        new_state = 0
    else:
        new_state = 1
    return new_state


def cloning_factor(s, state, dt):
    return np.exp(state * s)


DIM_ENSEMBLE = 400
observation_time = 400.

def main(s):

    t = 0
    C = 0.

    # initialization of the population
    # each clone in the ensemble is described by a 3-uple (time,dt,state) 
    #   time  = next time at which it will evolve
    #   dt = time since last evolution
    #   state = 0 or 1 = empty or occupied 

    ensemble = [(0., 0., 0.) for count in range(DIM_ENSEMBLE)]

    # orders the population into a Heap Queue
    heapq.heapify(ensemble)
 
    while t < observation_time:
        # we pop the first element of populat, which is always the next to evolve
        (t, dt, state) = heapq.heappop(ensemble)
        Y = cloning_factor(s, state, dt)
        # the copy we poped out is to be replaced by p copies
        y = int(Y + np.random.random()) 

        if y == 0:
            # one copy chosen at random replaces the current copy
            to_clone = choice(ensemble)
            heapq.heappush(ensemble, to_clone)
        elif y == 1:
            # the current copy is evolved without cloning
            new_state = flip(state)

            # interval until next evolution
            Deltat = np.random.gamma(1, escape(new_state))
            to_clone = (t + Deltat, Deltat, new_state)
            heapq.heappush(ensemble, to_clone)
        else: # p>1 : make y clones ; population size becomes N+y-1 ; remove y-1 clones uniformly 
            pcount = y
            while pcount > 0:
                pcount -= 1
                new_state = flip(state)
                Deltat = np.random.gamma(1, escape(new_state)) # interval until next evolution
                toclone = (t + Deltat, Deltat, new_state)
                heapq.heappush(ensemble, toclone)

            # we first chose uniformly the p-1 distinct indices to remove, among the N+p-1 indices 
            listsize = DIM_ENSEMBLE + y - 1
            indices = sample(range(listsize), y - 1) 
            # the list of indices to remove is sorted from largest to smallest, so as to remove largest indices first
            indices.sort(reverse = True)
            for i in indices:
                heapq_remove(ensemble, i)

        C += np.log((DIM_ENSEMBLE + Y - 1.) / DIM_ENSEMBLE)

    return C / t


if __name__ == '__main__':
  main()
