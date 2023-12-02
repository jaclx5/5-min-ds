from .alignment import Alignment
from .algorithm import *

class AlgorithmBruteForce(Algorithm):
    def run(self, aln:Alignment, max_steps:int):
        """
        Run at most `max_steps` steps of the brute force algorithm, or until it finds the solution.
        At the end of the run a tree representing a state of the algorithm is produced and can be
        graphycally represented.
        """
        aln.reset()

        solution = None
        expanded = aln

        i = 0
        for i in range(max_steps):
            # get the first node to expand
            expanded = aln.get_node_to_expand()

            if expanded:
                expanded.expand()

            else:
                # exits when no more expansion is possible
                # in the brute force case, by definition, the best leaf
                # at the end of the full expansion is the solution
                solution = aln.get_solution()

                assert solution is not None, "No solution found check algorithm for correctness!"

                solution.color = COLOR_SOLUTION_BOX
                break
                
        # colour green the latest expanded node
        if expanded:
            expanded.color = COLOR_EXPANDED_BOX

        return solution is not None, i
