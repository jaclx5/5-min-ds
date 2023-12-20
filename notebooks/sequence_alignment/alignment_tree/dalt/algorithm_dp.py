from .alignment import Alignment
from .algorithm import *

class AlgorithmDynamicProgramming(Algorithm):
    def run(self, aln:Alignment, max_steps:int):
        """
        Run at most `max_steps` steps of the brute force algorithm, or until it finds the solution.
        At the end of the run a tree representing a state of the algorithm is produced and can be
        graphycally represented.
        """
        aln.reset()

        # keeps track of the best alignment for each position coordinates
        scoreboard = {aln.coords: aln}

        solution = None
        expanded = aln

        if max_steps == 0:
            i = 0
        else:    
            for i in range(max_steps):
                # get the best non expanded node so far
                to_expand = aln.get_best_node_to_expand()
    
                if to_expand:
                    # get the aligment with the best score in the same (i, j) coordinates as the next
                    # to expand alignment
                    best_aln_coords = scoreboard.get(to_expand.coords, None)
    
                    # compare the node to expand with the best already expanded
                    # for the same (i, j) position
                    if best_aln_coords is not None and to_expand.score < best_aln_coords.score:
                        # if the node is no better just "ignore" it, i.e. mark it as expanded
                        # without actually expanding it
                        to_expand.color = COLOR_IGNORED_BOX
                        ignore = True
                    else:
                        # if the node is the same or better, update the score and expand it
                        scoreboard[to_expand.coords] = to_expand
                        expanded = to_expand
                        ignore = False
                        
                    to_expand.expand(ignore=ignore)
    
                else:
                    # by definition, the best leaf at the end of the expansion must be the solution
                    solution = aln.get_solution()
                    break

            i += 1

        # colour yellow each best score in their position
        for node in scoreboard.values():
            node.color = "#FF00FF"

        # colour green the latest expanded node
        if expanded:
            expanded.color = COLOR_EXPANDED_BOX

        # colour red the best unexplored node so far
        if solution:
            solution.color = COLOR_SOLUTION_BOX
        else:
            # colour red the best unexplored node so far
            # which will be the green box in the a step
            best_non_expanded = aln.get_best_node_to_expand()

            if best_non_expanded:
                best_non_expanded.color = COLOR_BEST_BOX

        return solution is not None, i