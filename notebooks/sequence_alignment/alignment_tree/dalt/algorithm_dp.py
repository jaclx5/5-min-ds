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

        # keeps track of the best score for each position coordinates
        scoreboard = {aln.coords: aln.score}

        solution = None
        to_expand = aln
    
        i = 0
        for i in range(max_steps):
            # get the best non expanded node so far
            to_expand = aln.get_best_node_to_expand()

            if to_expand:
                # get the (i, j) position of the node to expand and its score
                coords = to_expand.coords
                score = to_expand.score

                # compare the node to expand with the best already expanded
                # for the same (i, j) position
                if score < scoreboard.get(coords, -math.inf):
                    # if the node is no better just "kill" it, i.e. mark it as expanded
                    # without actually expanding it
                    kill = True
                    #an_expanded.set_tag(color="#000000")
                else:
                    # if the node is the same or better, update the score and expand it
                    scoreboard[coords] = score
                    kill = False
                    
                to_expand.expand(kill=kill)

            if expanded is None:
                # by definition, the best leaf at the end of the expansion must be the solution
                solution = aln.get_solution()
                solution.color = COLOR_SOLUTION_BOX
                break

                
            self._count_steps += 1

        # colour green the latest node to be expanded
        if an_expanded:
            an_expanded.set_tag(color="#00ff00")

        # colour red the best unexplored node so far
        if not self._solved:
            an_start.get_best_non_expanded().set_tag(color="#ff0000")

        self._an_start = an_start