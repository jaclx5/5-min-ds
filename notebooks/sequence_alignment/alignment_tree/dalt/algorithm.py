import math

from .alignment import Alignment

COLOR_EXPANDED_BOX = "#00ff00"
COLOR_SOLUTION_BOX = "#0000ff"
COLOR_BEST_BOX = "#ff0000"
COLOR_IGNORED_BOX = "#000000"
COLOR_BEST_FROM_SET_BOX = "#FF00FF"

class Algorithm:
    """
    Abstract base class for all algorihtms
    """
    def run(self, aln:Alignment, max_steps):
        """
        Updates the alignment `aln` with the state of the algorithm after `max_steps` or until a
        solution is found.

        Returns:
            tuple (bool, int): True if a solution was found, the number of steps executed.
        """
        assert False, 0     # pragma: no cover

