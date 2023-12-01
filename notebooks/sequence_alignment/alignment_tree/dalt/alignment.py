from enum import Enum, auto

from .node import Node


class Operation(Enum):
    MATCH = 1
    GAP_UP = 2
    GAP_DOWN = 3


class AlignmentNode(Node):
    """
    This class represents a Node in the alignment tree (i.e. a partial alignment step of an
    alignment algorithm)

    Each instance of this class is either a root of a sub tree containing the partial alignments
    that derive from it or the leaf node of the tree containing a complete alignment from each no
    other alignment can be derived.

    This class takes care of:
        - Keeping the information about:
            - The sequences to be aligned.
            - The scoring scheme.
            - The textual representation of the alignment step.
            - The current score of this step.
            - The next position of each sequence to be consumed.
        - Generating its own children alignments.

    See the super class `node.Node` for more information.

    Args:
        seq1, seq2 (str): Sequences to be aligned.
        vmatch, vmismatch, vgap (int): Values of the scoring scheme.
        ops (list[int]): List of operations performed so far.

    Private Attributes:
        _score (int): Score of the current step.
        _i, _j: Next positions of the first and second sequences to be consumed.
    """

    def __init__(self, seq1: str, seq2: str, vmatch: int, vmismatch: int, vgap: int, ops: list[int]=[]):
        super().__init__()

        """
        Computes the score and the textual representation of the alignment
        """
        self._seq1 = seq1
        self._seq2 = seq2
        self._vmatch = vmatch
        self._vmismatch = vmismatch
        self._vgap = vgap
        self._ops = ops

        self._apply_ops()

    def _apply_ops(self):
        # "start" is the name by default for the empty alignment
        mask, mseq1, mseq2 = "start", "", ""
        
        self._score, i, j = 0, 0, 0

        if self._ops:
            # compute mask & score
            mask = ""

            for op in self._ops:
                match op:
                    case Operation.MATCH if self._seq1[i] == self._seq2[j]:
                        ms1, m, ms2, inc_score, inc_i, inc_j = self._seq1[i], "|", self._seq2[j], self._vmatch, 1, 1

                    case Operation.MATCH if self._seq1[i] != self._seq2[j]:
                        ms1, m, ms2, inc_score, inc_i, inc_j = self._seq1[i], "x", self._seq2[j], self._vmismatch, 1, 1

                    case Operation.GAP_UP:
                        ms1, m, ms2, inc_score, inc_i, inc_j = "-", "-", self._seq2[j], self._vgap, 0, 1

                    case Operation.GAP_DOWN:
                        ms1, m, ms2, inc_score, inc_i, inc_j = self._seq1[i], "-", "-", self._vgap, 1, 0

                    case _: # pragma: no cover
                        raise Exception(f"Invalid operation {op}")

                mseq1 += ms1
                mask += m
                mseq2 += ms2
                i += inc_i
                j += inc_j

                self._score += inc_score

        self._i, self._j = i, j 

        """
        Finally builds the textual representation of the AlignmentNode which is a three line string
        containing:
            - The first sequence with the operations applied.
            - The masked representation of the operations + the score
            - The second sequence with the operations applied.
        """
        self.text = f"{mseq1}\n{mask} ({self.score})\n{mseq2}"

    def _get_score(self):
        return self._score

    score = property(fget=_get_score, doc="Computed score of the alignment.")

    def _can_consume_seq1(self):
        return self._i < len(self._seq1)

    def _can_consume_seq2(self):
        return self._j < len(self._seq2)

    def _can_consume(self):
        return self._can_consume_seq1() and self._can_consume_seq2()

    def child_alignment_factory(self, op):
        if op == Operation.GAP_DOWN and self._can_consume_seq1() or \
           op == Operation.GAP_UP and self._can_consume_seq2() or \
           op == Operation.MATCH and self._can_consume():
            return Alignment(self._seq1, self._seq2, self._vmatch, self._vmismatch, self._vgap, self._ops + [op])
        else:
            return None

    def is_solution(self):
        return not self._can_consume_seq1() and not self._can_consume_seq2()


class Alignment(AlignmentNode):
    """
    This class extends the AlignmentNode class with some convenience methods for algorithm
    implementation.

    The user of this class usually creates a single instance which will be the root of the full
    alignment tree and then uses the available methods to control the algorithm.

    Args:
        see `AlignmentNode` class.

    Private Attributes:
        see `AlignmentNode` class.
        _expanded (bool): Indicates if the node was already expanded.

    """
    def __init__(self, seq1: str, seq2: str, vmatch: int, vmismatch: int, vgap: int, ops: list[int]=[]):
        super().__init__(seq1, seq2, vmatch, vmismatch, vgap, ops)

        self._expanded = False

    def reset(self):
        super().reset()
        
        self._expanded = False

    def expand(self, kill=False):
        if self._expanded:
            assert False, f"Trying to re-expand an already expanded node!\n----\n{self.text}\n-----"

        self._expanded = True

        if kill:
            return

        for op in (Operation.GAP_DOWN, Operation.MATCH, Operation.GAP_UP):
            child = self.child_alignment_factory(op)

            if child is not None:
                self.add_child(child)

    def get_node_to_expand(self):
        """
        Returns the first child that can be expanded.

        Returns:
            Alignment: the alignment corresponding to the children to be expanded.
        """
        if not self._expanded and not self.is_solution():
            return self
        else:
            for child in self._children:
                to_expand = child.get_node_to_expand()

                if to_expand is not None:
                    return to_expand

        return None
            
    def get_best_node_to_expand(self):
        """
        Returns the child with the highest score, from all children of the sub tree that can be expanded.

        Returns:
            Alignment: the alignment corresponding to the best children.
        """
        if self._expanded:
            best = list(filter(lambda _: _, [child.get_best_node_to_expand() for child in self._children]))

            if best:
                return max(best, key=lambda child: child.score)
            else:
                return None

        elif self.is_solution():
            # solutions nodes cannot be expanded
            return None

        else:
            return self

    def get_best_leaf(self):
        """
        Returns the leaf alignment node with the highest score, from all leafs of the sub tree.

        Returns:
            Alignment: the alignment corresponding to the best leaf. 
        """
        if self._children:
            return max([child.get_best_leaf() for child in self._children], key=lambda child: child.score)
        else:
            return self

    def get_solution(self):
        """
        Returns the first solution node starting in a width-first search to the tree.
        """
        if self.is_solution():
            return self
        else:
            for child in self._children:
                solution = child.get_solution()

                if solution:
                    return solution

        return None
