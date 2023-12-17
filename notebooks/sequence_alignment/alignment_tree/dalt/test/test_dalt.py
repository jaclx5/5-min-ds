# to get coverage run:
# $ coverage run test_dalt.py; coverage html

import sys

# adding parent folder to the system path
sys.path.insert(0, '../..')
 
from dalt.alignment import Alignment
from dalt.simulation import Simulation
from dalt.algorithm_bf import AlgorithmBruteForce

MATCH = 2
MISMATCH = -1
GAP = -2

#
# First test
#
aln = Alignment("AB", "AC", MATCH, MISMATCH, GAP)

best_leaf = aln.get_best_leaf()
assert aln == best_leaf

best_non_expanded = aln.get_best_node_to_expand()
assert aln == best_non_expanded

solution = aln.get_solution()
assert solution is None

aln.expand()

best_leaf = aln.get_best_leaf()
assert best_leaf.id == "*.2"

best_non_expanded = aln.get_best_node_to_expand()
assert best_non_expanded.id == "*.2"

best_non_expanded.expand()

solution = aln.get_solution()
assert solution.id == "*.2.2"

# try to expand `aln` to test if a reexpansion is avoided
try:
    aln.expand()
except Exception as e:
    pass

aln.draw_text()

# expand the solution (it will try to expand but will not add new nodes)
solution.expand()
assert len(solution._children) == 0

# because we expanded the solution the best non expanded is "*.2.1"
best_non_expanded = aln.get_best_node_to_expand()
assert best_non_expanded.id == "*.2.1"

# mark it as "expanded" without actually expanding it
best_non_expanded.expand(kill=True)

#
# First test
#
# define the alignment sequence and scoring scheme
aln = Alignment("AB", "AX", 2, -1, -1)

# get an instance of the algorithm to use
algo = AlgorithmBruteForce()

# get an instance of the simulation engine
s = Simulation(aln, algo)

s.frame(max_steps=20).save("x.png")

s.movie(max_steps=10).save(".")

print(s.get_steps())

s.movie(max_steps=10, start_step=1).save(".", image_name="other_step_$STEP$.png")