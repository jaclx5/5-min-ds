# Sequence Alignments

This repository contains the companion code for the [How to Compare Text with Sequence Alignment
Algorithms](https://jaclx5.github.io/sequence_alignments_1) blog post series.

The `alignment_tree` folder contains the code used to generate the alignment tree images for the
blog. Inside it you can find:

- The `dalt` (from __d__rawing __al__ignment __t__rees acronym) package that exposes an API to
  draw the trees and explore alignment alogrithms.
- The `demo.ipynb` notebook with the code to generate the tree images for the post.

## `dalt` package

The package is organized in a series of abstraction layers to simplify the drawing process. It
exposes the following classes:
- `Canvas` (`canvas.py`): Takes care of actually producing the images abstracting from the upper
  layers the complexity of the SVG description. The class exposes a simple API with the `add_box`
  and `add_link` methods to position and link the boxes in the canvas, and the method `image` that
  draws whole think on an image.
  When adding boxes to the `Canvas` a column and row must be specified. Columns correspond to the
  depth of the box in the tree and rows the horizontal position of the box.
  Links will allways draw a line between the right most edge of the first box to the left edge of
  the second box.
- `Node` (`node.py`): Represents the nodes of a tree that can be drawn in the `Canvas`. Each new
  `Node` must be added to the tree as a child of an existing `Node`. The `draw` method is called to
   generate the final image of the tree.
- `AlignmentNode` (`alignment.py`): Sub class of `Node`. It represents a step (or partial alignment)
  in the algorithm. It contains the methods that allow the generation of new steps of an alignment
  by expanding the current step.
- `Alignment` (`alignment.py`): Sub class of `AlignmentNode`. It contains some convenience methods
  to traverse the alignment tree and find specific partial alignments (e.g. the best allignment at
  a given stge of the algorithm, the best non explored alignments, ...). This greatly simplifies
  the development of the aligorithms.
- `Algorithm` (`algorithm.py`): Abstract class that must be inherited by all algorithms to be
  tested. It exposes a single method `run` that must be implemented by the concrete algorithms.
  The methods should advance the algorithm up to the specied step.
- `Simulation` (`simulation.py`): Exposes the `frame` and `movie` methods. The former takes an
  alignment and an algorithm and runs it a number of steps returning an image of the final tree.
  The later generates all the frames from step 1 until a predefined number of steps.

## The Algorithms

__It's important to note that the algorithms in this package ARE NOT efficient and ARE NOT intended
to be used in any practical way. They were developed uniquely to illustrate the concepts in this
blog and to generate graphical representations of alignments for pedagogic purposes.__

- `AlgorithmGreedy` (algorithm_greedy.py): Implements the greedy algorithm as described in the
  [Third Post](https://jaclx5.github.io/sequence_alignments_3) of the series. It allways explore
  the best non expanded node and stops as soon as it finds a solution.
  
- `AlgorithmBruteForce` (algorithm_greedy.py): Implements the Brute Force algorithm as described in
  the [Third Post](https://jaclx5.github.io/sequence_alignments_3) of the series. It explores all
  possible alignments, not very practical indeed!

  