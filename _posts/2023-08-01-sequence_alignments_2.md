---
title: How to Compare Text with Sequence Alignment Algorithms - Part 2
layout: post
author: jaclx5
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

<style>
    div#DNA {
        font-family:monospace;
        font-size: 20px;
    }
    .seqname {
        color: black;
    }
    .A {
        color: red;
    }
    .C {
        color: green;
    }
    .G {
        color: orange;
    }
    .T { 
        color: blue;
    }
</style>

_In the previous post of this series I described the biological context and appplications of sequence alignment algorithms as a motivation for it's development. In this second post I will go deep into the technical definitions and implementation details.<<rephrase>>_

In the first post of this series I discussed the original biological problems that sequence alignments came to solve. Now we will go deeper into the technical description of alignments.

> Check also the [companion notebook](add companion address) (more info at the [end of the post](#code)). 

I this second post I will focus on "Optimal Sequence Alignments". But first lets start by formally defining Sequence Alignment.


# What is an Alignment?

An alignment between two sequences of letters $$A$$ and $$B$$ is a succession of "changes" on $$A$$ that will transform it on $$B$$. To simplify the analysis of alignment algorithms we consider only atomic "changes", i.e., the possible changes can only be one of the following ($$A_i$$ corresponds to the letter in the $$i^{th}$$ position of $$A$$):
- Replace any element $$A_i$$ by any other letter of the vocabulary.
- Delete any element  $$A_j$$.
- Insert an element after or before position $$k$$.

As a simple example lets take two very short sequences "POINTER" and "PUNTER".

[[change the background to make it visually nicer]]
[[check here for inspiration: https://medium.com/@clem.boin/creating-a-minimal-kernel-development-setup-using-qemu-and-archlinux-987896954d84]]

<pre>
<code class="python">
    A = POINTER
    B = PUNTER
</code>
</pre>

A possible alignment could be the obtained by performing the following two changes:

- Replace $$A_2$$ by "U"
- Delete $$A_3$$

The resulting alignment would be:

<pre>
<code class="python">
    A = POINTER
         x-
    B = PU-NTER
</code>
</pre>

Note that although the alignment above seems somehow "natural", there's nothing special about it. The following alternative list of changes:

- Delete $$A_2$$
- Insert an "U" before position $$2$$
- Delete $$A_3$$

Produce an equally valid alignment:

<pre>
<code class="python">
    A         = PO-INTER
                 ---
    B         = P-U-NTER
</code>
</pre>

In the extreme (and clearly silly) case one could consider changing all letters:

- Delete $$A_1$$ to $$A_7$$
- Insert an "P" in position $$1$$
- Insert an "U" after position $$1$$
- Insert an "N" after position $$2$$
- Insert an "T" after position $$3$$
- Insert an "E" after position $$4$$
- Insert an "R" after position $$5$$

Which would produce an also valid (although not very useful as we will see in a moment) alignment:

<pre>
<code class="python">
    A         = POINTER------
                -------------
    B         = -------PUNTER
</code>
</pre>

_In the examples above we added the symbol «-» in the position of the deleted/inserted letters to facilitate the comparison. Those «-», of course, are not part of the sequence._

From the examples bellow we can easily intuit that there are a huge number of possible alignments between two sequences, in fact for those tiny sequences we can expect almost 20000 possible alignments ([check here for an exact formula](end of the post)). The question that arises is how to choose from all those possible alignments the _best_ one?

There is not a single answer to the question. The _best_ alignment depends on what do we want to achieve with it. Intuitively we would prefer the alignment resulting from the smallest number of changes, which in our example would be the first proposed alignment. However, real world applications have specific goals and/or restrictions. Sometimes we aim to minimize insertions or deletions. In other cases we want to avoid specific letters substitutions, in the protein alignment field a lot of work has been done to determine which substitutions are _better_ than others ([see here](https://en.wikipedia.org/wiki/Substitution_matrix)).

The discussion on how to choose the best alignment lead us to the definition of __optimal alignment__


# Optimal Alignment

In order to choose the optimal alignment from the set of all possible alignments we need find a way to assign some kind of score to each alignment so we can compare them all against each other. A simple way to define a score is to evaluate each position of the alignment and assign some value depending on the aligned letters in that position.

lets assume that $$A_i$$ and $$B_i$$ represent the $i^{th}$ position of the __aligned__ sequences $$A$$ and $$B$$. It's easy to see that we have three possible pairing situations:

- A __MATCH__  if positions $$A_i$$ and $$B_i$$ are identical (represented by the _pipe_ character "|").
- A __MISSMATCH__ if positions $$A_i$$ and $$B_i$$ are different (represented by an "x").
- A __GAP__ if either $$A_i$$ or $$B_i$$ contain a "-" after the alignment (conveniently represented by an "-").

As a simple example lets take the first alignment of the example above:

<pre>
<code class="python">
    A = POINTER
        |x-||||
    B = PU-NTER
</code>
</pre>

This alignment consisted on:
- 5 MATCHEs in positions 1, 4, 5, 6 and 7.
- 1 MISSMATCH in position 2
- 1 GAP in position 3

If now we assign an individual and arbitrary value to each pairing say: MATCH=3, MISSMATCH=-1, GAP=-2 and sum all the pairings we obtain a total score for the whole alignment.

The alignment would then score $$(5 \times 3) + (1 \times -1) + (1 \times -2) = 12$$

The same way the for the remaining two examples we would have scores of $$(5 \times 3) + (3 \times -2) = 9$$:

<pre>
<code class="python">
    A         = PO-INTER
                |---||||
    B         = P-U-NTER
</code>
</pre>

And $$(13 \times -2) = -26$$:

<pre>
<code class="python">
    A         = POINTER------
                -------------
    B         = -------PUNTER
</code>
</pre>

In this example the first alignment is clearly the best one of the three (as we would probably expect).

A few observations are worth mentioning at this point.

First, it's easy to see that playing with the values assigned to each pair we can change the final scores dramatically. So the choice of values must have in mind the practical objectives of the alignment. When comparing biological sequences (DNA or Proteins) the goal is to "reconstruct" the evolutionary history of each position of the sequences (nucleotide or amino acid), thus a lot of research was dedicated to infer those values (see the [substitution matrices](link to Wikipedia) for more information). In our case, when comparing text values that promote matches and penalize more or less equally mismatches and gaps tend to work well.

Second, the absolute value of each alignment's score is not important by itself. We are only interested in the maximum score of all alignments. That does not mean that the values as meaning less. Coming back to biological sequences, the values assigned to each match/mismatch are [[check in the refernce here - the log prob of a mutation....]]

Finally, and more important for our discussion here, we still don't know how to __find__ the best score. [[end this sentence and motivate the following]]

# Alignment Algorithms

## Brute Force Approach

Now that we know how to compute a numerical score for any assignment e need to find a way to search for the best overall alignment, a.k.a. the optimal alignment.

Lets start by illustrating a brute force approach to optimal alignment search:

[[graphical element]]

As we can see in each step of the process we have three potential actions:

- Align the next letter from both sequences (either with a MATCH or MISMATCH).
- Align the next letter from $$A$$ sequence adding a GAP to sequence $$B$$.
- Align the next letter from $$B$$ sequence adding a GAP to sequence $$A$$.

After each of the steps the process repeats until one of the sequences is exhausted and GAPs are added to the remaining letters of the other sequence.

As it progresses, the algorithm recursively builds a ternary tree containing in its leaves the complete alignments. When no more expansion are possible we just have to look for the maximum score alignment, we are done!

Of course, due to the huge number of possible alignments, the brute force approach is completely impractical for all but the smallest sequences. However it provides us with a intuition for the dynamic programming algorithm.

## Dynamic Programming

Dynamic Programming was proposed in the 1950s by Richard Bellman (the choice of the name "Dynamic Programming" has a curious motivation as told by [Bellman himself in his autobiography](#Bellman1984)), it is a general optimization method applicable to a large class of problems well beyond sequence alignment. It's an elegant and useful method which description is out of the scope of this post but which is worth to [explore further](https://en.wikipedia.org/wiki/Dynamic_programming).

To get to the intuition behind dynamic programming let change slightly or brute force algorithm. Instead of expanding all possible combinations of actions we will expand only the most "promising" ones:

[[create a chart with the order of exploration expanding only the best nodes at each moment]]

This way we expand the ternary tree by "layers", expanding always the highest score node until both sequences have been completely aligned. This will avoid the expansion of numerous useless nodes.

Note that the main intuition of the dynamic programming algorithm is that once we obtain the optimal alignment of $$i$$ letters from sequence $$A$$ and $$j$$ letters from sequence $$B$$ we only have to care with the remaining sequences $$A_{k>i}$$ and $$B_{l>j}$$. So in each step $$i, j$$ we reduce our problem to the positions after this step.

To understand a little better lets focus on the first letters of both sequences. There are 3 ways to align $$A_1$$ and $$B_1$$ them:

A-
-A

A
A

-A
A-

Clearly, the second one is the best way to align those letters (which would not be the case of they were different and the GAP score was higher that the MISMATCH score!). Once we settle with this alignment we just have to care about aligning $$A_{i>1}$$ and $$B_{i>1}$$.

Formally we express this with the following equation:

[[add the dynamic programming bracket]]

## Needleman-Wunsch Algorithm

[[ate the end briefly refer their relation with the Levenshtein distance]]

The first proposed application of dynamic programming to pairwise sequence alignments is due to Needleman and Wunsch in 1970 with the, appropriately named, Needleman-Wunsch global alignment algorithm.

The cleverness of the [[NW]] algorithm is that the layer-by-layer approach for developing the alignment tree was replaced by a matrix that is filled left-to-right, top-to-bottom with the best score at any given moment:

[[show the matrix]]

Formally the value of each position in the [[NW]] matrix is defined by the following formula:

[[add a formula]]

When filling up the [[NW]] matrix, we need to keep track of both the score in each position but also the position from which we choose to come so we can backtrack and rebuild the alignment from the bottom-right corner of the matrix. A different backtrack path will generate a different alignment:

[[Show two matrices with one different backtrack and the respective resulting alignments]]

It's easy to understand now how such a small matrix can encode the enormous number of possible alignments:

[[Show graphically how different nodes of the tree represent different paths and each path represents a specific alignment]]


The [[NW]] algorithm is very easy to implement (as we will see in the next post) and pretty effective for relatively short sequences [[provide some benchmark here]].

The reason why the [[NW]] algorithm is called a Global Alignment algorithm has to do with the fact that it is meant to align sequences of similar lengths. If we try to align a shorter sequence against a longer one the results are not very good as we can see in the example bellow:

<pre>
<code class="python">
this is a longer sequence that includes a shorter sequence inside it
---|-----------------------|----------------|||--|||||||||----------
---s-----------------------h----------------ort-- sequence----------
</code>
</pre>

It's easy to see what happens here. As the algorithm explores all the possible paths there is no preference for aligning the letters from the shorter sequence close together as from a score point of view there's no difference in matching a letter from the shorter sequence with any similar letter from the longer one. In this example the "s" and the "h" are aligned with the first occurring "s" and "h" of the longer sequence.

Some changes to the [[NW]] algorithm have been proposed to handle this situation, for example having an extra penalty for starting a gap which would force the algorithm to group the shorter sequence letters. The best solution, however, is to apply a local alignment algorithm.

## Smith-Waterman Algorithm

In 1981 by Smith and Waterman published a local alignment algorithm named, creatively, Smith-Waterman algorithm which is an adaptation of the [[NW]] algorithm designed specifically to capture the best alignment of a small sequence inside a longer one.

The [[SW]] algorithm is based on the same matrix building mechanism of the [NW]] with two major differences:

1. When adding up the values of MATCHES, MISMATCHES and GAPS the final score of a given position never falls bellow zero.
2. The backtracking operation starts from the maximum position of the matrix (instead of the bottom-rightmost position) and stops ate the first time it encounters a zero (instead of going all the way to the top-leftmost position)

The following example illustrates the [[SW]] matrix building mechanism:

[[add a SW matrix]]

Formally the value of each position in the [[SW]] matrix is defined by the following formula:

[[add a formula]]

Now if we apply the [[SW]] to the previous two sequences we obtain a nice local alignment:

[[validate!!]]
<pre>
<code class="python">
this is a longer sequence that includes a shorter sequence inside it
------------------------------------------|||||--|||||||||----------
------------------------------------------short-- sequence----------
</code>
</pre>

# Caveats of the Optimal Alignment

[[mention the $$O(m x n)$$ complexity of the algorithm which make impractical for aligning too large sequences]]

That's why we need alternative algorithms to align really big sequences.

In the next post we will show how to implement a Python version of both the NW and SW algorithm and after that we will explore alternatives to Optimal Alignment.

For now we are done!

# End Notes

## Exact formula for the number of pair wise sequence alignments

In [this article](#Torres2004) the authors propose the following formula for the total number of possible alignments between a pair of sequences of sizes $$m$$ and $$n$$:

$$f(m, n) = \sum^{min(m, n)}_{k=0} {2^k {m \choose k} {n \choose k}}$$

<pre>
<code class="python">
    from math import comb
    
    def f(m, n):
        return sum([(2**k) * comb(m, k) * comb(n, k) for k in range(0, min(m, n))])

    print(f(7, 6))
</code>
</pre>



# References


<a name="Torres2004">[Torres2004]</a> Torres A., Cabada A., Nieto J. (2004). __An Exact Formula for the Number of Alignments Between Two DNA Sequences__. DNA sequence, 14, 427-430.

<a name="Bellman1984">[Bellman1984]</a>Bellman RE. (1984) Eye of the Hurricane: An Autobiography. World Scientific, Singapore.

Needleman, Saul B. & Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. 48 (3): 443–53. doi:10.1016/0022-2836(70)90057-4. PMID 5420325.

Smith, Temple F. & Waterman, Michael S. (1981). "Identification of Common Molecular Subsequences" (PDF). Journal of Molecular Biology. 147 (1): 195–197. CiteSeerX 10.1.1.63.2897. doi:10.1016/0022-2836(81)90087-5. PMID 7265238.

[ref] Eddy, S., "What is dynamic programming?", 2004, Nature Biotechnology, 22 (7), 909-910


==== THIRD POST ===

- **_Dot Plot_ methods (graphic and MUMMER)**

<a name="code"/>
# Links and Software

You can reproduce all the results from this post using the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/wordle/wordle.ipynb).

I also made available the package [pywordlesover](https://github.com/jaclx5/pywordlesolver) which allows you to solve WORDLE puzzles, to play WORDLE on the command line and benchmark the different solution strategies discussed in this post. The code is free to download and play with the code.

If you have questions or ideas feel free to share them in the comments box. 


https://mummer.sourceforge.net/

==== FOURTH POST ===


- **_K-word_ methods (FASTA and BLAST)**
==== FIFTH POST ===

# Beyond Biological Sequences

Just before finishing, one word regarding the application of sequence alignment algorithms to non Biological sequences. It is pretty easy to see that if we replace DNA nucleotides by the letters of any alphabet (or by any tokens) these algorithms can be applied as well.

Of course a few choices need to be made if we want to apply SA algorithms to general sequences (as we will see in the future) but the algorithms can be applied with no changes.
check - for future posts:
text_algorithms-crochemore-1994.pdf
diff algorithms