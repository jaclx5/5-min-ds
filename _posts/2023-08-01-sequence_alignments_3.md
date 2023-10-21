---
title: How to Compare Text with Sequence Alignment Algorithms - Part 3
layout: post
author: jaclx5

greedy1:
    slides:
        - slide:
            image: /images/sequence_alignments/greedy1/step_00.png
            caption: The algorithm starts with an "empty" alignment, which happens to be the trivially best so far.

        - slide:
            image: /images/sequence_alignments/greedy1/step_01.png
            caption: It expands the "empty" alignment (GREEN box), by applying the three possible operations<br/>
                The (P/P) alignment, RED box in the middle, is the best one so far, and will be expanded in the next step.
        - slide:
            image: /images/sequence_alignments/greedy1/step_02.png
            caption: After expanding (P/P), all of the obtained partial alignments have smaller scores.
        - slide:
            image: /images/sequence_alignments/greedy1/step_03.png
            caption: After expanding (PO/PU) we get three partial alignments with score == 1.<br/>
                     The algorithm chooses (arbitrarily) the one from the top, in a "depth first"-like search.
        - slide:
            image: /images/sequence_alignments/greedy1/step_04.png
            caption: As it gets no improvement, it continues expanding the alignments with score of 1...
        - slide:
            image: /images/sequence_alignments/greedy1/step_05.png
            caption: ...and it does it again.
        - slide:
            image: /images/sequence_alignments/greedy1/step_06.png
            caption: Still no luck, the best partial alignment so far has a score == 0.
        - slide:
            image: /images/sequence_alignments/greedy1/step_07.png
            caption: Due to the matching "N"s, the best partial alignment at this point has a score of 3...
        - slide:
            image: /images/sequence_alignments/greedy1/step_08.png
            caption: ...which keeps increasing with the subsequent matches...
        - slide:
            image: /images/sequence_alignments/greedy1/step_09.png
            caption: ...and increases again...
        - slide:
            image: /images/sequence_alignments/greedy1/step_10.png
            caption: ...until it "consumes" both full sequences, obtaining the final alignment in the BLUE box (POINTER/P-UNTER), which is our solution.

greedy2:
    name: greedy
    slides:
        - slide:
            image: /images/sequence_alignments/greedy2/step_00.png
            caption: Again we start with the empty alignment.
        - slide:
            image: /images/sequence_alignments/greedy2/step_01.png
            caption: The A/A alignmnt is clearly the choice at this point.
        - slide:
            image: /images/sequence_alignments/greedy2/step_02.png
            caption: Likewise, matching the "B"s greatly increases the score of the best partial alignment.
        - slide:
            image: /images/sequence_alignments/greedy2/step_03.png
            caption: Although the C/X mismatch decreases the overall score, ABC/ABX is still the best alignment so far.
        - slide:
            image: /images/sequence_alignments/greedy2/step_04.png
            caption: Now, the algorithm backtracks due to the accumulation of gaps and mismatches.
        - slide:
            image: /images/sequence_alignments/greedy2/step_05.png
            caption: Here it will try out the second alignment with a score of 4...
        - slide:
            image: /images/sequence_alignments/greedy2/step_06.png
            caption: ...with not much luck as the overall best score continues to decrease...
        - slide:
            image: /images/sequence_alignments/greedy2/step_07.png
            caption: ...and decrease...
        - slide:
            image: /images/sequence_alignments/greedy2/step_08.png
            caption: ...and decrease...
        - slide:
            image: /images/sequence_alignments/greedy2/step_09.png
            caption: ...and decrease...
        - slide:
            image: /images/sequence_alignments/greedy2/step_10.png
            caption: ...and decrease...
        - slide:
            image: /images/sequence_alignments/greedy2/step_11.png
            caption: ...and decrease.
        - slide:
            image: /images/sequence_alignments/greedy2/step_12.png
            caption: Finally, due to the long sequences of gaps<br /> ABC--/ABXAB becomes the best alignment.
        - slide:
            image: /images/sequence_alignments/greedy2/step_13.png
            caption: Which leads to a premature solution.
---

[[REREAD THE WHOLE POST AND CHECK:
- Improve the grredy alignment diagrams, caption and check the consistency with the explanation in the text.

- Consistent use of pronouns:
    "I" - if my option when writing the post.
    "we" - if it the community the subject.
    "you" - if I'm directly talking to the reader

- Consistent naming of the concept.
]]

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

_In the last post I formally defined what is a sequence alignment, and described a way to quantify how "good" an alignment is. Now I will start devising the algorithms to find the best (a.k.a. optimal) alignment between two sequences._

> Check also the [companion notebook](add companion address) (more info at the [end of the post](#code)). 

# Alignment Algorithms

If you have been following the previous posts of this series I hope you have learned quite bit about alignments. Lets recap a few things we learned so far:
- An alignment between two sequences $$A$$ and $$B$$ is a sequence of operations over $$A$$ that turn it into $$B$$;
- A huge number of distinct alignments exist between $$A$$ and $$B$$;
- It is possible to assign a score of how "good" each of those alignments is;
- And, finally, the alignment (or alignments) with the highest score is called __optimal alignment__.

Now, our next challenge will be to __design an algorithm that finds the optimal alignment between any two given sequences__.



## Greedy Approach - A False Start

Our first attemp will be a straightforward, although _naïve_, way to search for the optimal alignment, by probing systematically all possible alignment operations in each step of the way, "greedily" choosing the next search direction according to the best score found so far.

The sequence of images bellow illustrates this __Greedy__ algorithm on our familiar sequences: <code class="python">POINTER</code> and <code class="python">PUNTER</code>.

In each step of the algorithm the __<span style="color:#ff0000;">red</span>__ box represents the best alignment obtained so far, the __<span style="color:#00ff00;">green</span>__ box represents the previously expanded node, and the __<span style="color:#0000ff;">blue</span>__ one represents the solution when we find it. You can click on the left and right arrows, or in the bottom circles to navigate the steps of the algorithm.


{% include slideshow.html slideshow=page.greedy1 %}

In the first step, the algorithm starts with an _empty_ alignment with an score of 0 (zero), and expands it by applying the three possible operations:

- Align the next letter from the top sequence adding a GAP to the bottom sequence, obtaining a score of -2.
- Align the next letters from both sequences, obtaining a score of 3 due to the MATCH.
- Align the next letter from the bottom sequence adding a GAP to the top sequence, obtaining a score of -2.

Next, it chooses the best node found so far (thus the __greedy__ in the algorithm's name) and repeats the same expansion process. The algorithm keeps expanding the best node, until it "exhausts" both sequences.

As it progresses, the algorithm recursively builds a ternary tree that contains in its leaves the many alternative alignments.

Notice that, due to the accumulation of gaps and mismatches, in some steps (e.g., step 4, 7), the next best alignment doesn't always result from the one previously expanded. In these cases the algorithm "backtracks" to explore other branches of the expansion tree.

Although this algorithm seems pretty effective at a first glance it has a major flaw. As it usualy occurs with Greedy algorithms, it does not allways find the best alignment.

This pair of sequences, <code class="python">ABC</code> and <code class="python">ABXABC</code>, is an example for which the Greedy algorithm fails:

{% include slideshow.html slideshow=page.greedy2 %}

As we can see in this example, the greedyness of the algorithm, and the fact that it stops as soon as it founds a possible solution, leads to a sub-optimal alignment:

<pre class="graybox">
<code class="python">
    ABC---
    ||x---
    ABXABC
</code>
</pre>

$$(2 \times 3) + (1 \times -1) + (3 \times -2) = -1$$

We know that the above alignment is not the optimal one because we can, with a little manual effort, find an alignment with a higher score:

<pre class="graybox">
<code class="python">
    AB---C
    ||---|
    ABXABC
</code>
</pre>

$$(3 \times 3) + (0 \times -1) + (3 \times -2) = 3$$

At least in the sequence alignment world, greed doesn't pay.

## Brute Force - Almost There

To overcame the limitation of the Greedy algorithm one would need to continue exploring the alignment space after having found a first solution, back tracking and exploring all remaining possibilities. This systematic search over all possibilities is commonly called a "__Brute Force__" strategy, i.e., to test all possibilities making sure the optimal alignment is never missed.

Mechanically, the both algorithms are pretty similar, the main difference between the Greedy and the Brute Force algorithm is the stopping condition. While the former stops at the first solution, the latter keeps searching until all solutions have been tested.

Of course, due to the huge number of possible alignments, the brute force approach is completely impractical for all but the smallest sequences.

The bravest readers may want to take a look, at their own risk, at the complete brute force expansion of the alignment between <code class="python">ABC</code> and <code class="python">ABXABC</code>:

<div align="center">
    <a href="/images/sequence_alignments/full_brute_force.png"><img src="/images/sequence_alignments/full_brute_force.png" height="500"/></a>
</div>

If you search patiently, you will spot the optimal alignment in a branch near the center of the image:

<div align="center">
    <img src="/images/sequence_alignments/full_brute_force_detail.png"/>
</div>


At his point you may be wondering: "_Why am I wasting your time with flawed and impractical algorithms?_", I have a couple of good reasons for that:
 
First, a simple algorithm is a good proof of concept to show that an algorithm to find optimal alignments indeed exists (which is not always the case for all problems), even if this simple algorithm is inefficient and impractical.

Second, and most important, understanding how those basic algorithms work allows us to gain a familiarity with the problem that will make it much easier to understand the real useful algorithms that we will discuss later.

Fortunately, we still can find the optimal alignment without having to explore all possible alignments. That's what we will see in the next section.

## Dynamic Programming

Dynamic Programming was proposed in the 1950s by Richard Bellman (the choice of the name "Dynamic Programming" has a curious motivation as told by [Bellman himself in his autobiography](#Bellman1984)), it is a general optimization method applicable to a large class of problems well beyond sequence alignment, really worth to be [further explored](https://en.wikipedia.org/wiki/Dynamic_programming).

To get the fundamental concept behind the __dynamic programming__ lets slightly change or brute force algorithm. Instead of expanding all possible combinations of actions we will expand only the most "promising" ones:

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






# Notes

__<a name="note1">[1]</a>__ Add some note here


# References

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
https://github.com/ashvardanian/Stringzilla?utm_source=tldrnewsletter

==== FOURTH POST ===


- **_K-word_ methods (FASTA and BLAST)**
==== FIFTH POST ===

# Beyond Biological Sequences

Just before finishing, one word regarding the application of sequence alignment algorithms to non Biological sequences. It is pretty easy to see that if we replace DNA nucleotides by the letters of any alphabet (or by any tokens) these algorithms can be applied as well.

Of course a few choices need to be made if we want to apply SA algorithms to general sequences (as we will see in the future) but the algorithms can be applied with no changes.
check - for future posts:
text_algorithms-crochemore-1994.pdf
diff algorithms