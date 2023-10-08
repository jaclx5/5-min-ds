---
title: How to Compare Text with Sequence Alignment Algorithms - Part 2
layout: post
author: jaclx5

greedy:
    name: greedy
    slides:
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_00.png
            caption: In the first step we explore the three operations expanding from the start "empty" alignment (GREEN box).<BR/> The best obtained alignment (P/P) is represented by the RED box.
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_01.png
            caption:  Expand the (P/P) alignment.
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_02.png
            caption:  Expand the (P/P) alignment.
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_03.png
            caption: Expand the (PO/PU) alignment.
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_04.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_05.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_06.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_07.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_08.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_09.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_10.png
            caption: Hello there 1
        - slide:
            image: /images/sequence_alignments/greedy/greedy_slide_resized_11.png
            caption: Hello there 1
---

[[REREAD THE WHOLE POST AND CHECK:
- IMprove the grredy alignment diagrams, caption and check the consistency with the explanation in the text.

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

_In the previous post I described the original biological context and questions that sequence alignments came to answer. In this post I will one step further into the formal definition of alignments_


> Check also the [companion notebook](add companion address) (more info at the [end of the post](#code)). 

# How to Define an Alignment?

If you had the opportunity to read my [previous post](/sequence_alignments_1), I hope you built a good intuition about what an alignment is and what is it good for. In this post I want to go a step further into analyzing sequence alignments and devising effective algorithms to build them automatically. To achieve this, intuition is important but not enough. We need to formally define sequence alignments and establish a way to quantify how "_good_" they are.

## Let's get formal

If $$A$$ and $$B$$ are two sequences of letters of an alphabet, an alignment between them is a succession of "operations" on $$A$$ that will make it identical to $$B$$.

Let $$A_i$$ the letter in the $$i^{th}$$ position of $$A$$. The operations required to perform an alignment can be:
- __MATCH__: Keep $$A_i$$ unchanged.
- __MISMATCH__: Replace $$A_i$$ by a letter of the alphabet.
- __GAP__, can be either:
    - Delete $$A_i$$, shifting all remaining letters to the left, i.e., $$\forall j : j > i, A_{j-1} \leftarrow A_{j}$$).
    - Insert a letter of the alphabet in position $$i$$, shifting all remaining letters to the right, i.e., $$\forall j : j > i, A_j \rightarrow A_{j+1}$$.

Lets take the two following short sequences as an example:

<pre class="graybox">
<code class="python">
    A = POINTER
    B = PUNTER
</code>
</pre>

> _In the examples bellow we represent MATCHES with a "\|", GAPs with a "-", and MISMATCHES with a "x"._

A possible alignment could be the obtained by performing the following two operations:

- Replace $$A_2$$ by "U"
- Delete $$A_3$$
- Keep all other positions unchanged.

<pre class="graybox">
<code class="python">
    A = POINTER
        |x-||||
    B = PU-NTER
</code>
</pre>

Note that although the alignment above seems somehow "natural", there's nothing special about it. The following alternative sequence of operations produces an equally valid alignment:

- Delete $$A_2$$
- Insert "U" in $$A_2$$
- Delete $$A_3$$
- Keep all other positions unchanged.

<pre class="graybox">
<code class="python">
    A         = PO-INTER
                |---||||
    B         = P-U-NTER
</code>
</pre>

A last extreme example (and quite silly by the way) would be to delete all letters from sequence $$A$$ and inserting each letter of sequence $$B$$ in their positions:

- Delete $$A_1 \times 7$$ (note that, due to the shift, $$A_1$$ iterates over each letter of $$A$$)
- Insert "P" in $$A_1$$
- Insert "U" in $$A_2$$
- Insert "N" in $$A_3$$
- Insert "T" in $$A_4$$
- Insert "E" in $$A_5$$
- Insert "R" in $$A_6$$

<pre class="graybox">
<code class="python">
    A         = POINTER------
                -------------
    B         = -------PUNTER
</code>
</pre>

## Some Alignments Are Better Than Others

From the examples above it is easy to imagine that a huge number of possible alignments between two sequences exists. In fact, for those two small sequences there are [almost 20000 possible alignments](#note1).

With so many possible alignments the question that arises is: __"How to find _the best alignment_ of all?"__

In order to answer this question we must compute a score for each alignment and choose, as the "[best alignment](#note2)", the one with the highest score. A simple way to define a score is to assign a value to each operation that generates the alignment.

As a simple example lets take the first of the above alignments:

<pre class="graybox">
<code class="python">
    A = POINTER
        |x-||||
    B = PU-NTER
</code>
</pre>

This alignment consisted on:
- 5 MATCHEs
- 1 MISSMATCH
- 1 GAP

If we assign an arbitrary value to each operation say: $$MATCH=3$$, $$MISSMATCH=-1$$, $$GAP=-2$$, and sum all these values we obtain a total score for the whole alignment:

$$(5 \times 3) + (1 \times -1) + (1 \times -2) = 12$$

The same way the for the remaining two examples we would have scores of:

<pre class="graybox">
<code class="python">
    A         = PO-INTER
                |---||||
    B         = P-U-NTER
</code>
</pre>

$$(5 \times 3) + (0 \times -1) + (3 \times -2) = 9$$

And:

<pre class="graybox">
<code class="python">
    A         = POINTER------
                -------------
    B         = -------PUNTER
</code>
</pre>

$$(0 \times 3) + (0 \times -1) + (13 \times -2) = -26$$

In this example, the first alignment is clearly the "_best_" of the three. Notice that, in this context, "_best_" means "_the one with a higher score_". 

Note that the score depends on values __arbitrarily__ assigned to each one operation. If we would make a different choice of values, say: $$MATCH=0$$, $$MISSMATCH=0$$, $$GAP=1$$, we would obtain different alignment scores:

- <code>Alignment 1</code>: $$(5 \times 0) + (1 \times 0) + (1 \times 1) = 1$$
- <code>Alignment 2</code>: $$(5 \times 0) + (0 \times 0) + (3 \times 1) = 3$$
- <code>Alignment 3</code>: $$(0 \times 0) + (0 \times 0) + (13 \times 1) = 13$$

With these new values, the 3rd alignment becomes the "_best_", even if it does not look like a good alignment.

The point here is that the concept of "_best alignment_" depends on what we want to do with it. In certain scenarios we would prefer the alignment with the smallest number of changes, which, in our example, would be the first one. In other cases, however, we would want to explicitly avoid certain letters substitutions. Each real world application of the alignment have different goals and, for each of these goals, a specific assignment of values must to be fine tuned.

An actual example is the comparison of biological sequences (DNA or Proteins). Here the goal is to "reconstruct" the evolutionary history of each position of the [original sequence](#note3)  of nucleotides (for DNA) or amino acids (for proteins). Computational Biologists dedicate a lot of effort to infer the most probable changes to DNA or proteins during the evolutionary process (also called "_mutations_"), and also to compute the pairing values that will produce the alignments that best reflect those probabilities. The result is a [__Substitution Matrix__](https://en.wikipedia.org/wiki/Substitution_matrix) in which each pair of nucleotide or amino acids substitutions have its specific value.

## Pause for Recap

A few observations are worth mentioning at this point:

First, the concept of "_best alignment_" is known in the technical literature as **_optimal alignment_**, according to the mathematical meaning of [optimization](https://en.wikipedia.org/wiki/Mathematical_optimization). This is the name we will use from now on.

Second, it's easy to see that playing with the values assigned to each pairing will impact the final alignment score dramatically. So, the choice of values must have in mind the practical objectives of the alignment. In this post I'm caring most about comparing text, so the values of the first example ($$MATCH=3$$, $$MISSMATCH=-1$$, $$GAP=-2$$) that promote matches and penalize a little less mismatches than gaps, tend to work well.

Third, the absolute value of each alignment's score is not important by itself. We are only interested in the alignment with the maximum score. This does not mean that the values are always meaningless. In the case of the [__Substitution Matrixes__](https://en.wikipedia.org/wiki/Substitution_matrix) mentioned above the values assigned to each pairing is related with the [probability of observing the mutation in nature](https://en.wikipedia.org/wiki/Substitution_matrix#Log-odds_matrices)


**** vvv CONTINUE HERE vvv
MAYBE SPLIT HERE THE POSTS AND MAKE THIS THE END OF THE SECOND ONE



# Alignment Algorithms

At this point we know everything about how to define and quantify alignments. We still don't know, though, __how to find the optimal alignment__ score. In the remaining of this post we will explore some algorithms that, given a pair of sequences and the values for each pairing operation, return the **_optimal alignment_** for those two sequences.

## Greedy Approach - A False Start

A straightforward, although _naïve_, way to search for the optimal alignment would be to systematically try the possible alignment operations, "greedily" choosing each next step according to the best score found so far.

The sequence of images bellow illustrates this __greedy__ algorithm on our familiar sequences: <code class="python">POINTER</code> and <code class="python">PUNTER</code>:

{% include slideshow.html slideshow=page.greedy %}

In the first step you can see that we started with an _empty_ alignment with an obvious score of 0 and expanded it by applying the three possible operations:

- Align the next letter from the top sequence adding a GAP to the bottom sequence, obtaining a score of -2.
- Align the next letters from both sequences, obtaining a score of 3 due to the MATCH.
- Align the next letter from the bottom sequence adding a GAP to the top sequence, obtaining a score of -2.

The __green__ box represents the expanded node and the __red__ box the best alignment obtained so far.

Now, we repeat the same process expanding the best node found so far, thus the __greedy__ nature of the algorithm. We keep repeating it until we "exhaust" both sequences.

As it progresses, the algorithm recursively builds a ternary tree containing in its leaves the many alternative alignments.

Notice that, due to the accumulation of gaps and mismatches, in some steps (e.g., step 4, 7), the next best alignment does not result from the one previously expanded. In theses cases the algorithm backtracks to explore other branches of the expansion tree.

As you may already guessed, from the title of this section, although this algorithm seems pretty effective, it is not guarantee to find the best alignment. As it usual with all greedy algorithms.

Check this simple counter example:

$$s1 = ABC$$
$$s2 = ABXABC$$

The greedyness of the algorithm forces it to explore what looks like a promising path, immediately aligning the three first letters of both sequences, which leads to the following sub-optimal alignment:

{% include slideshow.html slideshow=page.greedy2 %}

<pre class="graybox">
<code class="python">
    ABC---
    ||x---
    ABXABC
</code>
</pre>

$$(2 \times 3) + (1 \times -1) + (3 \times -2) = -1$$

It's easy to see in this case that the alternative alignment would be better:

<pre class="graybox">
<code class="python">
    AB---C
    ||---|
    ABXABC
</code>
</pre>

$$(3 \times 3) + (0 \times -1) + (3 \times -2) = 0$$

To find this alignment we would need to back track and explore all options which would be a "__Brute Force__" approach, i.e., test all possibilities to make sure we find the optimal alignment.

Of course, due to the huge number of possible alignments, the brute force approach would be completely impractical for all but the smallest sequences.

<div align="center">
    <img src="/images/sequence_alignments/greedy/full_sample.png" height="300"/>
</div>

Fortunately not everything are bad news. We still can find the optimal alignment without having to explore all possible alignments. With the intuition we acquired exploring the greedy approach and discussing the brute force approach it will be easier to understand the best approach so far to discover optimal alignments: The __Dynamic Programming__ algorithm.


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

__<a name="note1">[1]</a>__ The following exact formula for the total number of alignments between a pair of sequences of sizes $$m$$ and $$n$$ was derived by [Torres et al.](#Torres2004):

$$f(m, n) = \sum^{min(m, n)}_{k=0} {2^k {m \choose k} {n \choose k}}$$

The following Python code computes the number of alignments for two sequences of sizes $$7$$ and $$6$$:

<pre class="graybox">
<code class="python">from math import comb

def f(m, n):
    return sum([(2**k) * comb(m, k) * comb(n, k) for k in range(0, min(m, n))])

print(f(7, 6))

>>> 19377</code>
</pre>

__<a name="note2">[2]</a>__ This definition of "best" may seem a little artificial. I hope it will become clear, later in the post, that there is a more profound motivation for it.


__<a name="note3">[3]</a>__ In this case the original sequence is the _ancestor_ sequence of the ones being compared.

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