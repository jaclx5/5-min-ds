---
title: Using Sequence Alignment Algorithms to Compare Text - Part 1
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

_In the previous post I described the original biological context and questions that sequence alignments came to answer. In this post I will move one step further into the formal definition of alignments, before entering the alignment algorithms themselves_

# How to Define an Alignment?

If you had the opportunity to read my [previous post](/sequence_alignments_1), I hope you built a good intuition about what an alignment is and what is it good for. In this post I want to go a step further into analyzing sequence alignments before starting to discuss effective algorithms to build them automatically.

To achieve our goals, intuition is helps but it is not enough. We need to formally define sequence alignments and establish a way to quantify how "_good_" they are.

# Let's get formal

If $$A$$ and $$B$$ are two sequences of letters of an alphabet, an alignment between them is a succession of "operations" on $$A$$ that will make it identical to $$B$$.

Let $$A_i$$ the letter in the $$i^{th}$$ position of $$A$$. The operations required to perform an alignment can be:
- __MATCH__: Keep $$A_i$$ unchanged.
- __MISMATCH__: Replace $$A_i$$ by anohter letter of the alphabet.
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

A third extreme example (and quite silly by the way) would be to delete all letters from sequence $$A$$ and inserting each letter of sequence $$B$$ in their positions:

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

Note that until this moment we are just concerned in transforming $$A$$ into $$B$$, so, for now, any sequence of operations that achieves this goal is good enough.

# Some Alignments Are Better Than Others - The Optimal Alignment

From the examples above it is easy to imagine that the number of possible alignments between two sequences is huge. In fact, for sequences as small as the ones from our example there are [almost 20000 possible alignments](#note1).

With so many possible alignments it would be really nice to find __"the _best alignment_ of all"__. To define "_best_" we must first compute a score for each alignment, so we can compare them numerically and choose, as the "[best alignment](#note2)", the one with the highest score.

A practical way to define a score is to assign a value to each operation that generates the alignment and sum all these values.

Lets go back to the sequences of our previous example:

<pre class="graybox">
<code class="python">
    A = POINTER
        |x-||||
    B = PU-NTER
</code>
</pre>

This alignment consists on:
- 5 MATCHEs
- 1 MISSMATCH
- 1 GAP

Lets assign a value (totally arbitray in this case) to each operation:
- MATCH = 3
- MISSMATCH = -1
- GAP = -2

Summing all these values we obtain a total score for the whole alignment:

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

With the proposed scoring scheme, the first alignment is clearly the "_best_" of the three. Notice that, in this context, "_best_" means "_the one with a higher score_". 

Note that the score depends on values __arbitrarily__ assigned to each one operation. If we would make a different choice of values for our scoring scheme, say:
- MATCH = 0
- MISSMATCH = 0
- GAP = 1

We would obtain different alignment scores:

- <code>Alignment 1</code>: $$(5 \times 0) + (1 \times 0) + (1 \times 1) = 1$$
- <code>Alignment 2</code>: $$(5 \times 0) + (0 \times 0) + (3 \times 1) = 3$$
- <code>Alignment 3</code>: $$(0 \times 0) + (0 \times 0) + (13 \times 1) = 13$$

With these new values, the 3rd alignment becomes the "_best_", even if it does not look like a good alignment.

The point here is that the concept of "_best alignment_" depends on what we want to do with it, and the scoring scheme must be chosen according to our goals.

In certain scenarios we would prefer the alignment with the smallest number of changes, which, in our example, would be the first one. In other cases, however, we would want to explicitly avoid certain letter substitutions (e.g. we could have different scores for different pairs of matching letters). Each real world application of the alignment have a different goal and, for each goal, a specific assignment of values must be fine tuned.

If all of this seems too abstract, I would take as a real world example the comparison of biological sequences (DNA or Proteins). Here the goal is to "reconstruct" the evolutionary history of each position of the [ancestor sequence](#note3)  of nucleotides (for DNA) or amino acids (for proteins). Computational Biologists dedicate a lot of effort to infer the most probable changes to DNA or proteins during the evolutionary process (also called "_mutations_"), and also to compute the pairing values that will produce the alignments that best reflect those probabilities. The result is a [__Substitution Matrix__](https://en.wikipedia.org/wiki/Substitution_matrix) in which each pair of nucleotide or amino acids substitutions have its specific value.

A common scoring scheme to align coding DNA is the following substitution matrix:

| | A| C| G| T|
|-| -| -| -| -|
|A| 4|-2| 0|-2|
|C|-2| 4|-1| 0|
|G| 0|-2| 4|-2|
|T|-2| 0|-2| 4|

GAP = -1

Os matches são represenatos na diagonal principal
Todas as outras posições correspondem a mismmatches
A matriz é simétrica e os valores representam a probabilidade, observada na natureza de cada nucleotido ser mutado por outro [Hamady 2006](#Hamady2006)

Dar exemplo de duas pequenas sequencias de DNA...

# A Pause for Recap

A few observations are worth mentioning at this point:

First, the concept of "_best alignment_" is known in the technical literature as **_optimal alignment_**, according to the mathematical meaning of [optimization](https://en.wikipedia.org/wiki/Mathematical_optimization). This is the name we will use from now on.

Second, it's easy to see that playing with the values assigned to each pairing will impact the final alignment score dramatically. So, the choice of values must have in mind the practical objectives of the alignment. In this post I'm caring most about comparing text, so the values of the first example ($$MATCH=3$$, $$MISSMATCH=-1$$, $$GAP=-2$$) that promote matches and penalize a little less mismatches than gaps, tend to work well.

Third, the absolute value of each alignment's score is not important by itself. We are only interested in the alignment with the maximum score. This does not mean that the values are always meaningless. In the case of the [__Substitution Matrixes__](https://en.wikipedia.org/wiki/Substitution_matrix) mentioned above the values assigned to each pairing is related with the [probability of observing the mutation in nature](https://en.wikipedia.org/wiki/Substitution_matrix#Log-odds_matrices)


# Notes

__<a name="note1">[1]</a>__ The following exact formula for the total number of alignments between a pair of sequences of sizes $$m$$ and $$n$$ was derived by [Torres et al.](#Torres2004):

$$f(m, n) = \sum^{min(m, n)}_{k=0} {2^k {m \choose k} {n \choose k}}$$

The following Python code computes the number of alignments for two sequences of sizes $$m=7$$ and $$n=6$$:

<pre class="graybox">
<code class="python">from math import comb

def f(m, n):
    return sum([(2**k) * comb(m, k) * comb(n, k) for k in range(0, min(m, n))])

print(f(7, 6))

>>> 19377</code>
</pre>

__<a name="note2">[2]</a>__ This definition of "best" may seem a little artificial. I hope it will become clear, later in the post, that there is a more profound motivation for it.


__<a name="note3">[3]</a>__ The _ancestor_ sequence is the one that originated both sequences being aligned through mutation.

# References

<a name="Torres2004">[Torres2004]</a> Torres A., Cabada A., Nieto J. (2004). __An Exact Formula for the Number of Alignments Between Two DNA Sequences__. DNA sequence, 14, 427-430.

<a name="Hamady2006">[Hamady2006]</a> Hamady M., Betterton M., Knight R. (2006). __Using the nucleotide substitution rate matrix to detect horizontal gene transfer.__ BMC bioinformatics. 7:476.