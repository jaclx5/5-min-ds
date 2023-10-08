---
title: How to Compare Text with Sequence Alignment Algorithms - Part 2
layout: post
author: jaclx5
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

