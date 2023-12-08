---
title: Using Sequence Alignment Algorithms to Compare Text - Some Formal Definitions - Part 2
layout: post
author: jaclx5
---

_In [previous post](/sequence_alignments_1) I described the original biological context and questions that sequence alignments came to answer. This time I will move one step further into the __formal definition of alignments__. We need this formal definition before discussing the alignment algorithms themselves._

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

# How to Define an Alignment?

If you had the opportunity to read my [previous post](/sequence_alignments_1), I hope you got some intuition about what an alignment is and what is it good for. Intuition alone, however, is not enough to achieve our goals of understanding sequence alignment algorithms and being able to use them. In this post I will move beyond the previous hand waving description of alignments to establish two concepts:
1. A formal definition of sequence alignment;
1. A way to quantify how "_good_" each alignment is.


# Let's get formal

Here is our formal definition of alignment:

> If $$A$$ and $$B$$ are two sequences of letters of an alphabet, an alignment between them is: __a succession of "operations" on $$A$$ that will make it identical to $$B$$__.

Let $$A_i$$ and $$A'_i$$ represent the letters in the $$i^{th}$$ position of sequence $$A$$, before and after an operation, respectively. We will consider the following four operations on $$A_i$$:
- __MATCH__, which keeps $$A_i$$ unchanged, i.e., $$A_i = A'_i$$.
- __MISMATCH__, which replaces $$A_i$$ by anohter letter of the alphabet, i.e., $$A_i \neq A'_i$$.
- __INSERTION__, of a letter in $$A_i$$, shifting all remaining letters to the right, i.e., $$\forall j : j >= i, A_j = A'_{j+1}$$.
- __DELETION__, of $$A_i$$, shifting all remaining letters to the left, i.e., $$\forall j : j >= i, A'_{j} = A_{j+1}$$.

Note that both __INSERTIONS__ and __DELETIONS__ will generically referred as __GAPS__.

## Practical Examples

To make things less abstract we will take the two following short sequences:

{% highlight markdown %}
A = POINTER
B = PUNTER
{% endhighlight %}

and try a few ways to align them, i.e., to transform $$POINTER$$ into $$PUNTER$$ by applying the operations described above.

> From now on, in all examples, I'll "annotate" __MATCHES__ with a "\|", __MISMATCHES__ with a "x", and __GAPS__ with a "\-".

## Alignment 1

A possible alignment can be the obtained by performing the following seven operations:

- a __MATCH__ in $$A_1$$
- a __MISMATCH__ in $$A_2$$ replacing "O" $$\rightarrow$$ "U"
- a __DELETION__ of $$A_3$$
- Four __MATCHES__ in positions $$A_3$$ to $$A_6$$ (note that the index is updated after the deletion, that's why we refer to $$A_3$$ twice)

{% highlight markdown %}
A = POINTER
    |x-||||
B = PU-NTER
{% endhighlight %}

That's it! From this alignment we can see that the two words share the same letter $$P$$ and the suffix $$NTER$$ and only differ in two letters.

This doesn't look like much, in fact, aligning tiny sequences is pretty obvious. Don't give up yet, this is just a first step to help familiarize ourselves with the basic concepts before proceeding to longer and difficult to align sequences.

## Alignment 2

Note that although the alignment in the previous example seems somehow "natural", there is nothing particularly special about it. The following alternative sequence of operations produces an equally valid alignment:

- a __MATCH__ in position $$A_1$$
- a __DELETION__ of $$A_2$$
- an __INSERTION__ of "U" in $$A_2$$
- a __DELETION__ of $$A_3$$
- Four __MATCHES__ in positions $$A_3$$ to $$A_6$$

{% highlight markdown %}
A = PO-INTER
    |---||||
B = P-U-NTER
{% endhighlight %}

This is a very similar alignment where we avoid a __MISMATCH__ on position $$A_2$$.

## Alignment 3

An extreme example would be to delete all letters from sequence $$A$$ and insert each letter of sequence $$B$$ in their respective positions:

- __DELETION__ of $$A_1 \times 7$$ (again, the $$A_1$$ position iterates over all letters of $$A$$)
- a __INSERTION__ of "P" in $$A_1$$
- a __INSERTION__ of "U" in $$A_2$$
- a __INSERTION__ of "N" in $$A_3$$
- a __INSERTION__ of "T" in $$A_4$$
- a __INSERTION__ of "E" in $$A_5$$
- a __INSERTION__ of "R" in $$A_6$$

{% highlight markdown %}
A = POINTER------
    -------------
B = -------PUNTER
{% endhighlight %}

You can argue that this is a silly alignment and I would agree. However, for now, we are just concerned in transforming $$A$$ into $$B$$. Any sequence of operations that achieves this goal is a valid alignment, as good as any other. __For the moment we have no way to say if an alignment is better than the other__.

# The Optimal Alignment - Some Alignments are Better

From the three examples above it is easy to see that the number of possible alignments between two sequences is huge. In fact, for sequences as small as $$POINTER$$ and $$PUNTER$$ there are almost 20000 possible alignments <a name="1eton"/>[[1]](#note1).

With so many possible alignments it would be nice to have a way to choose the "__best__" alignment of all, a.k.a., the "__optimal__" alignment. One way to do this would be to define some numerical quantity for each alignment, that we shall call "__score__", and call "optimal" the one with the highest score.

A practical way to define the alignment's score is to assign a value to each of its operations and sum all these values.

As an example, let's assign one (totally arbitrary) value to each operation:

<table class="data-table">
    <tr><td><b>Operation</b></td><td><b>Value</b></td></tr>
    <tr><td>MATCH</td><td>3</td></tr>
    <tr><td>MISSMATCH</td><td>-1</td></tr>
    <tr><td>GAP</td><td>-2</td></tr>
</table>

Going back to our previous alignment 1:

{% highlight markdown %}
A = POINTER
    |x-||||
B = PU-NTER
{% endhighlight %}

This alignment consists on:
- $$5$$ MATCHEs $$\times 3$$
- $$1$$ MISSMATCH $$\times -1$$
- $$1$$ GAP $$\times -2$$

Summing all these values we obtain a total score for the whole alignment:

$$(5 \times 3) + (1 \times -1) + (1 \times -2) = 12$$

The same way, for alignment 2:

{% highlight markdown %}
A = PO-INTER
    |---||||
B = P-U-NTER
{% endhighlight %}

we have a score of:

$$(5 \times 3) + (0 \times -1) + (3 \times -2) = 9$$

And, for alignment 3:

{% highlight markdown %}
A = POINTER------
    -------------
B = -------PUNTER
{% endhighlight %}

a score of:

$$(0 \times 3) + (0 \times -1) + (13 \times -2) = -26$$

We conclude that, with the proposed scoring scheme, the alignment from example 1 is clearly the "_best_" of the three, and, according to our definition, it is the "__optimal__" alignment.

## Alternative Scoring Schemes

Note that the score depends on values __arbitrarily__ assigned to each operation. If we had made a different choice of values for our scoring scheme, say:

<table class="data-table">
    <tr><td><b>Operation</b></td><td><b>Alternative Value</b></td></tr>
    <tr><td>MATCH</td><td>0</td></tr>
    <tr><td>MISSMATCH</td><td>0</td></tr>
    <tr><td>GAP</td><td>1</td></tr>
</table>

We would obtain totally different alignment scores:

- Alignment 1: $$(5 \times 0) + (1 \times 0) + (1 \times 1) = 1$$
- Alignment 2: $$(5 \times 0) + (0 \times 0) + (3 \times 1) = 3$$
- Alignment 3: $$(0 \times 0) + (0 \times 0) + (13 \times 1) = 13$$

With this alternative scoring scheme, the alignment 3 becomes the optimal one, even if it does not look at all like a good alignment.

The point here is that the concept of "_optimal alignment_" may depend on what we want to achieve with the alignment, and the scoring scheme must be chosen according these objectives.

In certain scenarios we would prefer the alignment with the smallest number of changes and even smaller number of gaps, which, in our example, would be the first one. In other cases, however, we would want to explicitly avoid certain letter substitutions (e.g. we could have different scores for different pairs of matching letters).

Each application of sequence alignments has a particular goal and, for each goal, a specific scoring scheme must be fine tuned.

## A Biology Inspired Scoring Scheme

If this whole discussion regarding scoring schemes seems too arbitrary or abstract, I want to give you a real world example coming from my previous post: The comparison of biological sequences.

Here, the goal is to "reconstruct" the evolutionary history of each position of the ancestor sequence of nucleotides (for DNA) or amino acids (for proteins) <a name="2eton"/>[[2]](#note2).

To achieve this goal, Computational Biologists dedicate a lot of effort to infer the most probable mutations of nucleotides or amino acids during the evolutionary process. The scoring scheme that reflects the inferred mutation probabilities is represented in [__Substitution Matrices__](https://en.wikipedia.org/wiki/Substitution_matrix).

The following substitution matrix is a common scoring scheme used to align DNA sequences <a name="3eton"/>[[3]](#note3):

<table class="data-table">
    <tr>
        <td>&nbsp;</td>
        <td><span class="A"><b>A</b></span></td>
        <td><span class="C"><b>C</b></span></td>
        <td><span class="G"><b>G</b></span></td>
        <td><span class="T"><b>T</b></span></td>
    </tr>
    <tr><td><span class="A"><b>A</b></span></td><td> 4</td><td>-2</td><td> 0</td><td>-2</td></tr>
    <tr><td><span class="C"><b>C</b></span></td><td>-2</td><td> 4</td><td>-1</td><td> 0</td></tr>
    <tr><td><span class="G"><b>G</b></span></td><td> 0</td><td>-2</td><td> 4</td><td>-2</td></tr>
    <tr><td><span class="T"><b>T</b></span></td><td>-2</td><td> 0</td><td>-2</td><td> 4</td></tr>
    <tr><td><b>INDEL</b></td><td colspan="4">-1</td></tr>
</table>

Each entry in a substitution matrix is related to the probability of observing, in a given time span (usually in the order of thousands of generations), a mutation of the nucleotide in the row into the nucleotide in the column <a name="4eton"/>[[4]](#note4). Note that the diagonal of the matrix corresponds to the event of no mutation which is much more probable than a mutation.

From an alignment point of view, the diagonal corresponds to a __MATCH__ operation, while all other positions correspond to __MISSMATCHES__. Insertion and deletions of nucleotides, which also occur in nature and are commonly called "_indels_", have a specific score and correspond to the __GAP__ operation. See [Hamady 2006](#Hamady2006) for an interesting discussion on the subject.

To see how all of this works, let's take a look into the following alignment between two small sequences:

<div align="center" id="DNA">
<span class="seqname">Seq1&nbsp;</span><span class="G">GG</span><span class="C">C</span><span class="-">-</span><span class="G">GG</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||x-||</span><br/>
<span class="seqname">Seq2&nbsp;</span><span class="G">GG</span><span class="A">A</span><span class="T">T</span><span class="G">GG</span>
</div>

And now let's consider this alternative alignment between the same sequences: 

<div align="center" id="DNA">
<span class="seqname">Seq1&nbsp;</span><span class="G">GG</span><span class="-">-</span><span class="C">C</span><span class="G">GG</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||-x||</span><br/>
<span class="seqname">Seq2&nbsp;</span><span class="G">GG</span><span class="A">A</span><span class="T">T</span><span class="G">GG</span>
</div>

Given the evolutionary model represented by our substitution matrix above, __which one of these two alignments is the most probable one, from a biological point of view?__

In other words, in the ancestry of those two sequences, was there a mutation between "C" $$\leftrightarrow$$ "A" and an a _indel_ of the "T", or alternatively, an _indel_ of the "A" with a mutation between "C" $$\leftrightarrow$$ "T"?

Computing the scores we get respectively:

$$(4 \times 4) + (1 \times -2) + (1 \times -1) = 13$$

and

$$(4 \times 4) + (1 \times 0) + (1 \times -1) = 15$$

Clearly the second alignment is the most probable. It's easy to see why. The "C" $$\leftrightarrow$$ "T" mutation has a score of $$0$$, while the "C" $$\leftrightarrow$$ "A" mutation has a lower score of $$-2$$, thus the former one is the most probable one.

# A Pause for Recap

A few observations are worth recap at this point:

First, the idea of "__optimal alignment__" is closely related with the mathematical concept of [optimization](https://en.wikipedia.org/wiki/Mathematical_optimization), hence its name.

Second, it's easy to see that playing with the values assigned to each pairing will impact the final alignment score dramatically. The choices of values must have in mind the practical objectives of the alignment. Our goal in these posts will be to compare text, so the values of the first example ($$MATCH=3$$, $$MISSMATCH=-1$$, $$GAP=-2$$) that promote matches and penalize a little less mismatches than gaps, tend to work well.

Third, contrary to the biological example presented above, when aligning text, the absolute value of each alignment's score is not important by itself. We are only interested in the alignment with the maximum score.

# Next Steps

Now that we have a formal foundation to define and evaluate alignments, it's time to start exploring the algorithms to perform theses alignments automatically. This will be the subject of the next post.


# Notes

<a name="note1" />__[[1]](#1eton)__ The following exact formula for the total number of alignments between a pair of sequences of sizes $$m$$ and $$n$$ was derived by [Torres et al.](#Torres2004):

$$f(m, n) = \sum^{min(m, n)}_{k=0} {2^k {m \choose k} {n \choose k}}$$

The following Python code computes the number of alignments for two sequences of sizes $$m=7$$ and $$n=6$$:

{% highlight python %}
from math import comb
def f(m, n):
    return sum([(2**k) * comb(m, k) * comb(n, k) for k in range(0, min(m, n))])
print(f(7, 6))
>>> 19377
{% endhighlight %}

<a name="note2" />__[[2]](#2eton)__ The _ancestor_ sequence is the one that originated both sequences being aligned through mutation.

<a name="note3" />__[[3]](#3eton)__ For the pedantic, this matrix is a common scoring scheme for __coding__ DNA which would not be applicable to __non coding__ DNA. This discussion however is way beyond the scope of this post.

<a name="note4" />__[[4]](#4eton)__ In reality the substitution matrix represents the log odds score of the probability of the mutation over the probability of observing the original nucleotide, but this is a technical discussion way beyond the scope of our discussion. For more information see [here](https://en.wikipedia.org/wiki/Substitution_matrix#Log-odds_matrices).

# References

<a name="Torres2004">[\[Torres2004\]](https://pubmed.ncbi.nlm.nih.gov/15018352/)</a> Torres A., Cabada A., Nieto J. (2004). __An Exact Formula for the Number of Alignments Between Two DNA Sequences__. DNA sequence, 14, 427-430.

<a name="Hamady2006">[\[Hamady2006\]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-476)</a> Hamady M., Betterton M., Knight R. (2006). __Using the nucleotide substitution rate matrix to detect horizontal gene transfer.__ BMC bioinformatics. 7:476.