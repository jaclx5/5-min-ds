---
title: How to Compare Text with Sequence Alignment Algorithms - Part 4
layout: post
author: jaclx5

---

_In the [previous post](/sequence_alignments_3) we start exploring sequence alignment algorithms with the greedy and brute force approaches... now it's time for dynamic programming._

> Check also the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/tree/master/notebooks/sequence_alignment) (more info at the [end of the post](#code)). 

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

## Dynamic Programming

Dynamic Programming was proposed in the 1950s by Richard Bellman (the choice of the name "Dynamic Programming" has a curious motivation as told by [Bellman himself in his autobiography](#Bellman1984)). It is a general optimization method applicable to a large class of problems that go well beyond sequence alignment. It is really worth further [explore](https://en.wikipedia.org/wiki/Dynamic_programming).

The main intuition of the dynamic programming algorithm is that certain problems can be broken into a sequence of sub-problems, in a way that the optimal solution of the original problem can be found by recursively combining the optimal solution of each sub-problem, without the need to recompute them multiple times.

How does dynamic programming applies to sequence alignments?

Remember that there are many different ways to align the $$i$$ first elements of sequence $$A$$ with the $$j$$ first elements of sequence $$B$$ (the sub alignments belonging to the set $$Aln_{i, j}$$). The brute force algorithm explores all of them to find the optimal overall sub alignment.

The dynamic programming approach tells us that once we find the optimal sub-alignment from the set $$Aln_{i,j}$$ (denoted by $$max(Aln_{i,j})$$) we can disregard all other sub-alignments from that set ($$Aln_{i,j} \setminus \{max(Aln_{i,j})\}$$), saving a lot of useless computation. This is because we assume that the overall optimal alignment will only include optimal sub-alignments (which is kind of intuitive).

To make things more concrete lets look at the 6 first depth levels of the brute force alignment tree for sequences <code class="python">ABC</code> and <code class="python">ABXABC</code>:

<div align="center">
    <a href="/images/sequence_alignments/brute_force_compact.png"><img src="/images/sequence_alignments/brute_force_compact.png" height="1400"/></a>
</div>

The coloured nodes (__<span style="color:#0072B2;">blue</span>__ and __<span style="color:#fd6a02;">yellow</span>__) represent all the 13 partial alignments that consumed exactly 2 elements from each sequence (all the sub alignments from the set $$Aln_{2,2}$$), with the <span style="color:#0072B2;">blue</span> box representing the best one of those alignments, i.e., the __optimal sub-alignment of exactly 2 elements of each sequence__, denoted $$max(Aln_{2,2})$$.

According to the Dynamic Programming approach, the final optimal alignment includes __only__ the best sub-alignments, hence there's no need to continue expanding the <span style="color:#fd6a02;">yellow</span> boxes: The optimal solution will not be there <a name="1eton"/>[[1]](#note1).

With this move alone we will save the expansion of $$12 \times 29 =  348$$ nodes, a full $$47\%$$ of all $$743$$ nodes explored by the brute force algorithm. Not Bad!


## The Dynamic Programming Algorithm

Now that we know the basic concept behind Dynamic Programming, let's design a simple sequence alignment algorithm exploting these concepts and see how it actually performs comparing with our Brute Force approach.

The outline of out algorithm will be:
{% highlight markdown %}

1. In each step look for the best non expanded alignment so far [Aln_best_i/j].
2. If there is any other (already expanded) alignment that consumed exactly
   the same i elements of seq. A and j elements of B as [Aln_best_i/j],
   with a better score:
    - Discard [Aln_best_i/j].
3. Otherwise:
    - Expand [Aln_best_i/j] adding its children to the tree.
4. Go back to step 1 until all nodes have been expanded or discarded.
{% endhighlight %}

The sequence of images bellow illustrates some steps pf this algorithm applied to the previous sequences: <code class="python">ABC</code> and <code class="python">ABXABC</code>, respectively our $$A$$ and $$B$$ sequences.

As in the previous post, each image is a step of the algorithm and the boxes (or nodes) represent partial alignments observed at that step. The __<span style="color:#ff0000;">red</span>__ box represents the best "non expanded" alignment observed so far, the __<span style="color:#00ff00;">green</span>__ box represents the expanded node that originated the current state, and the __<span style="color:#0000ff;">blue</span>__ one will represent the solution whenever we find it. But now we have two new colors: __<span style="color:#000000;">black</span>__ represents all the nodes that we will not care to expand and  __<span style="color:#ff00ff;">pink</span>__ represents the best alignment that consumed the respective elements.

********* ADD THE MOVIE HERE *************************


In this example we show only the steps 8, 9 and 10. Until step 8 the Dynamic Programming algorithm works exaclty as the Greedy algorithm, but in the next step we will observe the big difference: The alignment (ABC-/AB-X) is the best non expanded alignment so far with score == 2. However, a previous expanded alignment (ABC/ABX) have consumed the same elements but with a higher score == 5. Hence, there's no point in pursuing the (ABC-/AB-X). We discard it (it's marked in black in step 9) and proceed with the exploration saving the work of exploring that whole branch.

As you can see in the image bellow, the full dynamic programming tree for this alignment takes $$73$$ steps a $$10\times$$ improvement over the $$743$$ steps of the Brute Force algorithm to find the same (hopefully) optimal solution:

<div align="center">
    <a href="/images/sequence_alignments/brute_force_compact.png"><img src="/images/sequence_alignments/full_dynamic_programming.png"/></a>
</div>






# Notes

<a name="note1" />__[[1]](#1eton)__ The solution is also not guaranteed to be a descendant of the blue node as there are other sibling nodes to explore.



# References

<a name="Bellman1984">[Bellman1984]</a>Bellman RE. (1984) Eye of the Hurricane: An Autobiography. World Scientific, Singapore.






<a name="code"/>
# The Code

You can reproduce all the images from this post using the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/demo.ipynb).

In the repository you will find the package [dalt](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/dalt) which allows you to try your own algorithms and see the resulting trees.

__It's important to note that the algorithms in this package ARE NOT efficient and ARE NOT intended to be used in any practical way. They were developed uniquely to illustrate the concepts in this post and to generate graphical representations of alignments for pedagogic purposes.__

If you have questions or ideas feel free to share them in the comments box. 