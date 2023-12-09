---
title: How to Compare Text with Sequence Alignment Algorithms - Dynamic Programming - Part 4
layout: post
author: jaclx5

---

_In the [previous post](/sequence_alignments_3) we start exploring sequence alignment algorithms with the greedy and brute force approaches... now it's time for dynamic programming._

> Check also the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/tree/master/notebooks/sequence_alignment) (more info at the [end of the post](#code)). 

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

## Dynamic Programming

Dynamic Programming was proposed in the 1950s by Richard Bellman (the choice of the name "Dynamic Programming" has a curious motivation as told by [Bellman himself in his autobiography](#Bellman1984). It is a general optimization method applicable to a class of problems that go beyond sequence alignment. For more information on the usage of Dynamic Programming in other fields it is worth to start exploring [here](https://en.wikipedia.org/wiki/Dynamic_programming).

The main intuition behind Dynamic Programming is that certain problems can be broken into a sequence of sub-problems, in a way that the optimal solution of the original problem can be found by recursively combining the optimal solution of each sub-problem, without the need to recompute each sub-problem more than once.

How does dynamic programming applies to sequence alignments?

Remember that there are many different ways to align the $$i$$ first letters of sequence $$A$$ with the $$j$$ first letters of sequence $$B$$, those are all alignments belonging to the set $$Aln_{i, j}$$.

While the Brute Force algorithm expands all of the alignments on $$Aln_{i, j}$$ to find the optimal overall sub alignment, the Dynamic Programming approach tells us that once we find the optimal sub-alignment from the set $$Aln_{i,j}$$, denoted by $$max_{score}(Aln_{i,j})$$, we can disregard all other sub-alignments from that set, i.e., we can ignore the set $$Aln_{i,j} \setminus \{max(Aln_{i,j})\}$$, saving a lot of useless computation. This is because it can be prooved that the overall optimal alignment will only include optimal sub-alignments (which is kind of intuitive). ((*** IT WILL NOT INCLUDE ALL OPTIMAL SUB-ALIGNMENST***, this is a sublety that may become clearer in a moment))

To make things more concrete lets look at the 6 first depth levels of the Brute Force alignment tree for sequences <code class="python">ABC</code> and <code class="python">ABXABC</code>:

<div align="center">
    <a href="/images/sequence_alignments/brute_force_compact.png"><img src="/images/sequence_alignments/brute_force_compact.png" height="1400"/></a>
</div>

The coloured nodes (__<span style="color:#0072B2;">blue</span>__ and __<span style="color:#fd6a02;">yellow</span>__) represent all the 13 partial alignments that consumed exactly 2 letters from each sequence (all the sub alignments from the set $$Aln_{2,2}$$), with the <span style="color:#0072B2;">blue</span> box representing the best one of those alignments, i.e., the __optimal sub-alignment of exactly 2 letters of each sequence__, denoted $$max_{score}(Aln_{2,2})$$.

According to the Dynamic Programming approach, the final optimal alignment includes __only__ best sub-alignments, hence there's no need to continue expanding the <span style="color:#fd6a02;">yellow</span> boxes: The optimal solution will not be there <a name="1eton"/>[[1]](#note1).

With this move alone we will save the expansion of $$12 \times 29 =  348$$ nodes, a full $$47\%$$ of all $$743$$ nodes expanded by the Brute Force algorithm. Not Bad!

## Sequence Alignment as a Graph Traversal


** We have $$i \times j$$ $$(Aln_{i,j})$$ sets

**Not all combinations of i and j values, but only the maximum score of each**

**IT DOES NOT INCLUDE ALL BEST SUB ALIGNMENTS, WE CAN SEE THAT THE PATH FROM $$max(Aln_{0,0})$$ to **

$$max_{score}(Aln_{0,0}) \rightarrow max_{score}(Aln_{1,1}) \rightarrow max_{score}(Aln_{2,2}) \rightarrow max_{score}(Aln_{2,3}) $$
$$max_{score}(Aln_{2,4}) \rightarrow max_{score}(Aln_{2,5}) \rightarrow max_{score}(Aln_{3,5})$$




## The Dynamic Programming Algorithm

Now that we know the basic concept behind Dynamic Programming, let's design a simple sequence alignment algorithm exploiting these concepts and see how it actually performs comparing with our Brute Force approach.

The outline of out algorithm will be:
{% highlight markdown %}

1. In each step look for the best non expanded alignment so far [Aln_best_i/j].
2. If there is any other (already expanded) alignment that consumed exactly
   the same i letters of seq. A and j letters of B as [Aln_best_i/j],
   with a better score:
    - Discard [Aln_best_i/j].
3. Otherwise:
    - Expand [Aln_best_i/j] adding its children to the tree.
4. Go back to step 1 until all nodes have been expanded or discarded.
{% endhighlight %}

The sequence of images bellow illustrates some steps pf this algorithm applied to the previous sequences: <code class="python">ABC</code> and <code class="python">ABXABC</code>, respectively our $$A$$ and $$B$$ sequences.

As in the previous post, each image is a step of the algorithm and the boxes (or nodes) represent partial alignments observed at that step. The __<span style="color:#ff0000;">red</span>__ box represents the best "non expanded" alignment observed so far, the __<span style="color:#00ff00;">green</span>__ box represents the expanded node that originated the current state, and the __<span style="color:#0000ff;">blue</span>__ one will represent the solution whenever we find it. But now we have two new colors: __<span style="color:#000000;">black</span>__ represents all the nodes that we will not care to expand and  __<span style="color:#ff00ff;">pink</span>__ represents the best alignment that consumed the respective letters.

********* ADD THE MOVIE HERE (change the code to show only from min_steps to max_steps *************************


In this example we show only the steps 8, 9 and 10. Until step 8 the Dynamic Programming algorithm works exaclty as the Greedy algorithm, but in the next step we will observe the big difference: The alignment (ABC-/AB-X) is the best non expanded alignment so far with score == 2. However, a previous expanded alignment (ABC/ABX) have consumed the same letters but with a higher score == 5. Hence, there's no point in expanding the (ABC-/AB-X) alignment. We discard it (it's marked in black in step 9) and proceed with the algorithm saving the computation required by the whole branch expansion.

As you can see in the image bellow, the full dynamic programming tree for this alignment takes $$73$$ steps a $$10\times$$ improvement over the $$743$$ steps of the Brute Force algorithm to find the same (hopefully) optimal solution:

<div align="center">
    <a href="/images/sequence_alignments/brute_force_compact.png"><img src="/images/sequence_alignments/full_dynamic_programming.png"/></a>
</div>






# Notes

<a name="note1" />__[[1]](#1eton)__ The solution is also not guaranteed to be a descendant of the blue node as there are other sibling nodes to expand.



# References

<a name="Bellman1984">[\[Bellman1984\]](https://archive.org/details/eyeofhurricaneau0000bell)</a> Bellman RE. (1984) Eye of the Hurricane: An Autobiography. World Scientific, Singapore.






<a name="code"/>
# The Code

You can reproduce all the images from this post using the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/demo.ipynb).

In the repository you will find the package [dalt](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/dalt) which allows you to try your own algorithms and see the resulting trees.

__It's important to note that the algorithms in this package ARE NOT efficient and ARE NOT intended to be used in any practical way. They were developed uniquely to illustrate the concepts in this post and to generate graphical representations of alignments for pedagogic purposes.__

If you have questions or ideas feel free to share them in the comments box. 