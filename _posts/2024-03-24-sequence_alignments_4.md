---
title: How to Compare Text with Sequence Alignment Algorithms - Dynamic Programming - Part 4
layout: post
author: jaclx5

dp1:
    slides:
        - slide:
            image: /images/sequence_alignments/dp/step_08.png
            caption: This is the same tree as the one<br/>
                     from the Greedy Algorithm at step 8.<br/>
                     The node (ABC-/AB-X) would be expanded next...
        - slide:
            image: /images/sequence_alignments/dp/step_09.png
            caption: ...however, the node (ABC/ABX) is the best<br/>
                     of the set Aln(3,3) so we can safely<br/>
                     ignore the node (ABC-/AB-X)...
        - slide:
            image: /images/sequence_alignments/dp/step_10.png
            caption: ...and, for the same reason,<br/>
                     we can also ignore the node (AB-C/ABX-).
---

_In the [previous post](/sequence_alignments_3) we started exploring sequence alignment algorithms with the Greedy and Brute Force approaches... now it's time for Dynamic Programming._

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

> Check also the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/tree/master/notebooks/sequence_alignment) (more info at the [end of the post](#code)). 

# Dynamic Programming

The Dynamic Programming technique was proposed in the 1950's by Richard Bellman as a general optimization method applicable to a broad class of problems <a name="1eton"/>[[1]](#note1). The choice of the name "Dynamic Programming" has a curious marketing-motivated story told by Bellman himself in his [autobiography](#Bellman1984).

The main intuition behind Dynamic Programming is the following: Certain problems can be broken into a sequence of sub-problems in a way that the optimal solution of the original problem can be found by recursively combining each sub-problem's optimal solution. This approach avoids the effort of recomputing each sub-problem more than once.


# Dynamic Programming in Sequence Alignments

The description above is quite abstract and some how confusing, I agree. To make it more understandable let's see how it applies Dynamic Programming to sequence alignments.

We saw in the previous post that there are many different ways to align the $$i$$ first letters of sequence $$A$$ with the $$j$$ first letters of sequence $$B$$. Those are all alignments belonging to the set $$Aln_{i,j}$$.

The Brute Force algorithm expands all alignments belonging to $$Aln_{i, j}$$ (for all combinations if $$i$$ and $$j$$) to find the optimal overall alignment.

The Dynamic Programming approach, however, tells us that when we find the optimal sub-alignment from the set $$Aln_{i,j}$$ - denoted by $$max_{score}(Aln_{i,j})$$ - we can disregard all other sub-alignments from that set in all subsequent computations, i.e., we only need to expand the alignments that start from $${max_{score}(Aln_{i,j})}$$ and ignore all othter alignments from the remaining set $$Aln_{i,j} \setminus {max_{score}(Aln_{i,j})}$$ <a name="2eton"/>[[2]](#note2), avoiding a lot of useless computation.

This works because it can be proved that all partial alignments leading to the optimal alignment are __optimal partial alignments__ themselves (which is kind of intuitive).

Just a small parenthesis in the the reasoning flow: I'm not saying that the optimal alignment is the sequence of __all__ optimal partial alignments. I'm just saying that the partial alignments that lead to the optimal alignment are all optimal. This is a subtlety that, hopefully, will become clearer in a moment.

To get a visual representation of these ideas take a look at the first six depth levels of the Brute Force alignment tree, for the sequences <code class="python">ABC</code> and <code class="python">ABXABC</code>:

<div align="center">
    <a href="/images/sequence_alignments/brute_force_compact.png"><img src="/images/sequence_alignments/brute_force_compact.png" height="1400"/></a>
</div>

The color nodes (__<span style="color:#0072B2;">blue</span>__ and __<span style="color:#fd6a02;">yellow</span>__) represent all the 13 partial alignments that consumed exactly 2 letters from each sequence, i.e., all the partial alignments from the $$Aln_{2,2}$$ set. The <span style="color:#0072B2;">blue</span> node is the best one of those partial alignments, i.e., the __optimal partial alignment of exactly 2 letters of each sequence__, denoted $$max_{score}(Aln_{2,2})$$.

Remember that the Brute Force algorithm will expand all of these 13 partial alignment and their descendants.

With the Dynamic Programming, however, we know that the optimal alignment will __only__ include __optimal partial alignments__. Hence, there is no need to expand any of the <span style="color:#fd6a02;">yellow</span> nodes. We may just ignore them (and their $$12 \times 29$$ children). The optimal alignment will not be a descendant of any of these branches. __We only have to expand the <span style="color:#0072B2;">blue</span> node__.

Note that, as we saw above, the optimal alignment may not descend from the <span style="color:#0072B2;">blue</span> node. It can descend from some other black node in the tree. But that we can't know at this stage, we surely need to expand it.

With this tree pruning move, allowed by the Dynamic Programming, we will avoid the expansion of $$12 \times 29 =  348$$ nodes, a full $$47\%$$ of all $$743$$ nodes expanded by the Brute Force approach. Not Bad!


## The Dynamic Programming Algorithm

Now that we understand the idea behind Dynamic Programming, let's design a simple sequence alignment algorithm exploiting it and see how it actually performs comparing with our Brute Force approach.

Here is the outline of the algorithm:

1. Start with a tree containing only the initial partial alignment, i.e., the one before any alignment: $$start \in Aln_{0,0}$$.

1. Search for the __best non expanded partial alignment__, and record the set $$Aln_{i,j}$$ to which it belongs. Let's call it $$BNEPA_{i,j}$$ (in the first iteration this will be the start alignment with $$i=0$$ and $$j=0$$).

1. Search for the __best expanded partial alignment__ from the set $$Aln_{i,j}$$. Let's call it $$BEPA_{i,j}$$ (in the first iteration this will be an empty set as we didn't expand any alignment yet).

1. If $$Score(BNEPA_{i,j}) >= Score(BEPA_{i,j})$$:
    - Expand the $$BNEPA_{i,j}$$ and add its children to the tree.

1. Otherwise:
    - Ignore the $$BNEPA_{i,j}$$ in all further searches.

1. Go back to step 2 until all nodes have been expanded or discarded.

1. The solution will be the best alignment that consumed all letters from $$A$$ and $$B$$.


The sequence of images bellow illustrates the steps 8, 9 and 10 of the Dynamic Programming algorithm applied to the previous sequences: <code class="python">ABC</code> and <code class="python">ABXABC</code>.

In the previous post we defined that the __<span style="color:#ff0000;">red</span>__ node corresponds to the best "non expanded" partial alignment (BNEPA), the __<span style="color:#00ff00;">green</span>__ one was the node expanded in the previous step, and the __<span style="color:#0000ff;">blue</span>__ one is the optimal alignment, whenever we find it. But now we have two new colors: __<span style="color:#ff00ff;">pink</span>__ that corresponds to the best expanded alignment of the respective set $$Aln_{i,j}$$, and __<span style="color:#000000;">black</span>__ corresponding to all the nodes that we will not expand thanks to the Dynamic Programming insight.

{% include slideshow.html slideshow=page.dp1 %}

Until step 8 (first frame) the Dynamic Programming algorithm works exactly as the Greedy algorithm.

In step 9 (second frame) the Dynamic Programming diverges from the Greedy: The $$BNEPA_{3,3}$$ $$(ABC-/AB-X)$$ with $$score=2$$ has a lower score than the $$BEPA_{3,3}$$ $$(ABC/ABX)$$ with $$score=5$$. Hence, there's no point in expanding the alignment $$(ABC-/AB-X)$$ and we discard it (marking it with a __black__ box) and avoid all the effort associated with the respective expansion.

In step 10 (third frame) the same logic is applied to the new $$BNEPA_{3,3}$$ $$(AB-C/ABX-)$$ that also has a lower score ($$score=2$$) than the $$BEPA_{3,3}$$ $$(ABC/ABX)$$ with $$score=5$$. Again, we discard it, mark it with a __black__ box, and save a lot of work.

We can summarize the Dynamic Programming principle as: "Ignore any node that is not the best alignment of the respective set".

As you can see in the image bellow, the full Dynamic Programming tree for this alignment takes _only_ $$73$$ steps. This is a $$10\times$$ improvement over the $$743$$ steps of the Brute Force algorithm:

<div align="center">
    <img src="/images/sequence_alignments/full_dynamic_programming.png"/>
</div>

Note, as expected, that all the "ancestors" of the optimal alignment (the __<span style="color:#0072B2;">blue</span>__ node) are the best partial alignments of each corresponding set:

$$(/) \in Aln_{0,0}, score=0$$

$$(A/A) \in Aln_{1,1}, score=3$$

$$(AB/AB) \in Aln_{2,2}, score=6$$

$$(AB-/ABX) \in Aln_{2,3}, score=4$$

$$(AB--/ABXA) \in Aln_{2,4}, score=2$$

$$(AB---/ABXAB) \in Aln_{2,5}, score=0$$

$$(AB---C/ABCABC) \in Aln_{3,6}, score=3$$

Also note that, as we mentioned before, not all optimal partial alignments are part of the optimal alignment "path" (e.g. $$(ABC/ABX) \in Aln_{3,3}, score=5$$ is in a different branch of the tree that does not lead to the optimal alignment).

# Next Steps: Needleman-Wunsch

The algorithm presented above, although simple to understand, is neither elegant or efficient. Hopefully, it provides a pedagogical understaning of the Dynamic Programming applied to sequence alignment, which should make it much easier to understand the real deal in terms of Dynamic Programming Sequence Alignment: the Needleman-Wunsch algorithm, which we will discuss in the next post.


# Notes

<a name="note1" />__[[1]](#1eton)__ For more information on the usage of Dynamic Programming in other fields it is worth to take a look [here](https://en.wikipedia.org/wiki/Dynamic_programming).

<a name="note2" />__[[2]](#2eton)__ The operator "$$\setminus$$" in the formula $$X \setminus Y$$ means: "_all elements of $$X$$ except $$Y$$_".

# References

<a name="Bellman1984">[\[Bellman1984\]](https://archive.org/details/eyeofhurricaneau0000bell)</a> Bellman RE. (1984) Eye of the Hurricane: An Autobiography. World Scientific, Singapore.

<a name="code"/>
# The Code

You can reproduce all the images from this post using the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/demo.ipynb).
