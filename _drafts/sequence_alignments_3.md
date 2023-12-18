---
title: How to Compare Text with Sequence Alignment Algorithms - The First Algorithms - Part 3
layout: post
author: jaclx5

greedy1:
    slides:
        - slide:
            image: /images/sequence_alignments/greedy1/step_00.png
            caption: The algorithm starts with an empty alignment, which happens to be the trivially best so far.
        - slide:
            image: /images/sequence_alignments/greedy1/step_01.png
            caption: It expands the "empty" alignment (GREEN box), by applying the three possible operations<br/>
                The (P/P) alignment (RED box) is the best one so far, and will be expanded in the next step.
        - slide:
            image: /images/sequence_alignments/greedy1/step_02.png
            caption: After expanding (P/P), all of the obtained partial alignments have smaller scores.
        - slide:
            image: /images/sequence_alignments/greedy1/step_03.png
            caption: After expanding (PO/PU) we end up with three partial alignments with score == 1.<br/>
                     The algorithm chooses (arbitrarily) the one from the top, in a "depth first"-like search.
        - slide:
            image: /images/sequence_alignments/greedy1/step_04.png
            caption: As it gets no improvement, it will continue expanding the alignments with score == 1...
        - slide:
            image: /images/sequence_alignments/greedy1/step_05.png
            caption: ...and it does it again.
        - slide:
            image: /images/sequence_alignments/greedy1/step_06.png
            caption: Still no luck, the best partial alignment so far has a score == 0.
        - slide:
            image: /images/sequence_alignments/greedy1/step_07.png
            caption: Due to the matching "N"s, we find a partial alignment with score == 3...
        - slide:
            image: /images/sequence_alignments/greedy1/step_08.png
            caption: ...the score keeps increasing with the subsequent matches...
        - slide:
            image: /images/sequence_alignments/greedy1/step_09.png
            caption: ...and increases again...
        - slide:
            image: /images/sequence_alignments/greedy1/step_10.png
            caption: ...until it "consumes" both full sequences, obtaining a final alignment (POINTER/P-UNTER)<br/>
                     which is the algorithm's solution (BLUE box).

greedy2:
    name: greedy
    slides:
        - slide:
            image: /images/sequence_alignments/greedy2/step_00.png
            caption: Again we start with the empty alignment.
        - slide:
            image: /images/sequence_alignments/greedy2/step_01.png
            caption: The partial alignment (A/A) is clearly the best choice at this point.
        - slide:
            image: /images/sequence_alignments/greedy2/step_02.png
            caption: Likewise, matching the "B"s greatly increases the best score so far.
        - slide:
            image: /images/sequence_alignments/greedy2/step_03.png
            caption: Although the C/X mismatch decreases the overall score, the alignment (ABC/ABX) is still the best one so far.
        - slide:
            image: /images/sequence_alignments/greedy2/step_04.png
            caption: In this step the algorithm backtracks due to the accumulation of gaps and mismatches.
        - slide:
            image: /images/sequence_alignments/greedy2/step_05.png
            caption: Here it will try out the second alignment with a score == 4...
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
            caption: Finally, The unassuming alignment (ABC--/ABXAB)<br />
                     ends up as the best alignment to expand.
        - slide:
            image: /images/sequence_alignments/greedy2/step_13.png
            caption: Which leads to a (premature) solution.
---

_In this post, with the help of the formal definitions and the scoring scheme presented in the [previous post](/sequence_alignments_2), I will start exploring some algorithms to automatically find this optimal alignment between two sequences._

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

> Check also the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/tree/master/notebooks/sequence_alignment) (more info at the [end of the post](#code)). 

# Alignment Algorithms

If you have been following the previous posts of this series you have learned quite a bit about alignments. Let's recap a few things we learned so far:
- An alignment between two sequences $$A$$ and $$B$$ is a sequence of operations over $$A$$ that turns it into $$B$$;
- There are a huge number of possible distinct alignments between $$A$$ and $$B$$;
- It is possible to compute a __score__ of how "good" each of those alignments is;
- And, finally, we call the alignment with the highest score the __optimal alignment__.

Now, our next challenge will be to __design an algorithm that finds the optimal alignment between any two given sequences__.

## A Little Bit of Notation 

Before moving on, let's add some notation to help us describe the algorithms to be developed.

An alignment algorithm, as per algorithm definition, is a set of mechanical rules applied, sequentially, step by step. I will illustrate the algorithms in this post with ternary trees like the one bellow:

<div align="center">
    <img src="/images/sequence_alignments/notation_1.png" height="250"/>
</div>

Each tree is a snapshot of a given step of the algorithm. The nodes in the tree correspond to partial alignments being evaluated by the algorithm. Each tree depicts all the nodes (partial alignments) evaluated up to that step, either expand (with children) or not.

Each node in the tree contains the following extra information:

<div align="center">
    <img src="/images/sequence_alignments/notation_2.png" height="250"/>
</div>

We define $$Aln_{i, j}$$ as the set of all alignments that consumed $$i$$ letters from sequence $$A$$ and $$j$$ letters from sequence $$B$$.

Note, as an example, the two following nodes, corresponding to two distinct partial alignments, with distinct sequences of operations:

<div align="center">
    <img src="/images/sequence_alignments/notation_3.png" height="100"/>
</div>

Both of them, although distinct, consumed $$3$$ letters from sequence $$A$$ and $$2$$ letters from sequence $$B$$. Hence, we can say that both alignments belong to the $$Aln_{3, 2}$$ set.

Finally, certain special nodes will be depicted with colored boxes. The __<span style="color:#ff0000;">red</span>__ node represents the best "non expanded" alignment evaluated so far, the __<span style="color:#00ff00;">green</span>__ node represents the expanded node that originated the current state, and the __<span style="color:#0000ff;">blue</span>__ one represents the solution alignment whenever we find it. Note that the a solution requires that all letters of both sequences have been "consumed", thus, all partial alignments are only intermediate steps of an algorithm until a solution is reached.

Now, the algorithms.

## Greedy Approach - A False Start

Our first attempt to build an alignment algorithm will be a straightforward, although _naïve_, way to search for the optimal alignment. This algorithm systematically expands all possible alignment operations in each step of the way, "greedily" choosing the next search direction according to the best score found so far.

The images bellow illustrate this __greedy algorithm__ applied to our [already familiar sequences](/sequence_alignments_2#external1):

{% highlight markdown %}
A = POINTER
B = PUNTER
{% endhighlight %}

Click on the left and right arrows, or in the bottom circles to navigate the steps of the algorithm.

{% include slideshow.html slideshow=page.greedy1 %}

In the first step, the algorithm starts with an _empty_ alignment ($$\text{empty} \in Aln_{0,0}$$) with a score of zero, and "expands" it by applying the three possible operations:

- __INSERT__ in $$A_1$$, obtaining the partial alignment $$\text{(P/-)} \in Aln_{1,0}$$ with a score of -2.
- __MATCH__ between $$A_1$$ and $$B_1$$, obtaining $$\text{(P/P)} \in Aln_{1,1}$$ with a score of 3.
- __DELETION__ of $$A_i$$, obtaining $$\text{(-/P)} \in Aln_{0,1}$$ with a score of -2.

Next, it chooses the best node found so far (thus the __greedy__ in the algorithm's name) and repeats the same expansion process. The algorithm keeps "expanding" the best node, until it finds a solution (i.e. an alignment that consumes all letters of both sequences).

As it progresses, the algorithm recursively builds a ternary tree that contains in its leaves many possible alignments.

Note that, due to the accumulation of gaps and mismatches, in some steps (e.g. steps 4 and 7), the next best alignment to be expanded (red box) isn't always a child of previously expanded one (green box). In these cases the algorithm "backtracks" to expand other branches of the tree.

Although the algorithm seems to be, at a first glance, pretty effective, it has a major flaw. __As it usually occurs with greedy algorithms, it often fails to find the optimal alignment__.

Take, for example, this other pair of sequences, <code class="python">ABC</code> and <code class="python">ABXABC</code>:

{% include slideshow.html slideshow=page.greedy2 %}

As we can see in this example, the greediness of the algorithm, and the fact that it stops as soon as it founds the first solution, leads to the following solution $$(\text{ABC---/ABXABC}) \in Aln_{3,6}$$ __which is not the optimal one__:

{% highlight markdown %}
    ABC---
    ||x---
    ABXABC
{% endhighlight %}

$$(2 \times 3) + (1 \times -1) + (3 \times -2) = -1$$

How do I know this is not the optimal alignment? Because, with a little manual effort, I can find the alignment $$(\text{AB---C/ABXABC}) \in Aln_{3,6}$$ which has a higher score:

{% highlight markdown %}
    AB---C
    ||---|
    ABXABC
{% endhighlight %}

$$(3 \times 3) + (0 \times -1) + (3 \times -2) = 3$$

Of course that, for all but the shortest alignments, it is impossible to manually verify if a specific solution is the optimal one or not. The point I am making here is that, with this simple counter-example, I can prove that the greedy approach some times fails.

At least in the sequence alignment world, greed doesn't pay of.

## Brute Force - Almost There

One way to overcame the limitations of the greedy algorithm is to continue exploring the alignment space, even after having found a first solution, by backtracking and exploring all remaining possibilities. This systematic search over all possibilities is commonly called a "__Brute Force__" strategy, i.e., to test all possibilities and make sure that the optimal alignment is never missed.

Mechanically, both algorithms are pretty similar, the main difference between the greedy and the brute force algorithms is the stopping condition. While the former stops at the first solution, the latter keeps searching until all solutions have been tested.

Of course, due to the huge number of possible alignments, the brute force approach is completely impractical for all but the smallest sequences.

The bravest readers may want to take a look, at their own risk, at the complete brute force expansion of the alignment between <code class="python">ABC</code> and <code class="python">ABXABC</code>:

<div align="center">
    <a href="/images/sequence_alignments/full_brute_force.png"><img src="/images/sequence_alignments/full_brute_force.png" height="500"/></a>
</div>

If you search patiently, you will spot the optimal alignment, mentioned above, in a branch near the center of the image:

<div align="center">
    <img src="/images/sequence_alignments/full_brute_force_detail.png"/>
</div>

# Next Steps

At his point you may ask: "_Why are we wasting time with flawed and impractical algorithms?_"

That's a good question and I have a couple of good answers:
 
First, a simple algorithm, even if it is inefficient and impractical, is a good proof of concept to show that an algorithm to find optimal alignments indeed exists (which is not always the case for every problem).

Second, and most important, understanding how a basic algorithm works allows us to gain some familiarity with the problem, which will make it much easier to understand the really useful algorithms that we will discuss later.

Fortunately, __we still can find the optimal alignment without having to expand all possible alignments__. The solution is inspired in a technique called __Dynamic Programming__ and that will be the main subject of my next post.

<a name="code"/>

# The Code

You can reproduce all the images from this post using the [companion notebook](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/demo.ipynb).

In the repository you will find the package [dalt](https://github.com/jaclx5/jaclx5.github.io/blob/master/notebooks/sequence_alignment/alignment_tree/dalt) which allows you to try your own algorithms and see the resulting trees.

__It is important to note that the algorithms in this package ARE NOT efficient and ARE NOT intended to be used in any practical way. They were developed uniquely to illustrate the concepts in this post and to generate graphical representations of alignments for pedagogic purposes.__

If you have questions or ideas feel free to share them in the comments box. 