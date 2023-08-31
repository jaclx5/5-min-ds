---
title: How to Compare Text with Sequence Alignment Algorithms - Part 1
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

_In this series of posts I discuss a few approaches on how to apply Bioinformatics sequence alignment tools to general non-exact text comparison and search problems._

Every once in a while a "fuzzy" (non-exact) text comparison problem pops up during my busy days. Problems like:

- Searching databases with non-exact text queries;
- Matching text fragments containing mistakes and errors;
- Looking up names in abbreviated forms or riddled with typos;
- Comparing documents against templates.

Recurrently appear on mine or in my colleagues' projects. 

When I ask, as I often do, "_Why don't we use a sequence alignment algorithm here?_" the usual reaction ranges from blank stares to raised eyebrows.

Although sequence alignment algorithms are essential tools in any Bioinformatician's toolbox, it is my experience that most of my fellow data scientists seem to be unaware of them, which always surprises me. Sequence alignment algorithms have been around for more than 50 years, they can be pretty useful to tackle more generic text comparison problems, and several commonly used tools (e.g. the _git diff_ tool and the _Levenshtein_ distance) apply sequence alignment in some a way. Yet, a lot of people don't know them.

It happens so often that I have to explain sequence alignment algorithms to my suspicious colleagues, that I decided to write this series of posts so I can easily point them to. Another selfish reason to write these posts is that I have wanted for a long time to collect in a single place a documented code repository of the algorithms mentioned here so I don't need to re-code them over and over.

I hope you find them helpful too.

In this first post I will present the original Biological motivation that led the to the development of Sequence Alignment algorithms.

In the following posts I will present the technical details of pair-wise "optimal alignment" algorithms. Explain what is an "optimal alignment", why it is such a computationally hard problem and how Dynamic Programming approaches solve it (at least partially), and go in depth into the inner workings of the **Needleman-Wunsch** and **Smith-Watterman** algorithms and its implementation.

Later on I will present non optimal alignment algorithms (e.g. **dot plot** and **k-word** approaches for aligning big sequences) and discuss some creative ways of applying alignment algorithms to generic text and full documents.

# The Biological Context

## The DNA is a (kind of) string

All known cellular life forms use DNA molecules as a medium for storage of genetic information and its transmission to the next generation of cells. DNA molecules consist of a long, linear chain of nucleotides. Each nucleotide contains one of four "bases": adenine, cytosine, guanine and thymine (represented respectively by the letters A, C, T and G). Of course, DNA is much more complex and interesting than this rough description and the curious reader have plenty of high quality [on line](https://www.genome.gov/genetics-glossary/Deoxyribonucleic-Acid]) content to explore further.

The string bellow is a tiny excerpt of the human DNA sequence in the standard nucleotide color code: 

<div align="center" id="DNA">
<span class="seqname">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GG</span><span class="C">CC</span><span class="A">A</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="A">A</span><span class="C">CCC</span><span class="T">T</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><span class="C">CCC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T<br/></span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="T">TT<br/></span><span class="G">GG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="C">CC</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="G">GGG</span><span class="A">AA</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="G">GGG</span><span class="T">T</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span><br/>
</div>

For the curious, this sequence corresponds to the [Human 5S ribosomal RNA (5S rRNA)](https://rfam.org/family/RF00001) gene.

The DNA molecules inside each individual cell of any organism contain the chemical instructions for the growth and functioning of the cell. The number of nucleotides of a DNA molecule can go from a few tens of thousands ($$\sim 10^3$$), in bacteria, to several billions ($$\sim 10^9$$), in complex multicellular organisms. As a reference the human DNA molecule contains around 3 billion nucleotides ($$3 \times 10^9$$).

## Copy with Errors

As any reliable medium for information storage and transmission, DNA molecules must be faithfully copied from one generation (of cell) to the next. To ensure a reliable replication, cells evolved a [diverse and hugely complex mechanism](#Pray2008) for checking and correcting errors.

In spite of all the efforts to ensure low error rates during copy, random mistakes do happen, and the resulting copy of the DNA molecule turns out different from the original one. In this case the child cell contains a different copy from its parent's DNA. Those differences, called __mutations__, can be small, affecting a single nucleotide change, or huge, affecting large chunks of the DNA molecule.

Most of the time mutations occur in regions of the DNA with little or no impact for the descendant cell. Those are called [_neutral mutations_](#Kimura1983).

Less frequently, mutations occur in important regions, with bad or even catastrophic consequences for the next generation cell rendering it less fit of non viable altogether. Those are usually referred to as _deleterious mutations_.

On rare occasions, however, mutations can provide some kind of advantage and make the child cell more "successful" than it's "sister" cells not affected by that mutation. By "successful" here one show consider "more able to produce further copies of itself". We would call them _beneficial mutations_.

In the later case those "happy" mutations become part of the genetic heritage of the new cell and will be transmitted to its offspring.

After many generations, with the accumulation of numerous mutations, the resulting organism will be very different from it's distant cousins that descend from the original sister cells that didn't suffer the original mutation. So different in fact that a new species appear.

The diagram bellow tries to illustrate this process in a simplified way:

<div align="center">
    <img src="/images/sequence_alignments/mutations.png" width="700"/>
</div>

In the diagram each column corresponds to one generation, i.e., one replication step of an imaginary "_species A_" cell. The parent cell (left most column) divides itself into four children that inherit its sequence. During the replication process the top most children suffers a random deleterious mutation ($$T3 \rightarrow A3$$) and dies leaving no descendency. The second children suffers no mutation, the third one suffers a neutral mutation $$(G3 \rightarrow C3)$$ and the fourth one suffers a benefficial mutation ($$A1 \rightarrow T1$$) that provides some advantage in the current environment of the population. The second, third and fourth children go on reproducing and eventually accumulating new random mutations. With time, the N-th great grand children of the bottom cell accumulate so many mutations and become so different from their cousins (remember that the sequence _genotype_ defines some of the visible characteristics of the organism _phenotype_) that they establish a new "_species A_".

Bear in mind that this is a huge and mostly flawed simplification. Real life is a lot more complex and interesting. First, the neutral, deleterious or beneficial nature of a mutation depends on the context in which the cell exists (environment, ecosystem, competition, ...). When the context changes what was once beneficial or neutral can become deleterious and vice versa. Second, the definition of species is often not clear and it's dificult of not impossible to relate specific mutations with the emergence of a new species. And finally, even if we refer to mutations as random, the selective pressure will keep some mutations in the population and eliminate others guiding the evolution according to what in a given moment works best. For a fascinating non-technical discussion on all things related with the subject see Stephen Jay Gould's "[__Wondferful Life__](https://en.wikipedia.org/wiki/Wonderful_Life_(book))", Richard Dawkins' "[__The Selfish Gene__](https://en.wikipedia.org/wiki/The_Selfish_Gene)" and "[__Blind Watchmaker__](https://en.wikipedia.org/wiki/The_Blind_Watchmaker)" and Daniel C. Dennet's "[__Darwin's Dangerous Idea__](https://en.wikipedia.org/wiki/Darwin%27s_Dangerous_Idea)".

The relevant point, for this post, is that __random mutations, more than just errors in the DNA replication process, are a source of biological diversity that allow Nature to explore the possibilities of life__. 


## DNA Mutations are a Tool for Biologists

From the Biologist's point of view, mutations, are also an important source of information. Bellow are some examples of how DNA mutations are used in Biology to better understand life:

### A Window into History

Mutations allow biologists to peek into the past and infer ancestral sequences from the current living organisms.

All species on earth, from bacteria to elephants, are genetically related. They all (including we humans, of course) descend from some ancestral organism conveniently referred to as "Last Universal Common Ancestor", [LUCA](https://en.wikipedia.org/wiki/Last_universal_common_ancestor). Unfortunately, the DNA sequence of the LUCA and of all extinct intermediate ancestors is forever lost, except for some very [special cases](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100745/). All we can access today is the DNA sequence of the currently living species with their similarities and differences.

The sequence alignment shown bellow corresponds to a highly conserved region of DNA of all domains of life, the 5S Ribosomal RNA. We can see that all species share this similar sequence (with some mutations). In particular the two positions indicated by stars at the bottom line of the diagram represent nucleotides conserved in all those species. Probably cells containing a different nucleotide in that position are not very successful.

<div align="left" id="DNA">
<span class="seqname">Human&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><br/>
<span class="seqname">Chimp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="-">-</span><span class="T">TT</span><span class="G">G</span><span class="A">AA</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">TT</span><span class="G">GGG</span><br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><br/>
<span class="seqname">Chicken&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="-">--</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">GG</span><span class="C">C</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><br/>
<span class="seqname">Zebra&nbsp;Fish&nbsp;&nbsp;</span><span class="T">TT</span><span class="A">A</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="-">-</span><span class="A">A</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><br/>
<span class="seqname">Arabidopsis&nbsp;</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="C">CC</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GG</span><span class="T">T</span><span class="A">AA</span><span class="G">GGG</span><span class="T">T</span><span class="G">G</span><span class="C">C</span><span class="T">TT</span><span class="G">GGG</span><br/>
<span class="seqname">Letuce&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="C">CC</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="T">TT</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="G">G</span><span class="C">C</span><span class="T">TT</span><span class="G">GGG</span><br/>
<span class="seqname">Rice&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="C">CC</span><span class="G">G</span><span class="A">AA</span><span class="G">G</span><span class="T">TT</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="G">G</span><span class="C">C</span><span class="T">TT</span><span class="G">GGG</span><br/>
<span class="seqname">Nematode&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">CC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="G">GG</span><span class="C">C</span><span class="A">AA</span><span class="G">G</span><span class="T">TT</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="T">TT</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><br/>
<span class="seqname">Yeast&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">CC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="A">AA</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="G">G</span><span class="T">TTT</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="G">GG</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="A">A</span><br/>
<span class="seqname">E.&nbsp;coli&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="C">CC</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="A">AA</span><span class="G">G</span><span class="T">T</span><span class="G">G</span><span class="A">AAA</span><span class="C">C</span><span class="G">G</span><span class="C">CC</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><br/>
<span class="seqname">Deinococcus&nbsp;</span><span class="C">CC</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="-">-</span><span class="T">T</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="G">G</span><span class="A">AAA</span><span class="C">C</span><span class="A">A</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><br/>
</div>

Computational biologists and phylogeny specialists developed [amazingly smart techniques and models](#Yang2006) to infer the probable DNA of ancestor cells, trace the mutations that led to the existing species DNA sequences and establish phylogenetic relationships between current species.


### Identification of Individuals

Even the less attentive movie or TV series spectator knows that a DNA sample is enough to find a murder. It goes without saying that the DNA sample must be sequenced and somehow compared with existing sequence databases. The fact that each individual carries specific mutations renders possible the identification of an individual.


### Understanding the Relationship Between Genetics and Cell Functions

It is common knowledge that all physical characteristics, function and behavior (a.k.a the phenotype) of an organism is in some way or another determined or affected by the genetic sequence of the DNA. However, the causal relationship between each specific part of the DNA and the cell phenotype remains largely unknown and is still a subject of intense study, by comparing DNA mutations across organisms and phenotypic differences.

One example of this research are the [Genome-wide association studies](#Uffelmann2021) that compare large number of mutations observed in many individual genomes to find statistical correlations with interesting characteristics of those individuals.

# Back to Sequence Alignments

Most of the research work involved in the examples above and in any subject related with DNA sequence analysis requires some kind of __Sequence Alignment__.

To see why is that, take a look at the following complete sequences of the 5S Ribosomal RNA gene from Mouse and Chicken:

<div align="left" id="DNA">
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GG</span><span class="C">CC</span><span class="A">A</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="A">A</span><span class="C">CCC</span><span class="T">T</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><span class="C">CCC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GG</span><span class="C">CC</span><span class="A">A</span><span class="T">T</span><span class="C">CCC</span><span class="A">A</span><span class="C">CCC</span><span class="T">T</span><span class="G">GGG</span><span class="T">T</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="C">CC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">TT</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><span class="G">G</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;</span><span class="*">**********</span><span class="&nbsp;">&nbsp;</span><span class="*">********</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">***************</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><br/>
<br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="T">TT</span><span class="G">GG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="C">CC</span><span class="G">G</span><span class="C">CC</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">GG</span><span class="C">C</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="T">TTT</span><span class="G">GGG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="C">CC</span><span class="A">AA</span><span class="G">G</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">**************</span><span class="&nbsp;">&nbsp;</span><span class="*">****</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><br/>
<br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="T">TT</span><span class="G">GG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="C">CC</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="G">GGG</span><span class="A">AA</span><span class="T">T</span><span class="A">A</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">GG</span><span class="C">C</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="T">TTT</span><span class="G">GGG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><br/>
<br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="C">CC</span><span class="G">GGG</span><span class="T">T</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="G">G</span><span class="A">AA</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GGG</span><span class="T">TT</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="G">G</span><span class="T">TTTT</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;</span><span class="*">*</span><br/>
</div>

Comparing both sequences we see that the first 100 nucleotides seem pretty similar with a 65% match, the remaining 65 nucleotides though are much less similar with 21%. In average the whole gene has a 27% match which is equivalent to comparing two random sequences.

How can it be? If these genes are highly conserved we would expect a much higher similarity. Indeed with a some adjustments, i.e., adding a few gaps here and there in the sequence we obtain the following alignment:

<div align="left" id="DNA">
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GG</span><span class="C">CC</span><span class="A">A</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="A">A</span><span class="C">CCC</span><span class="T">T</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><span class="C">CCC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GG</span><span class="C">CC</span><span class="A">A</span><span class="T">T</span><span class="C">CCC</span><span class="A">A</span><span class="C">CCC</span><span class="T">T</span><span class="G">GGG</span><span class="T">T</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="C">CC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="-">-</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;</span><span class="*">**********</span><span class="&nbsp;">&nbsp;</span><span class="*">********</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">***************</span><span class="&nbsp;">&nbsp;</span><span class="*">******</span><br/>
<br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><span class="-">-</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="T">TT</span><span class="G">GG</span><span class="-">-</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="-">-</span><span class="C">CC</span><span class="G">G</span><span class="C">CC</span><span class="A">AA</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">GG</span><span class="C">C</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="T">TTT</span><span class="G">GGG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="C">CC</span><span class="A">AA</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">******</span><span class="&nbsp;">&nbsp;&nbsp;</span><span class="*">********</span><span class="&nbsp;">&nbsp;</span><span class="*">***********</span><span class="&nbsp;">&nbsp;</span><span class="*">****</span><span class="&nbsp;">&nbsp;</span><span class="*">********</span><span class="&nbsp;">&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;</span><span class="*">****</span><br/>
<br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><span class="-">-</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="T">TT</span><span class="G">GG</span><span class="-">-</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="-">-</span><span class="C">CC</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="G">G</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">GG</span><span class="C">C</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="T">TTT</span><span class="G">GGG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="C">CC</span><span class="T">T</span><span class="G">G</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">******</span><span class="&nbsp;">&nbsp;&nbsp;</span><span class="*">********</span><span class="&nbsp;">&nbsp;</span><span class="*">***********</span><span class="&nbsp;">&nbsp;</span><span class="*">****</span><span class="&nbsp;">&nbsp;</span><span class="*">********</span><span class="&nbsp;">&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;</span><span class="*">****</span><br/>
<br/>
<span class="seqname">Mouse&nbsp;&nbsp;&nbsp;</span><span class="G">GG</span><span class="A">AA</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="G">GGG</span><span class="T">T</span><span class="-">-</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span><br/>
<span class="seqname">Chicken&nbsp;</span><span class="G">GG</span><span class="A">AA</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GGG</span><span class="T">TT</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="G">G</span><span class="T">TTTT</span><br/>
<span class="seqname">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*******</span><span class="&nbsp;">&nbsp;</span><span class="*">**</span><span class="&nbsp;">&nbsp;</span><span class="*">*</span><span class="&nbsp;">&nbsp;</span><span class="*">*****</span><span class="&nbsp;">&nbsp;&nbsp;&nbsp;&nbsp;</span><span class="*">*</span><br/>
</div>

Now, we get an 80% match on the whole gene.

That's what alignments are important for in Bioinformatics, to allow Biologists to probe the real evolutionary relationship between sequences, namely to understand (among other things):

- Which nucleotides didn't change (are conserved) since the last common ancestor [1](#note1).
- Which nucleotides were replaced by other nucleotides.
- Which nucleotides were inserted/deleted [2](#note2).
- How evolutionary distant two sequences (and species) are.

And this is not a small feat!

## Types of Sequence Alignments

Now that we understand the biological motivation and the importance of sequence alignments the natural follow up question would be "_How do we perform an alignment between two sequences?_"

At first glance, aligning two sequences seems a pretty simple operation. Put one sequence over another and shift them around adding a space here and there until you get a good match. This apparent simplicity, however, hides some elusive complexities. For instance: How many alignments should one try until obtaining a good alignment? How can we be sure that the obtained alignment is really a good one?

[[provide an example - simple alignments are pretty obvious but sometimes there are diferences how do we distinguish mutations from deletion followed by an inserition?]]

A way to answer this questions is to define a way to quantify the _how good_ each alignment is and to come up with an algorithm that is sure to find _good_ alignments.

For this we have good and bad news. The good news is that indeed we can define the _goodness_ of any alignment and we can even to devise an algorithm to find the best alignment of all, i.e., the "Optimal Alignment". The bad news is that the problem of finding the optimal alignment between a pair of sequences is a very hard problem [[3]](#note3) which renders optimal alignment algorithms impractical for too big sequences.

The balance between obtaining the optimal alignment for limited size sequences or approximate alignment for big sequences, led to the development of several types of alignment algorithms:

- **Optimal alignment algoritms**
- **_Dot Plot_ algorithms**
- **_K-word_ algorithms**

On the upcoming posts we will explore each one of these types of algorithms in detail.

> __Before finishing let me mention that in this series of posts I am only considering pairwise sequence alignments. Multiple sequence alignment algorithms, like [Clustal](http://www.clustal.org/), also exist and are pretty useful, but we will not explore them__.



# References

- <a name="Pray2008"/>Pray, L. [__DNA Replication and Causes of Mutation__](https://www.nature.com/scitable/topicpage/dna-replication-and-causes-of-mutation-409/), 2008, Nature Education 1 (1):214

- <a name="Kimura1983"/>Kimura, M. [__The neutral theory of molecular evolution__](https://en.wikipedia.org/wiki/Neutral_theory_of_molecular_evolution), 1983, Cambridge University Press.

- <a name="Yang2006"/>Yang Z. [__Computational Molecular Evolution__](https://academic.oup.com/book/34842), 2006, Oxford University Press.

- <a name="Uffelmann2021"/>Uffelmann E, Huang QQ, Munung NS, de Vries J, Okada Y, Martin AR, Martin HC, Lappalainen T and Posthuma D, 2021, [__Genome-wide association studies__](https://www.nature.com/articles/s43586-021-00056-9), Nature Reviews Methods Primers volume 1, 59.

# Notes

__<a name="note1">[1]</a>__ Note that whenever we find the same nucleotide in a similar (homologous) position in two sequences we cannot exclude the possibility of the occurrence of reversible mutations.

__<a name="note2">[2]</a>__ In reality we cannot distinguish, from sequence alone, an insertion in one sequence's ancestor from a deletion in another sequence's ancestor.

__<a name="note3">[3]</a>__ Optimal alignment algorithms have a time and memory complexity of $$\mathcal{O}(n \times m)$$.