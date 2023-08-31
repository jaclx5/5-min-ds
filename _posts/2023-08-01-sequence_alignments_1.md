---
title: Text Comparing with Sequence Alignment Algorithms - Part 1
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

==== SECOND POST ===

Add this ti the first post.

Update the end of the first post with more info

> **TLDR**: If you are not interested in the biological or in Bioinformatics in general and want to go directly into the algorithms, you can jump directly to the [second post](xxx).


In the first post of this series I discussed the original biological problems that sequence alignments came to solve. Now we will go deeper into the technical description of alignments.

> Check also the [companion notebook](add companion address) (more info at the [end of the post](#code)). 

I this second post I will focus on "Optimal Sequence Alignments". But first lets start by formally defining Sequence Alignment.


# What is an Alignment?

An alignment between two sequences of letters $$A$$ and $$B$$ is a succession of "changes" on $$A$$ that will transform it on $$B$$. To simplify the analysis of alignment algorithms we consider only atomic "changes", i.e., the possible changes can only be one of the following ($$A_i$$ corresponds to the letter in the $$i^{th}$$ position of $$A$$):
- Replace any element $$A_i$$ by any other letter of the vocabulary.
- Delete any element  $$A_j$$.
- Insert an element after or before position $$k$$.

As a simple example lets take two very short sequences "POINTER" and "PUNTER".

[[change the background to make it visually nicer]]
[[check here for inspiration: https://medium.com/@clem.boin/creating-a-minimal-kernel-development-setup-using-qemu-and-archlinux-987896954d84]]

<pre>
<code class="python">
    A = POINTER
    B = PUNTER
</code>
</pre>

A possible alignment could be the obtained by performing the following two changes:

- Replace $$A_2$$ by "U"
- Delete $$A_3$$

The resulting alignment would be:

<pre>
<code class="python">
    A = POINTER
         x-
    B = PU-NTER
</code>
</pre>

Note that although the alignment above seems somehow "natural", there's nothing special about it. The following alternative list of changes:

- Delete $$A_2$$
- Insert an "U" before position $$2$$
- Delete $$A_3$$

Produce an equally valid alignment:

<pre>
<code class="python">
    A         = PO-INTER
                 ---
    B         = P-U-NTER
</code>
</pre>

In the extreme (and clearly silly) case one could consider changing all letters:

- Delete $$A_1$$ to $$A_7$$
- Insert an "P" in position $$1$$
- Insert an "U" after position $$1$$
- Insert an "N" after position $$2$$
- Insert an "T" after position $$3$$
- Insert an "E" after position $$4$$
- Insert an "R" after position $$5$$

Which would produce an also valid (although not very useful as we will see in a moment) alignment:

<pre>
<code class="python">
    A         = POINTER------
                -------------
    B         = -------PUNTER
</code>
</pre>

_In the examples above we added the symbol «-» in the position of the deleted/inserted letters to facilitate the comparison. Those «-», of course, are not part of the sequence._

From the examples bellow we can easily intuit that there are a huge number of possible alignments between two sequences, in fact for those tiny sequences we can expect almost 20000 possible alignments ([check here for an exact formula](end of the post)). The question that arises is how to choose from all those possible alignments the _best_ one?

There is not a single answer to the question. The _best_ alignment depends on what do we want to achieve with it. Intuitively we would prefer the alignment resulting from the smallest number of changes, which in our example would be the first proposed alignment. However, real world applications have specific goals and/or restrictions. Sometimes we aim to minimize insertions or deletions. In other cases we want to avoid specific letters substitutions, in the protein alignment field a lot of work has been done to determine which substitutions are _better_ than others ([see here](https://en.wikipedia.org/wiki/Substitution_matrix)).

The discussion on how to choose the best alignment lead us to the definition of __optimal alignment__


# Optimal Alignment

In order to choose the optimal alignment from the set of all possible alignments we need find a way to assign some kind of score to each alignment so we can compare them all against each other. A simple way to define a score is to evaluate each position of the alignment and assign some value depending on the aligned letters in that position.

lets assume that $$A_i$$ and $$B_i$$ represent the $i^{th}$ position of the __aligned__ sequences $$A$$ and $$B$$. It's easy to see that we have three possible pairing situations:

- A __MATCH__  if positions $$A_i$$ and $$B_i$$ are identical (represented by the _pipe_ character "|").
- A __MISSMATCH__ if positions $$A_i$$ and $$B_i$$ are different (represented by an "x").
- A __GAP__ if either $$A_i$$ or $$B_i$$ contain a "-" after the alignment (conveniently represented by an "-").

As a simple example lets take the first alignment of the example above:

<pre>
<code class="python">
    A = POINTER
        |x-||||
    B = PU-NTER
</code>
</pre>

This alignment consisted on:
- 5 MATCHEs in positions 1, 4, 5, 6 and 7.
- 1 MISSMATCH in position 2
- 1 GAP in position 3

If now we assign an individual and arbitrary value to each pairing say: MATCH=3, MISSMATCH=-1, GAP=-2 and sum all the pairings we obtain a total score for the whole alignment.

The alignment would then score $$(5 \times 3) + (1 \times -1) + (1 \times -2) = 12$$

The same way the for the remaining two examples we would have scores of $$(5 \times 3) + (3 \times -2) = 9$$:

<pre>
<code class="python">
    A         = PO-INTER
                |---||||
    B         = P-U-NTER
</code>
</pre>

And $$(13 \times -2) = -26$$:

<pre>
<code class="python">
    A         = POINTER------
                -------------
    B         = -------PUNTER
</code>
</pre>

In this example the first alignment is clearly the best one of the three (as we would probably expect).

A few observations are worth mentioning at this point.

First, it's easy to see that playing with the values assigned to each pair we can change the final scores dramatically. So the choice of values must have in mind the practical objectives of the alignment. When comparing biological sequences (DNA or Proteins) the goal is to "reconstruct" the evolutionary history of each position of the sequences (nucleotide or amino acid), thus a lot of research was dedicated to infer those values (see the [substitution matrices](link to Wikipedia) for more information). In our case, when comparing text values that promote matches and penalize more or less equally mismatches and gaps tend to work well.

Second, the absolute value of each alignment's score is not important by itself. We are only interested in the maximum score of all alignments. That does not mean that the values as meaning less. Coming back to biological sequences, the values assigned to each match/mismatch are [[check in the refernce here - the log prob of a mutation....]]

Finally, and more important for our discussion here, we still don't know how to __find__ the best score. [[end this sentence and motivate the following]]

# Alignment Algorithms

## Brute Force Approach

Now that we know how to compute a numerical score for any assignment e need to find a way to search for the best overall alignment, a.k.a. the optimal alignment.

Lets start by illustrating a brute force approach to optimal alignment search:

[[graphical element]]

As we can see in each step of the process we have three potential actions:

- Align the next letter from both sequences (either with a MATCH or MISMATCH).
- Align the next letter from $$A$$ sequence adding a GAP to sequence $$B$$.
- Align the next letter from $$B$$ sequence adding a GAP to sequence $$A$$.

After each of the steps the process repeats until one of the sequences is exhausted and GAPs are added to the remaining letters of the other sequence.

As it progresses, the algorithm recursively builds a ternary tree containing in its leaves the complete alignments. When no more expansion are possible we just have to look for the maximum score alignment, we are done!

Of course, due to the huge number of possible alignments, the brute force approach is completely impractical for all but the smallest sequences. However it provides us with a intuition for the dynamic programming algorithm.

## Dynamic Programming

Dynamic Programming was proposed in the 1950s by Richard Bellman (the choice of the name "Dynamic Programming" has a curious motivation as told by [Bellman himself in his autobiography](#Bellman1984)), it is a general optimization method applicable to a large class of problems well beyond sequence alignment. It's an elegant and useful method which description is out of the scope of this post but which is worth to [explore further](https://en.wikipedia.org/wiki/Dynamic_programming).

To get to the intuition behind dynamic programming let change slightly or brute force algorithm. Instead of expanding all possible combinations of actions we will expand only the most "promising" ones:

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

# End Notes

## Exact formula for the number of pair wise sequence alignments

In [this article](#Torres2004) the authors propose the following formula for the total number of possible alignments between a pair of sequences of sizes $$m$$ and $$n$$:

$$f(m, n) = \sum^{min(m, n)}_{k=0} {2^k {m \choose k} {n \choose k}}$$

<pre>
<code class="python">
    from math import comb
    
    def f(m, n):
        return sum([(2**k) * comb(m, k) * comb(n, k) for k in range(0, min(m, n))])

    print(f(7, 6))
</code>
</pre>



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

==== FOURTH POST ===


- **_K-word_ methods (FASTA and BLAST)**
==== FIFTH POST ===

# Beyond Biological Sequences

Just before finishing, one word regarding the application of sequence alignment algorithms to non Biological sequences. It is pretty easy to see that if we replace DNA nucleotides by the letters of any alphabet (or by any tokens) these algorithms can be applied as well.

Of course a few choices need to be made if we want to apply SA algorithms to general sequences (as we will see in the future) but the algorithms can be applied with no changes.
check - for future posts:
text_algorithms-crochemore-1994.pdf
diff algorithms