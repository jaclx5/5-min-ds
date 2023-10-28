---
title: Using Sequence Alignment Algorithms to Compare Text - Part 1
layout: post
author: jaclx5
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

<style>
    #DNA {
        font-family:monospace;
        font-size: 20px;
    }
    .seqname {
        color: black;
    }
    .A {
        color: orange;
    }
    .C {
        color: red;
    }
    .G {
        color: blue;
    }
    .T { 
        color: green;
    }
</style>

_In this series of post I will discuss how to apply Bioinformatics sequence alignment algorithms to general text comparison and search problems, the ones we typically encounter on document processing and NLP projects._

Every once in a while, in the middle of a busy work day, a "fuzzy" (non-exact) text comparison problem pops up:

- How to search a set of documents or data set rows with a non-exact text query?
- How to match text fragments containing mistakes and errors?
- How to look up names in abbreviated form or riddled with typos form a list?
- How to compare documents against templates?

Questions like these appear time and time again on our data science projects.

When I suggest to my colleagues: "_Why don't we use a sequence alignment algorithm here?_", as I often do, the usual reactions range from blank stares to raised eyebrows.

Although sequence alignment algorithms are essential tools in any Bioinformatician's toolbox, it is my experience that many of my fellow data scientists seem to be unaware of them, which always surprises me. Sequence alignment algorithms have been around for more than 50 years and they can be pretty useful to tackle more generic text comparison problems. Several commonly used tools (e.g. the _git diff_ tool) and techniques  (e.g. the _Levenshtein_ distance) apply ideas similar to the ones used in sequence alignment in some way. Yet, many data practitioners still don't know them.

It happens so often that I have to explain sequence alignment algorithms to my suspicious colleagues that I decided to write this series of posts to point them to.

# Roadmap

In this first post I will present the __original Biological motivation__ that led to the development of Sequence Alignment algorithms.

In the following posts I will __formally define sequence alignments__, explain what is an __optimal alignment__, and present the technical details of __pair-wise optimal alignment algorithms__. I will also discuss why searching for an optimal alignment is such a computationally expensive problem, and how __Dynamic Programming__ give us a huge help on solving (at least partially) it. To illustrate the Dynamic Programming approach I will dig into the inner workings of the **Needleman-Wunsch** and **Smith-Watterman** algorithms and their implementations. Finally, I plan to discuss non-optimal alignment algorithms (e.g. **dot plot** and **k-word** approaches for aligning big sequences) and show some examples of how to apply alignment algorithms to generic text and full documents.

I hope you stay with me and find this subject interesting and, maybe, helpful.

# The Biological Context

## The DNA is a (kind of) string

All known cellular life forms use DNA molecules as a medium for storage and transmission of genetic information to the next generation of cells. DNA molecules consist of a long, linear chain of nucleotides. Each nucleotide contains one of four "bases": Adenine, Cytosine, Guanine and Thymine (represented, respectively, by the letters <span id="DNA" class="A">A</span>, <span id="DNA" class="C">C</span>, <span id="DNA" class="G">G</span> and <span id="DNA" class="T">T</span>).

The DNA molecule inside each individual cell of any organism contains the chemical instructions for the growth and functioning of the cell. The number of nucleotides in a DNA molecule can go from a few tens of thousands ($$\sim 10^3$$), in bacteria, to several billions ($$\sim 10^9$$), in complex multicellular organisms. As a reference the human DNA molecule contains around 3 billion nucleotides ($$3 \times 10^9$$).

The string bellow is a tiny excerpt of the human DNA sequence, corresponding to the [Human 5S ribosomal RNA (5S rRNA)](https://rfam.org/family/RF00001) gene:

<div align="center" id="DNA">
<span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="G">GG</span><span class="C">CC</span><span class="A">A</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="A">A</span><span class="C">CCC</span><span class="T">T</span><span class="G">G</span><span class="A">AA</span><span class="C">C</span><span class="G">G</span><span class="C">C</span><span class="G">G</span><span class="C">CCC</span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">G</span><span class="T">T</span><span class="C">C</span><span class="T">T<br/></span><span class="G">G</span><span class="A">A</span><span class="T">T</span><span class="C">C</span><span class="T">T</span><span class="C">C</span><span class="G">GG</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="A">AA</span><span class="G">G</span><span class="C">C</span><span class="A">A</span><span class="G">GGG</span><span class="T">T</span><span class="C">C</span><span class="G">GGG</span><span class="C">CC</span><span class="T">T</span><span class="G">GG</span><span class="T">TT</span><span class="A">A</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="C">C</span><span class="T">TT<br/></span><span class="G">GG</span><span class="A">A</span><span class="T">T</span><span class="G">GGG</span><span class="A">A</span><span class="G">G</span><span class="A">A</span><span class="C">CC</span><span class="G">G</span><span class="C">CC</span><span class="T">T</span><span class="G">GGG</span><span class="A">AA</span><span class="T">T</span><span class="A">A</span><span class="C">CC</span><span class="G">GGG</span><span class="T">T</span><span class="G">G</span><span class="C">C</span><span class="T">T</span><span class="G">G</span><span class="T">T</span><span class="A">A</span><span class="G">GG</span><span class="C">C</span><span class="T">T</span>&nbsp;<br/>
</div>

DNA is, obviously, much more complex and interesting than this rough description and the curious reader has plenty of high quality [on line](https://www.genome.gov/genetics-glossary/Deoxyribonucleic-Acid]) content to explore further.

## Copy with Errors

Any medium of storage and transmission of information needs to be reliable. DNA is not an exception. DNA molecules must be faithfully copied from a cell to its descendancy, over and over again. To ensure the reliable replication of DNA, cells evolved [several hugely complex mechanisms](#Pray2008) of error checking and correction.

In spite of all the efforts to ensure low error rates during copy, random mistakes do happen, and the resulting copy of the DNA molecule may turn out different from the previous generation one. When errors occur, the child cell receives a different copy from its parent's DNA. Those differences, called __mutations__, range from small, affecting a single nucleotide change, to huge, affecting large chunks of the DNA molecule.

Most of the time mutations occur in regions of the DNA with little or no impact for the descendant cell. Those are called [__neutral mutations__](#Kimura1983).

Less frequently, mutations occur in important regions, with bad or even catastrophic consequences for the cell, rendering it less fit to generate offspring, or even killing the cell. Those are usually referred to as __deleterious mutations__.

On even rarer occasions a mutation can provide some kind of advantage rendering the child cell more "successful" than its "sister" cells not affected by the mutation. In this context "successful" means "more able to produce further copies of itself". We call these __beneficial mutations__.

Any mutation that doesn't kill the cell becomes part of its genetic heritage, and it is transmitted to its offspring.

After many generations, with the accumulation of numerous mutations, the resulting organism will be very different from its distant cousins that descend from the common ancestor. So different in fact that a new species may appear.

The diagram bellow illustrates this process in a simplified way:

<div align="center">
    <img src="/images/sequence_alignments/mutations.png" width="700"/>
</div>

In the diagram the _X-axis_ represents time and each column corresponds to one generation, i.e., one replication step of an imaginary "_species A_". The parent cell (the single sequence in the left most column) divides itself into four children that inherit its sequence. During the replication process the top most children suffers a random deleterious mutation ($$T3 \rightarrow A3$$) and dies leaving no descendants. The second children suffers no mutation, the third one suffers a neutral mutation $$(G3 \rightarrow C3)$$ and the fourth one suffers a beneficial mutation ($$A1 \rightarrow T1$$) that provides some advantage in the current environment of the population. The second, third and fourth children go on reproducing and eventually accumulate more random mutations. With time, the $$N^{th}$$ great grand children of the fourth cell accumulate so many mutations, and become so different from their cousins, that they constitute a new "_species B_".

Bear in mind that this is a huge and flawed simplification. Real life is a lot more complex and interesting. First, the neutral, deleterious or beneficial nature of a mutation depends on the context in which the cell exists (environment, ecosystem, competition, ...). When the context changes what once was beneficial may become deleterious and vice versa. Second, the definition of species is often not clear, and it is difficult, if not impossible, to relate any specific mutation with the emergence of a new species. Finally, even if mutations occur randomly, the selective pressure will keep some of them in the population and eliminate others, _guiding_ the evolution according to what, in a given moment, "works better". For a fascinating non-technical discussion of the subject I strongly recommend Stephen Jay Gould's "[__Wondferful Life__](https://en.wikipedia.org/wiki/Wonderful_Life_(book))", Richard Dawkins' "[__The Selfish Gene__](https://en.wikipedia.org/wiki/The_Selfish_Gene)" and "[__Blind Watchmaker__](https://en.wikipedia.org/wiki/The_Blind_Watchmaker)" and Daniel C. Dennet's "[__Darwin's Dangerous Idea__](https://en.wikipedia.org/wiki/Darwin%27s_Dangerous_Idea)".

## DNA Mutations are a Tool for Biologists

From the Biologist's point of view, mutations are an important source of information, that helps us to understand life.

### A Window into History

Mutations allow biologists to peek into the past to infer the DNA sequences of the ancestors of currently living organisms.

All species on earth, from bacteria to elephants, are genetically related. They (we) all descend from some ancestral organism, conveniently referred to as "Last Universal Common Ancestor" or [LUCA](https://en.wikipedia.org/wiki/Last_universal_common_ancestor). Unfortunately, the DNA sequence of LUCA, and of all extinct intermediate species between LUCA and today's species, is forever lost (with the exception of a very few [special cases](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100745/)). All we can observe today is the DNA sequence of the currently living species with their similarities and differences.

The alignment bellow contains sequences of a fragment of the 5S Ribosomal RNA gene (the one previously shown) for 12 species of all domains of life:

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

The 5S rRNA gene is a highly conserved region of DNA. From the above alignment it is easy to see that all species share a similar sequence. In particular, the two positions, indicated by stars at the bottom line of the diagram, represent nucleotides conserved in all those species. Probably, cells containing a different nucleotide in that position are not very successful.

Relying on these type of alignments, and comparing the similarities and mutations between sequences, computational biologists and phylogeny specialists developed some [amazingly smart techniques and models](#Yang2006) to infer the probable DNA of ancestor cells, to trace the mutations that led to the existing DNA sequences and to establish phylogenetic relationships between current species.

### Identification of Individuals

Even the less attentive spectator of popular movie or TV series knows that a DNA sample is enough to find the murder.

Discounting some fictional exaggeration, it is true that sequencing a DNA sample and comparing it with existing sequence databases allows for the identification of individuals.

The fact that each individual carries specific mutations renders possible this identification with a reasonable degree of confidence.

### Understanding the Relationship Between Genetics and Cell Functions

It is common knowledge that all physical characteristics, function and behavior of an organism (a.k.a phenotype) is in some way determined or affected by its genetic sequence. However, the causal relationship between each specific nucleotide or region of the DNA and the organism's phenotype remains largely unknown and is still a subject of intense research.

Examples of this type of research are the [Genome-wide association studies](#Uffelmann2021) that compare large number of mutations observed in many individual genomes to find statistical correlations with interesting characteristics of those individuals.

Those are just a few examples of the possibilities opened by abstracting DNA as a string of letters (nucleotides), and comparing the differences (mutations) occurring on those strings.

# Back to Sequence Alignments

From everything that was said before, it is important to retain that __random mutations, more than just errors in the DNA replication process, are a source of biological diversity which allows Nature to explore the possibilities of life__. __Mutations are kept in the _written_ record of the DNA sequence and can be compared with those of other individuals and species__ and, finally, __sequence alignments are the key tool that allows this comparison__.

To see why alignments are so important, take a look at the following, __non aligned__, complete sequences of the (yet again) 5S Ribosomal RNA gene from Mouse and Chicken:

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

Comparing both sequences, _as is_ (i.e. without any alignment), we observe that the first 100 nucleotides seem pretty similar with 65% of nucleotides match. The remaining 65 nucleotides, however, are much less similar with only 21% match. Averaging over the whole gene only 27% of all nucleotides seem to be conserved between the two species.

This is too low for such a conserved sequence! 27% is the percentage of similarity one would expect from two random sequences.

How can it be? If these genes are highly conserved we expect a much higher similarity. Indeed with a some adjustments, i.e., adding a few gaps here and there in both sequences, we obtain the following alignment:

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

Now, we can see that the real similarity is around 80% on the whole gene. What we did here was to align both genes, in order to maximize their similarity

This is what alignments are used for in Bioinformatics: To allow Biologists to compare sequences in a meaningful way, and to probe the real evolutionary relationship between sequences, namely to understand (among other things):

- Which nucleotides are conserved (didn't change) since the last common ancestor [[1]](#note1).
- Which nucleotides were replaced by other nucleotides.
- Which nucleotides were inserted/deleted in any of the sequences [[2]](#note2).
- How evolutionary distant two sequences (and species) are.

And this is not a small feat!

## Types of Sequence Alignments

Now that we understand the biological motivation and the importance of sequence alignments the natural follow up question is "_How do we perform an alignment between two sequences?_"

At first glance, aligning two sequences seems a pretty simple operation. Put one sequence over another and shift them around, adding a space here and there, until you get a good match. This apparent simplicity, however, hides an elusive complexity. For instance: How many alignments should one try until obtaining a good alignment? How can we be sure that the obtained alignment is really a good one?

To answer these questions we first need to establish a way to quantify _how good_ each alignment is, and after that to come up with an algorithm that is sure to find the _best_ alignment.

I have good and not so good news.

The good news is that we can easily quantify the _goodness_ (pun no intended) of any alignment. We can even identify the best alignment of all, i.e., the "Optimal Alignment".

The less good news is that the challenge of finding the optimal alignment between a pair of sequences is a computationally hard problem [[3]](#note3), which renders optimal alignment algorithms impractical for really big sequences (e.g. genome size sequences).

The balance between obtaining the optimal alignment for limited size sequences or an approximate alignment for big sequences, led to the development of several types of alignment algorithms:

- **Optimal alignment algoritms**, which find the best alignment but only for small to medium sized sequences.
- **_Dot Plot_ algorithms**, which find similar regions in big sequences.
- **_K-word_ algorithms**, which search for similar sequences on huge databases of sequences.

On the upcoming posts we will explore each one of these types of algorithms in detail [[4]](#note4).

# Next Steps

I hope I could wet you appetite for this subject. In the next post we will formally define alignments and discover how to quantify how good they are.


# Notes

__<a name="note1">[1]</a>__ Note that whenever we find the same nucleotide in a similar (homologous) position in two sequences we cannot exclude the possibility of the occurrence of reversible mutations, i.e., a sequence of mutations that, over time, made both positions the same again:$$G \rightarrow T \rightarrow G$$.

__<a name="note2">[2]</a>__ In reality we cannot distinguish, from sequence alone, if the historical mutation was an insertion in one sequence or a deletion in the other sequence.

__<a name="note3">[3]</a>__ Optimal alignment algorithms have a time and memory complexity of $$\mathcal{O}(n \times m)$$.

__<a name="note4">[4]</a>__ It is important to mention that in this series of posts I am only considering pairwise sequence alignments. Multiple sequence alignment algorithms, like [Clustal](http://www.clustal.org/), also exist and are pretty useful, but we will not explore them,at least for now.


# References

- <a name="Pray2008"/>Pray, L. [__DNA Replication and Causes of Mutation__](https://www.nature.com/scitable/topicpage/dna-replication-and-causes-of-mutation-409/), 2008, Nature Education 1 (1):214

- <a name="Kimura1983"/>Kimura, M. [__The neutral theory of molecular evolution__](https://en.wikipedia.org/wiki/Neutral_theory_of_molecular_evolution), 1983, Cambridge University Press.

- <a name="Yang2006"/>Yang Z. [__Computational Molecular Evolution__](https://academic.oup.com/book/34842), 2006, Oxford University Press.

- <a name="Uffelmann2021"/>Uffelmann E, Huang QQ, Munung NS, de Vries J, Okada Y, Martin AR, Martin HC, Lappalainen T and Posthuma D, 2021, [__Genome-wide association studies__](https://www.nature.com/articles/s43586-021-00056-9), Nature Reviews Methods Primers volume 1, 59.