Raw code:

```{r}
library("GenomicRanges")
exons <- GRanges("seq", IRanges(start = c(1, 1, 13), end = c(5, 8, 15)))
exons

## Results in 2 reduced exons. Cannot get the counts for exons 1 or 2.
reduce(exons)

## Results in 3 disjoint exons. The sum of disjoint exon 1 and 2 is equal to exon 2.
disjoin(exons)
```

Output:

```{r}
exons
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]      seq  [ 1,  5]      *
  [2]      seq  [ 1,  8]      *
  [3]      seq  [13, 15]      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Results in 2 reduced exons. Cannot get the counts for exons 1 or 2.
reduce(exons)
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]      seq  [ 1,  8]      *
  [2]      seq  [13, 15]      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Results in 3 disjoint exons. The sum of disjoint exon 1 and 2 is equal to exon 2.
disjoin(exons)
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]      seq  [ 1,  5]      *
  [2]      seq  [ 6,  8]      *
  [3]      seq  [13, 15]      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
