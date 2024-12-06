# Class 6: R functions
Mudit

My first function

``` r
add <- function(x, y=1, z = 0){
  x + y
}
```

Can I just use it?

``` r
add(1, 1)
```

    [1] 2

``` r
add(x=1, y=100)
```

    [1] 101

``` r
add(c(100, 1, 100), 1)
```

    [1] 101   2 101

``` r
add(10)
```

    [1] 11

``` r
add(1, 1, 1)
```

    [1] 2

``` r
GENERATE_DNA <- function(l = 1){
  
  
  bases <- c('A', 'T', 'G', 'C')
  sequence1 <- sample(bases, l, replace = TRUE)
  #paste0(sample(bases, l, replace = TRUE), collapse = "")
  sequence1 <- paste0(sample(bases, l, replace = TRUE), collapse = "")
  return(sequence1)
  
}
```

``` r
GENERATE_DNA(100)
```

    [1] "ATCGTCTGCGTAACGGGAAGTTTCGAGATTCTGGCGAGTGCTAGCGATTAAAGTACTGCAGTGGCGGAGTTGTGCCATCTAGACTTCCGCGGTTGTTAAT"

``` r
GENERATE_PROTEIN <- function(l = 1){
  
  amino <- unique(bio3d::aa.table$aa1)
  bases <- c('A', 'T', 'G', 'C')
  #sample(bases, l, replace = TRUE)
  #paste0(sample(bases, l, replace = TRUE), collapse = "")
  sequence1 <- paste0(sample(amino, l, replace = TRUE), collapse = "")
  return(sequence1)
  
}
```

``` r
GENERATE_DNA(100)
```

    [1] "GCACGACCAGCCTGCATATGGCCGGGCTCCTAGGTAGGGTGGGCATTCTCTCTTGACGTATGCCATGTCCTCCGCAAAACCACAAACGATCTACTGTTAC"

``` r
unique(bio3d::aa.table$aa1)[1:20]
```

     [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y"
    [20] "V"

``` r
View(bio3d::aa.table)
```

Generate random protein sequences of length 6 to 12

``` r
answer<-sapply(6:12, GENERATE_PROTEIN)
cat(paste(">id.", 6:12, "\n", answer, sep=''), sep="\n")
```

    >id.6
    IYDNCQ
    >id.7
    GPDAWDG
    >id.8
    SFTITCVN
    >id.9
    YXHFXPYDL
    >id.10
    GDQRKQTLFI
    >id.11
    MHPRQDKLSDD
    >id.12
    LQPGQRYVMFIA

``` r
#grade() <- function(x = c(1, 1, 1, 1)){
#  x<-sort(x)
  
#}
```

``` r
df <- read.csv("student_homework.csv")
```
