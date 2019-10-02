# luciferase
Tool for plotting luciferase reporter data. Thanks due to Joshua Chiou and Mei-lin Okino for inspiration and contributions.

## Installation
```sh
pip3 install luciferase
```
or
```sh
pip3 install --user luciferase
```

## Command-line interface for barplots

### Barplots of enhancer activity

A script called `luciferase-barplot` for creating bar plots from JSON-formatted
data is included. After installing `luciferase`, you can use it like this:
```sh
luciferase-barplot --title "plot title" example.json example.pdf 
```

See also the help message:
```sh
luciferase-barplot -h
```

Examples of luciferase reporter data in JSON format:
```json
{
  "Non-risk, Fwd": [8.354, 12.725, 8.506],
  "Risk, Fwd": [5.078, 5.038, 5.661],
  "Non-risk, Rev": [9.564, 9.692, 12.622],
  "Risk, Rev": [10.777, 11.389, 10.598],
  "Empty": [1.042, 0.92, 1.042]
}
```
```json
{
  "Alt, MIN6": [5.47, 7.17, 6.15],
  "Ref, MIN6": [3.16, 3.04, 4.34],
  "Empty, MIN6": [1.07, 0.83, 0.76],
  "Alt, ALPHA-TC6": [2.50, 3.47, 3.33],
  "Ref, ALPHA-TC6": [2.01, 1.96, 2.31],
  "Empty, ALPHA-TC6": [1.042, 0.92, 1.042]
}
```

The input JSON should contain either five, six, or twelve entries. If it
contains five entries, the bars of the resulting plot will have a 2-2-1 style.
If it contains six entries, the bars will have a 2-1-2-1 style. If twelve,
the syle will be as with six entries but doubled.

Significance indicators will be written above the bars: `***` if p<0.001,
`**` if p<0.01, `*` if p<0.05, `ns` otherwise.

Here is an example of a plot in the 2-1-2-1 style:

![example barplot](https://github.com/anthony-aylward/islet-cytokines-outline/raw/master/figure/rs3787186_luc/dex_vs_untreated.png)


### Barplots of allelic ratio

A second script called `luciferase-ratioplot` takes the same input data and
produces a comparative plot of allelic ratios:

```sh
luciferase-ratioplot --title "plot title" example.json example.pdf
```

For this script, the number of entries in the input JSON should be a multiple
of 3. The resulting plot shows the estimated allelic ratio of enhancer activity
with confidence intervals (95% by default). Here is an example input dataset
and plot:

```
{
  "Alt, dex": [44.6, 37.6, 37.7],
  "Ref, dex": [149.4, 99.7, 124.5],
  "Empty, dex": [1.1, 1.0, 0.9],
  "Alt, untreated": [19.7, 16.2, 18.3],
  "Ref, untreated": [33.2, 30.3, 33.3],
  "Empty, untreated": [1.0, 1.0, 1.0]
}
```
![example ratio plot](example/ratio.png)

## Meta-analysis

It may be that we have performed two or more experiments
(from separate minipreps) and wish to meta-analyze the results. As an example,
let's consider the results of two identical experiments on a regulatory
variant at the _SIX3_ locus: [SIX3-MP0](example/six3-mp0.json) and [SIX3-MP1](example/six3-mp1.json). First we'll plot both datasets separately:
```sh
luciferase-barplot six3-mp0.json six3-mp0.png --light-colors '#DECBE4' '#FED9A6' '#FBB4AE' --dark-colors '#984EA3' '#FF7F00' '#E41A1C' --title 'SIX3-MP0'
luciferase-barplot six3-mp1.json six3-mp1.png --light-colors '#DECBE4' '#FED9A6' '#FBB4AE' --dark-colors '#984EA3' '#FF7F00' '#E41A1C' --title 'SIX3-MP1'
```

<img src="example/six3-mp0.png" width="400"/><img src="example/six3-mp1.png" width="400"/>

We can see that the results are fairly consistent in character, but checking
the y-axis tells us that they are on different scales. Intuitively, we might
conclude from these results that there are allelic effects under all three
conditions. Ideally though, we would like to use all of the data at once for
one plot to get the most accurate conclusions about allelic effects.

We might simply combine the data into one dataset, (as [here](example/six3-meta-nobatch.json)) and plot it:
```sh
luciferase-barplot six3-meta-nobatch.json six3-meta-nobatch.png --light-colors '#DECBE4' '#FED9A6' '#FBB4AE' --dark-colors '#984EA3' '#FF7F00' '#E41A1C'
```

![meta-analysis without batch](example/six3-meta-nobatch.png)

The bar heights look reasonable, and the allelic effects appear clear from
looking at them, but all of the hypothesis tests returned non-significant
results. What gives?

The answer is that combining data from experiments with different scales
breaks the assumptions of the significance test (a t-test). To meta-analyze
these data in a useful way, we first need to re-normalize the two experiments
to put both of them on the same scale. `luciferase-barplot` will re-normalize
the data automatically if the dataset includes an additional entry ("Batch")
indicating the batch of each data point, as in this example:
[SIX3-META](example/six3-meta.json).
```json
{
  "Alt, untreated": [19.7, 16.2, 18.3, 6.5, 8.0, 4.4],
  "Ref, untreated": [33.2, 30.3, 33.3, 8.4, 13.6, 17.1],
  "Empty, untreated": [1.0, 1.0, 1.0, 1.1, 1.0, 0.9],
  "Alt, hi_cyt_noTNFA": [11.0, 8.8, 10.1, 3.2, 3.7, 3.3],
  "Ref, hi_cyt_noTNFA": [17.1, 16.7, 18.8, 7.6, 6.7, 5.5],
  "Empty, hi_cyt_noTNFA": [1.1, 0.9, 1.0, 1.1, 0.9, 1.0],
  "Alt, hi_cyt": [10.8, 10.9, 9.1, 3.1, 2.7, 4.0],
  "Ref, hi_cyt": [17.8, 16.1, 18.0, 7.7, 7.0, 7.1],
  "Empty, hi_cyt": [1.0, 1.0, 1.0, 1.0, 1.0, 1.1],
  "Batch": [0, 0, 0, 1, 1, 1]
}
```

Here is what the results look like when they're re-normalized to correct for
batch
```sh
luciferase-barplot six3-meta.json six3-meta.png --light-colors '#DECBE4' '#FED9A6' '#FBB4AE' --dark-colors '#984EA3' '#FF7F00' '#E41A1C' --title 'SIX3-META'
```

![meta-analysis with batch](example/six3-meta.png)
