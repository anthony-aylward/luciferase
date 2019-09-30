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
![example ratio plot](https://github.com/anthony-aylward/luciferase/raw/master/example/ratio.png)

## Meta-analysis of luciferase reporter data

