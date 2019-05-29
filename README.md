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

The input JSON should contain either five entries or six entries. If it contains
five entries, the bars of the resulting plot will have a 2-2-1 style. If it
contains six entries, the bars will have a 2-1-2-1 style.

Significance indicators will be written above the bars: `***` if p<0.001,
`**` if p<0.01, `*` if p<0.05, `ns` otherwise.

Here is an example of a plot in the 2-1-2-1 style:

![example barplot](https://github.com/anthony-aylward/islet-cytokines-outline/raw/master/figure/rs3787186_luc/dex_vs_untreated.png)
