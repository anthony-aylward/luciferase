# luciferase
Helper functions and scripts for luciferase data

## Installation
```
pip3 install luciferase
```

## Examples

A script called `reporter-barplot` for creating bar plots from JSON-formatted data is included. After installing `luciferase`, you can use it like this:
```
reporter-barplot --title "plot title" example.json example.pdf 
```

See also the help message:
```
reporter-barplot -h
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
