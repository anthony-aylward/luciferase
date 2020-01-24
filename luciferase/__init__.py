"""Helper functions and scripts for luciferase reporter data

Functions
---------
luciferase_barplot
    create a barplot from luciferase reporter data

Examples
--------
import luciferase
luc_data = {
    'Non-risk, Fwd': [8.354, 12.725, 8.506],
    'Risk, Fwd': [5.078, 5.038, 5.661],
    'Non-risk, Rev': [9.564, 9.692, 12.622], 
    'Risk, Rev': [10.777, 11.389, 10.598],
    'Empty': [1.042, 0.92, 1.042]
}
luciferase.luciferase_barplot(luc_data, 'rs7795896.pdf', title='rs7795896')
luc_data = {
    'Alt, MIN6': [5.47, 7.17, 6.15],
    'Ref, MIN6': [3.16, 3.04, 4.34],
    'Empty, MIN6': [1.07, 0.83, 0.76],
    'Alt, ALPHA-TC6': [2.50, 3.47, 3.33],
    'Ref, ALPHA-TC6': [2.01, 1.96, 2.31],
    'Empty, ALPHA-TC6': [1.042, 0.92, 1.042]
}
luciferase.luciferase_barplot(luc_data, 'min6-v-alpha.pdf', title='MIN6 v.Alpha')
luciferase.luciferase_ratioplot(luc_data, 'min6-v-alpha-ratio.pdf', title='MIN6 v.Alpha')
"""

from luciferase.luciferase import luciferase_barplot
from luciferase.ratioplot import luciferase_ratioplot
from luciferase.swarmplot import luciferase_swarmplot
