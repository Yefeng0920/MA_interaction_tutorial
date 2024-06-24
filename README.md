# Methodological guide: Meta-analytical approaches for analysing stressor effects (or treatment effects in general) and interaction

## A step-by-step illustration (v 1.0)

## Last update: June 2024


# Reproducibility claim

The tutorial was developed by Yefeng Yang, and can be independently reproduced by Coralie Williams and Shinichi Nakagawa. All version information involved in the tutorial are printed at the end (subsection Software and package versions).

# Preamble

Using meta-analytical approaches for modelling interactive effects is in its infancy, there is a risk that the statistical analysis issues in leading journals as we identified in a paper published in Nature Communications could become widely repeated in the future. See also issues we identified in a recent Nature paper.

> Hu, N., Bourdeau, P. E., & Hollander, J. (2024). Responses of marine trophic levels to the combined effects of ocean acidification and warming. Nature Communications, 15(1), 3400.

> Siviter, H., Bailes, E. J., Martin, C. D., Oliver, T. R., Koricheva, J., Leadbeater, E., & Brown, M. J. (2023). Addendum: Agrochemicals interact synergistically to increase bee mortality. Nature, 617(7960), E7-E9.

Therefore, we provide a step-by-step tutorial with analytical scripts on rigorously analyzing multiple stressor effects (or, treatment effects in general) and their interactions using meta-analytical approaches.

# Credit

If our paper and tutorial have helped you, please cite the following paper:

> Yefeng Yang, Coralie Williams, Kyle Morrison, Jacob Bishop, Jinming Pan, Malgorzata Lagisz, Shinichi Nakagawa. Statistical issues in multiple stressors effects and their interactions. EcoEvoRxiv, 2024.


## Structure

The repository contains 3 folders:

- `dat`

- `function`

- `script`

### `dat` folder

`OAW_interaction.csv`: example dataset comprising 486 observations from 162 fully factorial experiments examining the two main effects, and their interaction.

### `function` folder

Various computing and visualizing help functions.

### `script` folder

`tutorial.Rmd`: a step-by-step tutorial containing `R` scripts in the format of `.Rmd`.

`tutorial.html`: a step-by-step tutorial containing `R` scripts and interactive results in the format of `.html`.


## Licence

The files in this dataset are licensed under the Creative Commons Attribution 4.0 International License (to view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/).


## Correspondence regarding the analysis and code

- Dr. Yefeng Yang

Evolution & Ecology Research Centre (EERC), 
School of Biological, Earth and Environmental Sciences (BEES), 
The University of New South Wales, Sydney, Australia

Email: yefeng.yang1@unsw.edu.au

