This file is the documentation for the following paper

Title: Density dependence of species interactions and ecological stability

Authors: Kazutaka Kawatsu, Yutaka Osada, Yumiko Ishii, Reiji Masuda, Masakazu Shimada and Michio Kondoh

Written on 20211215 by K.Kawatsu.

Last update: 20220917.

'data' directory includes dataset for the analysis and the result

'fig' directory includes figure files appeared in the manuscript

'R' directory includes R source codes for the analysis and visualization

'rpkg' directory is an R package for this paper, which includes R and cpp source codes

Installation of rpkg
--------------------

R codes of this manuscript depends on the specific R package 'rpkg'.
rpkg contains cpp and header files, thus it must be compiled for your environment before usage.

Installation procedure:

if you use macOS or Linux environments, start a terminal app and execute the following command

    cd "path of your current workspace, which includes 'rpkg'"
    R CMD INSTALL rpkg

Or if you use Windows PC, run the following command in R console

    setwd("path of your current workspace, which includes 'rpkg'")
    devtools::check()
    install.packages("rpkg", repos = NULL, type = "source")

Description of rpkg
-------------------

rpkg contains simple R help documentation for each function.
Please refer to them for usage of each function

For example,

    ?rpkg::find_best_dim
    help(rpkg::find_best_dim)

