# Greta Dynamic Factor Analysis:
#

#### Source:
# https://mdscheuerell.github.io/gretaDFA/



####  Packages  ####
library(gmRi)
library(targets)
library(tidyverse)
library(greta)


## R pkgs
installed.packages()[c("greta", "tensorflow"), "Version"]

## tensorflow
# tensorflow::install_tensorflow()
tensorflow::tf_version()


####  Process  ####