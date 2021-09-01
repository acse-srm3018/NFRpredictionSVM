# Auto-characterization NFR using multi-output least squares support vector regression

[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


**Note:**
This project was published in [Arabian Journal of Geosciences]((https://link.springer.com/article/10.1007/s12517-021-06559-9)).

## Project description
Horizontal wells sense multiple points without need to drill additional vertical wells became economically viable after 1980s. Characterization of this system is required for estimation their production ability, management their behaviors, and designing enhanced recovery scenarios. Well testing signals are the most valuable source of information for characterization of hydrocarbon reservoirs having horizontal wells. Estimation of reservoir parameters as well as extent of reliability of this estimation directly depends on correct detection of reservoir interpretation model.  Indeed, it is traditional to detect reservoir model using well testing analysis prior to starting parameter estimation.
The multi-layer perceptron (MLP) neural network have widely used for pattern recognition in different scientific and engineering branches up to now. Therefore, different topologies for the MLP paradigms were examined for identifying the reservoir models from the pressure derivative curves. 2560 of pressure derivative graphs for six reservoir models are utilized for designing the most reliable MLP neural network model.


Considering some statistical and computational criteria, the multi-layer perceptron (MLP) neural network model with five hidden neurons trained with the scaled conjugate gradient are chosen as a best intelligent model for the considered task. This model provides total classification accuracy (TCA) of 98.3%, mean square errors (MSE) of 0.00725, and regression coefficient (R2) of 0.97332. Finally the developed model has been examined using real field data, noisy data and some data sets outside the range of train and test data. In all cases, the network was able to detect the correct reservoir model with the probability near to 90%.

## Data generation
The data was generated using PanSystem and then we used a QTBindCode using C to automatically generated data (pressure data). The codes can be found [here](https://github.com/acse-srm3018/ReservoirClassification/tree/main/GenerationData). 

After that we used a pressure derivative function to take derivative.


## Dataset
You can download dataset from [Excel file](https://github.com/acse-srm3018/NFRpredictionSVM/blob/main/Dataset.xlsx).


The dataset includes 500 data within 6 clases that will be used for training ANN models.

## Basic information

An overview of the files is provided below.
- `src_python/` contains all python code files
- `src_matlab/` contains all matlab code files
- `images/` contains all images used in the article.
- `dataset.xlsx` all generated pressure data that use as training and validation dataset.
- `AJGS.pdf` article in pdf format.
- `LICENSE.txt` is the MIT license.
- `README.md` contains basic information for the repository and detailed information for how to compile and reproduce the results.


## Installation

you can clone and open directories using

```
git clone https://github.com/acse-srm3018/NFRpredictionSVM
```

## Documentation

 The articles published and can be found [here](https://github.com/acse-srm3018/NFRpredictionSVM/blob/main/AJGS.pdf).

## References


