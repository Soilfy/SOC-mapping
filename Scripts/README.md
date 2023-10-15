# Soilify SOC Scripts

This repository contains scripts written for the **Soilify SOC Project**. The scripts are developed using Python and JavaScript, and they utilize various libraries such as Pandas, NumPy, and Scikit-learn.

## Data Preparation Script

The data preparation script is named `soc_data_prep.ipynb`. This Jupyter notebook contains code for cleaning and preprocessing raw data. It reads the raw data from a CSV file, performs data cleaning and preprocessing, and saves the cleaned data to a new CSV file. The notebook also includes exploratory data analysis (EDA) to help gain insights into the data.

## Modeling Script

The modelling script is named `soc_models.ipynb`. This Jupyter notebook contains code for building and evaluating machine learning models. It reads the cleaned data from the CSV file generated by the data preparation script, splits the data into training and testing sets, and trains various machine learning models on the training data. The notebook then evaluates the performance of the models on the testing data and selects the best model based on evaluation metrics.

## Notebooks Description

### osc_data_prep.ipynb

This notebook is organized into the following sections:

1. **Introduction**: Provides a brief introduction to the notebook and the data.

2. **Data Loading**: Contains code for loading raw data from a CSV file.

3. **Data Cleaning**: Includes code for data cleaning, such as handling missing values, removing duplicates, and correcting data types.

4. **Data Preprocessing**: Provides code for preprocessing the data, including feature engineering and scaling.

5. **Exploratory Data Analysis**: Contains code for visualizing and exploring the data.

6. **Data Saving**: Includes code for saving the cleaned data to a new CSV file.

### soc_models.ipynb

This notebook is organized into the following sections:

1. **Introduction**: Provides a brief introduction to the notebook and the data.

2. **Data Loading**: Contains code for loading the cleaned data from the CSV file generated by the data preparation script.

3. **Data Splitting**: Includes code for splitting the data into training and testing sets.

4. **Model Training**: Provides code for training several machine learning models on the training data.

5. **Model Evaluation**: Contains code for evaluating the performance of the models on the testing data.

6. **Model Selection**: Includes code for selecting the best model based on evaluation metrics.