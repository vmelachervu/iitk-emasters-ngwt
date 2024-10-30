from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score

DiabetesData = pd.read_csv('Diabetes.csv')

X = DiabetesData.iloc[:, [0, 1]].values
Y = DiabetesData.iloc[:, 2].values


