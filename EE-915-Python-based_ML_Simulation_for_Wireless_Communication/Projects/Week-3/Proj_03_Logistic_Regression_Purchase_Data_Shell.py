
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np


def sigmoid(x):
  return 1/(1 + np.exp(-x))

purchaseData = pd.read_csv('Purchase_Logistic.csv')

X = purchaseData.iloc[:, [2, 3]].values
Y = purchaseData.iloc[:, 4].values


