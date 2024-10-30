import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score


purchaseData = pd.read_csv('Purchase_Logistic.csv')

X = purchaseData.iloc[:, [2, 3]].values
Y = purchaseData.iloc[:, 4].values

