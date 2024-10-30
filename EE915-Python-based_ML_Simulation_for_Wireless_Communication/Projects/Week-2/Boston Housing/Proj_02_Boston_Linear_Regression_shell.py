import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import pandas as pd

BosData = pd.read_csv('BostonHousing.csv')
X = BosData.iloc[:,0:11]
y = BosData.iloc[:, 13] # MEDV: Median value of owner-occupied homes in $1000s

# Write code here


print('Train RMSE =', rmse)
print('Train R2 score =', r2)
print("\n")

# Write code here


print('Test RMSE =', rmse)
print('Test R2 score =', r2)
