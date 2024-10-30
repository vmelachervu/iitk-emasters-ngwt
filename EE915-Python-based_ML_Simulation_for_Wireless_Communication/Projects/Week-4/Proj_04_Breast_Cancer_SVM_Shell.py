from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC


bcancer = datasets.load_breast_cancer()

X = bcancer.data
Y = bcancer.target

