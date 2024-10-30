from sklearn import datasets
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.metrics import confusion_matrix



irisset = datasets.load_iris()
X = irisset.data
Y = irisset.target

pca = PCA(n_components=2)
Xp = pca.fit(X).transform(X)

GMM = GaussianMixture(n_components = 3)
GMM.fit(Xp)
y_predG = GMM.predict(Xp)

cmat = confusion_matrix(Y, y_predG)

plt.figure(1)
plt.scatter(Xp[:, 0], Xp[:, 1], c=Y, cmap='jet',s=10)
plt.suptitle('Original Clusters')
plt.grid(1,which='both')
plt.axis('tight')
plt.show()



plt.figure(2)
plt.scatter(Xp[:, 0], Xp[:, 1], c=y_predG, cmap='jet',s=10)
plt.suptitle('GMM Clusters')
plt.grid(1,which='both')
plt.axis('tight')
plt.show()

