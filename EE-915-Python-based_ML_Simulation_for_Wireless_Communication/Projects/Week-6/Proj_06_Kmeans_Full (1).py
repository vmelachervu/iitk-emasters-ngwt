import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs



X, y = make_blobs(n_samples=2500,centers=4, n_features=2,random_state = 10)

plt.figure()
plt.scatter(X[:, 0], X[:, 1], c=None, cmap='jet',s=10)
plt.suptitle('Original Data')
plt.grid(1,which='both')
plt.axis('tight')
plt.show()



