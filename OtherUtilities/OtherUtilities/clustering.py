import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D  
import seaborn as sns


def cal_sse_and_silhouette_score(data, n_clusters, random_state=42, plot=True):
    """
    calculate the silhouette score for the given data and number of clusters
    data: 2D array, shape (n_samples, n_features), normally has a shape as neurons x timebins
    n_clusters: int, number of clusters to try
    random_state: int, random state for reproducibility
    plot: bool, whether to plot the results
    """
    sse_scores = []
    silhouette_scores = []
    
    for k in range(2, n_clusters+1):
        kmeans = KMeans(n_clusters=k, random_state=random_state)
        kmeans.fit(data)
        sse_scores.append(kmeans.inertia_)  # SSE
        silhouette_scores.append(silhouette_score(data, kmeans.labels_))  # Silhouette score

    if plot:
        plt.figure(figsize=(6, 3))
        plt.subplot(1, 2, 1)
        plt.plot(np.arange(2,n_clusters+1,1), sse_scores, marker='o')
        plt.xlabel('Number of Clusters')
        plt.ylabel('SSE')
        plt.title('Elbow Method for Optimal K')
        plt.subplot(1, 2, 2)
        plt.plot(np.arange(2,n_clusters+1,1), silhouette_scores, marker='o')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Silhouette Score')
        plt.title('Silhouette Method for Optimal K')
        plt.tight_layout()
        plt.show()
        
    return sse_scores, silhouette_scores
        
        
def kmeans_clustering(data, n_clusters=5, random_state=42):
    """
    Perform KMeans clustering on the data and return the labels
    data: 2D array, shape (n_samples, n_features), normally has a shape as neurons x timebins
    n_clusters: int, number of clusters to use
    random_state: int, random state for reproducibility
    """
    model = KMeans(n_clusters=5, random_state=42)
    labels = model.fit_predict(data)
    
    return model, labels

def plot_number_of_neurons_in_clusters(labels, pal=sns.color_palette('deep', 11)):
    """
    Plot the number of neurons in each cluster
    labels: 1D array, shape (n_samples,), the cluster labels of the data
    pal: list, color palette for the clusters
    """
    plt.figure(figsize=(3,3))
    plt.pie(np.unique(labels, return_counts=True)[1], labels=np.unique(labels, return_counts=True)[0], autopct='%1.1f%%', colors=pal)
    plt.title('Number of Neurons in Clusters')
    
    return plt.show()


def get_lowdimensional_data(data, labels, n_dim=3, pal=sns.color_palette('deep', 11), plot=True):
    """
    Visualize the low-dimensional data using t-SNE
    data: 2D array, shape (n_samples, n_features), normally has a shape as neurons x timebins
    labels: 1D array, shape (n_samples,), the cluster labels of the data
    n_dim: int, number of dimensions to reduce to (2 or 3)
    pal: list, color palette for the clusters
    """

    # reduce the dimensionality of the data using t-SNE
    tsne = TSNE(n_components=n_dim, perplexity=10, random_state=42)
    tsne_data = tsne.fit_transform(data)

    if plot:
        fig = plt.figure(figsize=(3, 3))
        if n_dim==3:
            ax = fig.add_subplot(111, projection='3d')
        for i in range(np.unique(labels).shape[0]):
            ax.scatter(tsne_data[labels == i, 0], tsne_data[labels == i, 1], tsne_data[labels == i, 2], s=10, color=pal[i], label='Cluster {}'.format(i))
        # ax.legend(frameon=False)
        ax.set_xlabel('t-SNE Dimension 1')
        ax.set_ylabel('t-SNE Dimension 2')
        ax.set_title('t-SNE Visualization')
        ax.grid(False)
        plt.show()
    
    return tsne_data


