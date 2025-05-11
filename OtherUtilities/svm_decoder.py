import sklearn
from sklearn import svm
import numpy as np


def svm_decoder(array1, array2, dt=None, sample_num=50, train_ratio=0.8):
    """
    array1: neurons x trials x timebins
    array2: neurons x trials x timebins
    sample_num: the number of sampling times for bootstrap
    dt: is the length of sliding window to calculate the accuracy across time.
        When dt is None, the decoding is based on the whole time.
    train_ratio: the ratio of train data
   """
   
    # check the min number of trials
    min_trial_num = np.min([array1.shape[1], array2.shape[1]])
    print('min_trial_num: ', min_trial_num)
    
    accuracy_across_sampels = []

    for samplei in range(sample_num):
        
        # random select the same number of trials
        trial_index1 = np.random.choice(array1.shape[1], min_trial_num, replace=False)
        trial_index2 = np.random.choice(array2.shape[1], min_trial_num, replace=False)
        sample_array1 = array1[:, trial_index1, :] 
        sample_array2 = array2[:, trial_index2, :] 
        
        # Total time 
        T = array1.shape[2] # time bins
        assert T == array2.shape[2], 'The time bins of two arrays are not equal'
        
        # check the length of sliding window
        if dt is None:
            dt = T
        else:
            assert dt<=T, 'The sliding window is larger than the total time'
        
        accuracy_across_time = []

        for intervali in range(np.int32(T/dt)):
            X1 = []
            y1 = []
            X2 = []
            y2 = []
            for triali in range(min_trial_num):
                X1.append(sample_array1[:,triali, intervali*dt:(intervali+1)*dt].reshape(-1))
                y1.append(0)
                X2.append(sample_array2[:,triali, intervali*dt:(intervali+1)*dt].reshape(-1))
                y2.append(1)
            X = X1 + X2
            y = y1 + y2
            
            # shuffle the data
            random_index = np.random.permutation(len(X))
            random_X = np.asarray(X)[random_index]
            random_y = np.asarray(y)[random_index]
            
            # split the data into train and test
            train_X = random_X[:int(train_ratio*len(X))]
            train_y = random_y[:int(train_ratio*len(X))]
            test_X = random_X[int(train_ratio*len(X)):]
            test_y = random_y[int(train_ratio*len(X)):]
            
            # train the svm
            clf = svm.SVC()
            clf.fit(train_X, train_y)
            
            # test the svm
            test_y_predict = clf.predict(test_X)
            
            # calculate the accuracy
            accuracy_intervali = sum(test_y_predict == test_y)/float(len(test_y))
            accuracy_across_time.append(accuracy_intervali)
            
        accuracy_across_sampels.append(accuracy_across_time)
        
    return accuracy_across_sampels

