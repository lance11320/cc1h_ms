
import scipy.stats as stats
from scipy.stats import wilcoxon

def statistic_test(data1, data2):
    """
    Perform paired t test or wilcoxon test on two data sets.
    data1: 1D array, shape (n_samples,)
    data2: 1D array, shape (n_samples,)
    """
    assert data1.shape == data2.shape, 'The shape of two data sets are not equal'
    assert len(data1.shape) == 1, 'The shape of data sets should be 1D'
    assert len(data2.shape) == 1, 'The shape of data sets should be 1D'
    assert len(data1) == len(data2), 'The length of two data sets are not equal'
    assert len(data1) >= 5, 'The length of data sets should be larger than 5'
    assert len(data2) >= 5, 'The length of data sets should be larger than 5'
    
    print('normal distribution test')
    _, p_normal1 = stats.normaltest(data1)
    print('normal distribution p_value of data1:', p_normal1)
    _, p_normal2 = stats.normaltest(data2)
    print('normal distribution p_value of data2:', p_normal2)

    if p_normal1<0.05 or p_normal2<0.05:
        print('wilcoxon test')
        _, p_value = wilcoxon(dat1, data2, zero_method='zsplit')       
             
    elif p_normal1>=0.05 and p_normal2>=0.05:
        # use paired t test
        print('paired t test')
        _, p_value = stats.ttest_rel(data1, data2)
    
    return p_value