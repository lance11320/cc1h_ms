import scipy
import scipy.stats as stats
import numpy as np

def pair_test(data1, data2):
    """
    Perform paired t test or wilcoxon test on two paired datasets.
    data1: 1D array, shape (n_samples,)
    data2: 1D array, shape (n_samples,)
    """
    assert data1.shape == data2.shape, 'The shape of two data sets are not equal'
    assert len(data1.shape) == 1, 'The shape of data sets should be 1D'
    assert len(data2.shape) == 1, 'The shape of data sets should be 1D'
    assert len(data1) == len(data2), 'The length of two data sets are not equal'
    assert len(data1) >= 5, 'The length of data sets should be larger than 5'
    assert len(data2) >= 5, 'The length of data sets should be larger than 5'
    
    if len(data1)>=8 and len(data2)>=8:
        print('normal distribution test')
        _, p_normal1 = stats.normaltest(data1)
        print('normal distribution p_value of data1:', p_normal1)
        _, p_normal2 = stats.normaltest(data2)
        print('normal distribution p_value of data2:', p_normal2)
    else:
        print('normal distribution test')
        _, p_normal1 = stats.shapiro(data1)
        print('normal distribution p_value of data1:', p_normal1)
        _, p_normal2 = stats.shapiro(data2)
        print('normal distribution p_value of data2:', p_normal2)

    if p_normal1<0.05 or p_normal2<0.05:
        print('wilcoxon test')
        _, p_value = stats.wilcoxon(data1, data2, zero_method='zsplit')       
             
    elif p_normal1>=0.05 and p_normal2>=0.05:
        # use paired t test
        print('paired t test')
        _, p_value = stats.ttest_rel(data1, data2)
    
    return p_value

def independent_test(data1, data2, equal_variance=True):
    """
    Perform t test or mannwhitneyu test on two independent datasets.
    data1: 1D array, shape (n_samples,)
    data2: 1D array, shape (n_samples,)
    """
    assert len(data1.shape) == 1, 'The shape of data sets should be 1D'
    assert len(data2.shape) == 1, 'The shape of data sets should be 1D'
    assert len(data1) >= 5, 'The length of data sets should be larger than 5'
    assert len(data2) >= 5, 'The length of data sets should be larger than 5'
    
    if len(data1)>=8 and len(data2)>=8:
        print('normal distribution test')
        _, p_normal1 = stats.normaltest(data1)
        print('normal distribution p_value of data1:', p_normal1)
        _, p_normal2 = stats.normaltest(data2)
        print('normal distribution p_value of data2:', p_normal2)
    else:
        print('normal distribution test')
        _, p_normal1 = stats.shapiro(data1)
        print('normal distribution p_value of data1:', p_normal1)
        _, p_normal2 = stats.shapiro(data2)
        print('normal distribution p_value of data2:', p_normal2)

    if p_normal1<0.05 or p_normal2<0.05:
        print('mannwhitneyu test')
        _, p_value = stats.mannwhitneyu(data1, data2, zero_method='zsplit')       
             
    elif p_normal1>=0.05 and p_normal2>=0.05:
        # use  t test
        print(' t test')
        _, p_value = stats.ttest_ind(data1, data2, equal_var=equal_variance)
    
    return p_value


def cohens_d(a, b):
    # Cohen's d
    n1, n2 = len(a), len(b)
    # calculate pooled SD
    pooled_sd = np.sqrt(((n1 - 1)*a.var(ddof=1) + (n2 - 1)*b.var(ddof=1)) / (n1 + n2 - 2))
    return (a.mean() - b.mean()) / pooled_sd


def calculate_CI(data1, data2, ci=0.95):
    df = len(data1) + len(data2) - 2
    mean_diff = data1.mean() - data2.mean()
    pooled_se = np.sqrt(data1.var(ddof=1)/len(data1) + data2.var(ddof=1)/len(data2))
    ci_low, ci_high = stats.t.interval(
        ci, df=df, loc=mean_diff, scale=pooled_se)
    return ci_low, ci_high


def permutation_test (a, b, n_perm = 10000):
    # for small samples
    rng = np.random.default_rng()
    obs_diff = a.mean() - b.mean()
    combined = np.hstack([a, b])
    count = 0
    for _ in range(n_perm):
        rng.shuffle(combined)
        new_a = combined[:len(a)]
        new_b = combined[len(a):]
        if abs(new_a.mean() - new_b.mean()) >= abs(obs_diff):
            count += 1
    p_perm = (count + 1) / (n_perm + 1) 
    print('Permutation test p-value:', p_perm)
    return p_perm
