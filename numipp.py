#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:13:09 2023

@author: yitianli

IPP-GMM code with OOP
"""

import numpy as np
from scipy.stats import randint

class PanelModel:
    '''
    A panel data model with fixed effects.
    
    ...
    
    Attributes
    -------
    z_F1 : numpy array 
        z_F1 is the dataset. It has dimension (n, k+1, T) where
            n is the number of units, 
            T is the number of periods, 
            k is the number of regressors.
        For example, let y be the regressand with dimension (n, 1, T) and x be 
        the regressor with dimension (n, k, T). Then z_F1 can be constructed using
        z = np.concatenate((y, x), axis=1).
    estimate_theta : function
        A function for estimating theta, the parameter of interest. It takes a dataset,
        e.g., z_F1, as parameters and returns the estimate of theta, in the format of 
        a numpy array.
    
    Methods
    -------
    comp_theta_jack():
        Corrects the IPP asymptotic bias using the delete-one jackkinfe. 
        It returns a numpy array consists the first and second order bias-corrected 
        estimator.
    comp_theta_boot():
        Corrects the IPP asymptotic bias using the delete-one jackkinfe.
        It returns a numpy array consists the first, second, and third order 
        bias-corrected estimator.
    
    '''
    
    def __init__(self, z_F1, estimate_theta):
        self.z_F1 = z_F1    # z_F1 is the original data with dimension (n,k,T)
        self.estimate_theta = estimate_theta
        self.theta_hat = self.estimate_theta(z_F1)
        
    def comp_theta_jack(self):
        '''
        This function corrects the IPP asymptotic bias using the delete-one jackkinfe.
        It returns the first and second order bias-corrected estimators.
        
        Parameters
        -------
            No parameters needed. The data and the estimating fuction of theta
            should be declared when create the PanelModel object.

        Returns
        -------
        theta_jack : A numpy array consists the first and second order bias-corrected estimator.

        '''
        T = self.z_F1.shape[2]
        theta_t = np.zeros(T)
        theta_tt_sum = np.zeros(T)

        for t in range(T):
            zt = np.delete(self.z_F1, t, 2)     # delete the t-th period from z_F1
            theta_t[t] = self.estimate_theta(zt)
            theta_tt = np.zeros(T-1)
            for s in range(T-1):
                ztt = np.delete(zt, s, 2)       # delete the s-th period form zt
                theta_tt[s] = self.estimate_theta(ztt)
            theta_tt_sum[t] = theta_tt.sum()
            
        theta_J1 = T*self.theta_hat - (T-1)*theta_t.sum()/T
        theta_J2 = T**2*self.theta_hat/2 - (T-1)**2*theta_t.sum()/T +(T-2)**2*theta_tt_sum.sum()/T/(T-1)/2
        theta_jack = np.r_[theta_J1, theta_J2]
        return theta_jack

    
    def comp_theta_boot(self):
        '''
        This function corrects the IPP asymptotic bias using the delete-one jackkinfe.
        It returns the first and second order bias-corrected estimators.
        
        Parameters
        -------
            No parameters needed. The data and the estimating fuction of ttheta
            should be declared when create the PanelModel object.

        Returns
        -------
        theta_boot : A numpy array consists the first, second, and third order bias-corrected estimator.

        '''
        n, k, T = self.z_F1.shape
        B = C = D = 10      # 10 bootstrap samples is used to estimate the bias

        z_F2 = self.bootstrap_resample(self.z_F1, B)    # generate 10 bootstrap samples from z_F1

        theta_F2 = np.zeros(B)
        theta_F3 = np.zeros((B,C))
        theta_F4 = np.zeros((B,C,D))

        for b in range(B): #tqdm():
            theta_F2[b] = self.estimate_theta(z_F2[:, :, :, b])
            z_F3 = self.bootstrap_resample(z_F2[:, :, :, b], C)     # generate 10 bootstrap samples from z_F2[:,:,:,b] 
            for c in range(C):
                theta_F3[b, c] = self.estimate_theta(z_F3[:, :, :, c])
                z_F4 = self.bootstrap_resample(z_F3[:, :, :, c], D)     # generate 10 bootstrap samples from z_F3[:,:,:,c] 
                for d in range(D):
                    theta_F4[b, c, d] = self.estimate_theta(z_F4[:, :, :, d])

        # On the estimator
        theta_F1 = self.theta_hat
        theta_boot_1 = 2*theta_F1 - theta_F2.mean()
        theta_boot_2 = 3*theta_F1 - 3*theta_F2.mean() + theta_F3.mean()
        theta_boot_3 = 4*theta_F1 - 6*theta_F2.mean() + 4*theta_F3.mean() - theta_F4.mean()
        theta_boot = np.r_[theta_boot_1, theta_boot_2, theta_boot_3]
        return theta_boot
        
    @staticmethod
    def bootstrap_resample(z_F1, B):
        N, K, T = z_F1.shape
        
        z_F2 = np.zeros((N,K,T,B))
        for b in range(B):
            idx = randint.rvs(low=0, high=T, size=[N, T])
            idx = idx + np.arange(n).reshape(-1,1)*T
            idx = idx.reshape(-1,1)[:,0]
            for k in range(K):
                z_F2[:, k, :, b] = z_F1[:, k, :].reshape(-1,1)[idx].reshape(n, T)

        return z_F2

        

if __name__ == "__main__":
    from numpy.random import default_rng
    from scipy.stats import norm
    from scipy.stats import logistic
    from scipy.optimize import root
    
    def main_logit(n, T, rng):
        import warnings
        warnings.filterwarnings("ignore")
        
        beta_0 = 1
        alpha_0 = norm.rvs(size=[n, 1], random_state=rng)*np.sqrt(1) - 0.5
        x = (norm.rvs(size=[n, 1, T], random_state=rng) > 0).astype(int)
        u = logistic.rvs(size=[n, 1, T], random_state=rng)
        y = (alpha_0[:, :, np.newaxis] + beta_0*x + u > 0).astype(int)
        z = np.concatenate((y, x), axis=1)
        
        def estimate_beta(z):
            def estimate_alpha(beta, z):
                y = z[:, [0], :]
                x = z[:, 1:, :]
                z_00 = np.logical_and(y==0, x==0).sum(axis=(1,2))
                z_01 = np.logical_and(y==0, x==1).sum(axis=(1,2))
                z_10 = np.logical_and(y==1, x==0).sum(axis=(1,2))
                z_11 = np.logical_and(y==1, x==1).sum(axis=(1,2))
                
                a = np.exp(beta)*z_10 - np.exp(beta)*z_01 + z_11 - z_00
                alpha_hat = np.log((a + np.sqrt(a**2+4*np.exp(beta)*(z_00+z_01)*(z_10+z_11))) / (2*np.exp(beta)*(z_00+z_01)))
                return alpha_hat
            
            def comp_score_beta(beta, z):
                alpha_hat = estimate_alpha(beta, z)
                y = z[:, [0], :]
                x = z[:, 1:, :]
                z_01 = np.logical_and(y==0, x==1).sum(axis=(1,2))
                z_11 = np.logical_and(y==1, x==1).sum(axis=(1,2))
            
                score = z_11/(1 + np.exp(alpha_hat + beta)) - z_01*np.exp(alpha_hat + beta)/(1 + np.exp(alpha_hat + beta))
                score = np.nan_to_num(score).sum(axis=0)
                return score
        
            beta_guess = 1
            beta_hat = root(comp_score_beta, beta_guess, z).x        
            return beta_hat[0]  
        
        print("Results for Logit with one binary regressor")
        logit = PanelModel(z, estimate_beta)
        beta_jack = logit.comp_theta_jack()
        beta_boot = logit.comp_theta_boot()
        print(f"beta_MLE: {logit.theta_hat}")
        print(f"beta_jack 1st: {beta_jack[0]}, 2nd: {beta_jack[1]}")
        print(f"beta_boot 1st: {beta_boot[0]}, 2nd: {beta_boot[1]}, 3rd: {beta_boot[2]}")

    def main_probit(n, T, rng):
        beta0 = 1
        alpha0 = norm.rvs(size=n, random_state=rng)*np.sqrt(1) - 0.5
        x = (norm.rvs(size=[n, 1, T], random_state=rng) > 0).astype(int)
        u = norm.rvs(size=[n, 1, T], random_state=rng)
        y = (alpha0[:, np.newaxis, np.newaxis] + beta0*x + u > 0).astype(int)
        z = np.concatenate((y, x), axis=1)
        
        def estimate_beta(z):          
            def trans_sample(z, cases):
                y = z[:, [0], :]
                x = z[:, 1:, :]
                z_00 = ((y==0) & (x==0)).sum(axis=(1,2))
                z_01 = ((y==0) & (x==1)).sum(axis=(1,2))
                z_10 = ((y==1) & (x==0)).sum(axis=(1,2))
                z_11 = ((y==1) & (x==1)).sum(axis=(1,2))
                z_ab = np.array([z_00, z_01, z_10, z_11]).T
                
                n_m = np.zeros(cases.shape[0]).astype(int)    
                for m in range(cases.shape[0]):
                    n_m[m] = ((z_ab == cases[m]).sum(axis = 1) == 4).sum() 
                return n_m
            
            def comp_cases(T):
                cases = np.array([[a,b,c,d] for a in range(T+1) for b in range(T+1) for c in range(T+1) for d in range(T+1)])
                cases = cases[(cases.sum(axis=1) == T), :]
                elinimates = (cases.max(axis=1) == T) | (cases[:,0] + cases[:,1] == 0) | (cases[:,2] + cases[:,3] == 0) \
                    | (cases[:,1] + cases[:,3] == 0) | (cases[:,0] + cases[:,2] == 0)
                cases = cases[~elinimates, :]
                return cases
            
            def comp_m(u, threshold=-5):
                m = norm.pdf(u)/norm.cdf(u)
                m[u<threshold] = -u[u<threshold] + 1/(-u[u<threshold] + 2/(-u[u<threshold] + 3/(-u[u<threshold])))
                return m
            
            def comp_d_m(u):
                m = comp_m(u)
                return -u*m - m**2
              
            def comp_score_alpha(alpha, beta, cases):
                # cases, _ = trans_sample(z)
                z_00 = cases[:, 0]
                z_01 = cases[:, 1]
                z_10 = cases[:, 2]
                z_11 = cases[:, 3]
                score_alpha = - z_00*comp_m(-alpha) - z_01*comp_m(-alpha-beta) + z_10*comp_m(alpha) + z_11*comp_m(alpha+beta)
                return score_alpha
            
            def comp_score_beta(beta, alpha_hat, cases, n_m):
                # cases, n_m = trans_sample(z)
                z_01 = cases[:, 1]
                z_11 = cases[:, 3]
                score_beta = np.sum(n_m*(- z_01*comp_m(-alpha_hat-beta) + z_11*comp_m(alpha_hat+beta)))
                return score_beta#np.array([score_beta])
            
            def comp_score_all(parameters, cases, n_m):
                z_00 = cases[:, 0]
                z_01 = cases[:, 1]
                z_10 = cases[:, 2]
                z_11 = cases[:, 3]
                
                beta = parameters[0]
                alpha = parameters[1:]
                
                score_alpha = - z_00*comp_m(-alpha) - z_01*comp_m(-alpha-beta) + z_10*comp_m(alpha) + z_11*comp_m(alpha+beta)
                score_beta = np.sum(n_m*(- z_01*comp_m(-alpha-beta) + z_11*comp_m(alpha+beta)))
                score = np.hstack([score_beta, score_alpha])
                return score
            
            def estimate_theta_cases(cases, n_m):
                alpha_init = np.ones(cases.shape[0])
                beta_init = np.array([1])
                init = np.hstack([beta_init, alpha_init]).astype(float)
                # estimates = run_newton(comp_score_all, comp_delta_all, init, eps=EPS, args=(cases, n_m))
                estimates = root(fun=comp_score_all, x0=init, args=(cases, n_m)).x
                return estimates[0]
            
            
            cases = comp_cases(z.shape[2])
            n_m = trans_sample(z, cases)
            theta_hat = estimate_theta_cases(cases, n_m)
            return theta_hat
        
        print("Results for Probit with one binary regressor")
        probit = PanelModel(z, estimate_beta)
        beta_jack = probit.comp_theta_jack()
        beta_boot = probit.comp_theta_boot()
        print(f"beta_MLE: {probit.theta_hat}")
        print(f"beta_jack 1st: {beta_jack[0]}, 2nd: {beta_jack[1]}")
        print(f"beta_boot 1st: {beta_boot[0]}, 2nd: {beta_boot[1]}, 3rd: {beta_boot[2]}")
        
                       
    rng = default_rng()
    n = 10000
    T = 5
    main_logit(n, T, rng)
    main_probit(n, T, rng)