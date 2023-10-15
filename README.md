# Introduction
Since its initial proposal by Neyman and Scott in 1948, the Incidental Parameter Problem (IPP) has gained recognition and extensive discussion among statisticians and econometricians. Despite various proposed solutions, many practitioners tend to sidestep the IPP by introducing additional assumptions into their models, as seen in the widely used BLP model for demand estimation. There are two primary reasons for this preference. First, numerous solutions, particularly analytical ones, are specific to partic- ular models. For example, Andersen’s (1970) conditional likelihood offers a solution for Logit but not for Probit models. In some cases, no suitable solution may exist for a practitioner’s model. Second, even when appropriate methods are available, practitioners might find it time-consuming to adapt these meth- ods and integrate them into their programming packages. To address this issue, this paper introduces a Python module, numipp.py, which offers two numerical methods for correcting the incidental parameter bias in large-T panel models: the delete-one panel jackknife, as proposed by Hahn and Newey (2004) , and nonparametric bootstrap, a part of the author’s Ph.D. work. These methods are not model-specific and are easy to implement.

# The incidental parameter problem
Suppose we observe a random variable $z_{it}$ (possibly a vector) for cross-sectional units $i=1,\dots,n$ and time periods $t=1,\dots,T$. Let the finite-dimensional parameters $\theta\in\Theta$ be the parameter of interest, $\alpha_{i}$ be a scalar individual effect, i.e., the incidental parameters, and $f\left(z\mid\theta,\alpha\right)$ be the density (or probability mass) function of $z_{it}$ (possibly given strictly exogenous covariates, embedded in $z_{it}$).

Assuming that $z_{it}$ are independent across $i$ and $t$, the fixed effect maximum likelihood estimator (MLE) of $\theta$ is defined as

$$
\begin{align*}
\widehat{\theta} & =\arg\max _{\theta}\sum _{i=1}^{n}\sum _{t=1}^{T}\log f\left(z _{it}\mid\theta,\widehat{\alpha} _{i}\left(\theta\right)\right)\\
\widehat{\alpha} _{i} & =\arg\max _{\alpha}\sum _{t=1}^{T}\log f\left(z _{it}\mid\theta,\alpha _{i}\right).
\end{align*}
$$

That is, we 

- first use the $T$ observation for $i$ to compute $\widehat{\alpha} _{i}\left(\theta\right)$, the MLE of $\alpha _{i}$ as a function of $\theta$;
- plug $\widehat{\alpha} _{i}\left(\theta\right)$ into the likelihood function $f\left(z _{it}\mid\theta,\alpha _{i}\right)$ and get $f\left(z _{it}\mid\theta,\widehat{\alpha} _{i}\left(\theta\right)\right)$, known as the profiled likelihood;
- compute $\widehat{\theta}$ as the maximizer of $\sum _{i}  \sum _{t} \log f\left(z _{it}\mid\theta,\widehat{\alpha} _{i}\left(\theta\right)\right)$.

In short panels where $T$ is small, $\widehat{\theta} _{\text{MLE}}$ is often inconsistent, i.e., as $n\to\infty$ with $T$ fixed, we have $\widehat{\theta} _{\text{MLE}}\overset{p}{\to}\theta _{T}\neq\theta _{0}$, where $\theta _{0}$ is the true value of $\theta$. Intuitionally, when the panel is short, we only have $T$ observations to estimate each $\alpha _{i}$, leading the estimate of $\alpha _{i}$, $\widehat{\alpha} _{i}\left(\theta\right)$, inaccurate even if $n$ goes to infinity. Further, this inaccuracy transfers to the profiled likelihood function, and as a consequence, maximizing it gives an inconsistent $\widehat{\theta}$. For a formal explanation, see Hahn and Newey (2004). This issue is known as the incidental parameter problem (IPP), a statistical problem that was identified by Neyman and Scott (1948). 

# The delete-one panel jackknife

Hahn and Newey (2004) define the ``delete-one'' jackknife estimator defined as  

$$
\begin{align*}
\widehat{\theta} _{1}^{\mathrm{jack}} & =T\widehat{\theta}-(T-1)\frac{\sum _{t}\widehat{\theta} _{(t)}}{T},
\end{align*}
$$

where $\widehat{\theta} _{(t)}$ is $\widehat{\theta}$ computed without the observations from $t$-th time period. 

## Higher order correction

This method can also be applied to yield higher-order bias corrections, although not literally by iterating the delete-one jackknife. Extending the delete-one panel jackknife to the second order gives the delete-two panel jackknife estimator 

$$
\begin{align*}
\widehat{\theta} _{2}^{\mathrm{jack}} & =\frac{1}{2}T^{2}\widehat{\theta}-(T-1)^{2}\frac{\sum _{t=1}^{T}\widehat{\theta} _{(t)}}{T}+\frac{1}{2}(T-2)^{2}\frac{\sum _{t=1}^{T}\sum _{s=t+1}^{T}\widehat{\theta} _{(t,s)}}{T(T-1)/2},
\end{align*}
$$

where $\widehat{\theta} _{(t,s)}$ is $\widehat{\theta}$ computed without the observations from the $t$-th and $s$-th time periods. 

## Discussion

It is easy to see that $\widehat{\theta} _{1}^{\mathrm{jack}}$ reduces the incidental parameter bias of $\widehat{\theta}$. Recall that, as $n\to\infty$ with $T$ fixed, we have $\widehat{\theta}\overset{p}{\to}\theta _{T}$. Here, under suitable regularity conditions, $\theta _{T}$ can be expanded as

$$
\begin{align}
\theta _{T} & =\theta _{0}+\frac{B _{1}}{T}+\frac{B _{2}}{T^{2}}+\ldots+\frac{B _{k}}{T^{k}}+o(T^{-k})\qquad\textrm{as }T\to\infty, \qquad(1)
\end{align}
$$

where $B _{1},\ldots,B _{k}$ are constants (i.e., not depending on $T$) and $k$ is a positive integer. As $n\to\infty$, we have $\widehat{\theta} _{(t)}\overset{p}{\to}\theta _{T-1}$ and, therefore, $\widehat{\theta} _{1}^{\mathrm{jack}}\overset{p}{\to}\theta _{T}^{\mathrm{jack}}$ where 

$$
\begin{align*}
\theta _{T}^{\mathrm{jack}} & =T\theta _{T}-(T-1)\theta _{T-1}\\
 & =\theta _{0}+T\left(\frac{B _{1}}{T}+\frac{B _{2}}{T^{2}}+\ldots+\frac{B _{k}}{T^{k}}+o(T^{-k})\right)\\
 & \qquad-(T-1)\left(\frac{B _{1}}{T-1}+\frac{B _{2}}{(T-1)^{2}}+\ldots+\frac{B _{k}}{(T-1)^{k}}+o((T-1)^{-k})\right)\\
 & =\theta _{0}+\frac{B _{2}}{T^{2}}+o(T^{-2}),
\end{align*}
$$

provided that $k\ge2$. Hence the delete-one jackknife eliminates the bias term $B _{1}/T$.

For $\widehat{\theta} _{2}^{\mathrm{jack}},$ it can be shown that, as $n\to\infty$ with $T$ fixed, $\widehat{\theta} _{2}^{\mathrm{jack}}\overset{\mathrm{p}}{\to}\theta _{0}+O(T^{-3})$, that is, the second-order jackknife eliminates the terms $B _{1}/T$ and $B _{2}/T^{2}$ from Eq(1). 

# Nonparametric bootstrap

Bootstrap can be used to correct the IPP bias as well. The application of parametric bootstrap on IPP is discussed in  Pace and Salvan (2006) and Kim and Sun (2016). Here, I introduce using nonparametric bootstrap. 

Let $F _{0}$ denote the true distribution of $z$, given $\alpha _{0}=(\alpha _{1,0},\ldots,\alpha _{n,0})$, i.e., $z\sim F _{0}$ and $\mathbb{E} _{F _{0}}$ denote the expectations under $F _{0}$. Let $\widehat{F}(z)$ be the empirical distribution defined as follows by $z$. For $z^{\*}\sim\widehat{F}(z)$, $z^{\*}$ has the same dimension, $T\times n$, as $z$, and $z _{it}^{\*}$ is independent across $i$ and independent and identically distributed
across $t$ with

$$
\mathrm{Pr} _{\widehat{F}(z)}[z _{it}^{\*}=z _{is}]=\frac{1}{T}\qquad(i=1,\ldots,n;t,s=1,\ldots,T).
$$

That is, sampling from $\widehat{F}(z)$ amounts, for each $i$, to resampling with replacement from the elements of $z _{i}$ to form $z _{i}^{\*}:=(z _{i1}^{\*},\ldots,z _{iT}^{\*})'$, and then $z^{\*}:=(z _{1}^{\*},\ldots,z _{n}^{\*})$. Note that this notation comprises nested bootstrap sampling: if $z^{\*}$ is a bootstrap sample (i.e., drawn from $\widehat{F}(z)$), then drawing from $\widehat{F}(z^{\*})$ gives a nested bootstrap sample, say $z^{\*\*}$, and so on.

The bias of $\widehat{\theta}$ is 

$$
\begin{align*}
\mathrm{bias}(\widehat{\theta}) & =\mathbb{E} _{F _{0}}[\widehat{\theta}(z)-\theta _{0}],
\end{align*}
$$

where we write $\widehat{\theta}$ as $\widehat{\theta}(z)$ to indicate that the estimator is based on the data $z$. By the bootstrap principle, on replacing $F _{0}$ with $\widehat{F}(z)$ and $z$ with $z^{\*}\sim\widehat{F}(z)$, we obtain the ideal first-order bootstrap estimate of the bias as 

$$
\begin{align*}
\widehat{\mathrm{bias}}(\widehat{\theta}) & =\mathbb{E} _{\widehat{F}(z)}[\widehat{\theta}(z^{\*})-\widehat{\theta}(z)]\\
 & =\mathbb{E} _{\widehat{F}(z)}[\widehat{\theta}(z^{\*})]-\widehat{\theta}.
\end{align*}
$$

Hence 

$$
\begin{align*}
\widehat{\theta} _{1}^{\mathrm{boot}} & =\widehat{\theta}-\widehat{\mathrm{bias}}(\widehat{\theta})\\
 & =2\widehat{\theta}-\mathbb{E} _{\widehat{F}(z)}[\widehat{\theta}(z^{\*})]
\end{align*}
$$

is an ideal first-order bootstrap bias-corrected estimator of $\theta _{0}$.

## Higher order correction

The bias of $\widehat{\theta} _{1}^{\mathrm{boot}}$ is 

$$
\begin{align*}
\mathrm{bias}(\widehat{\theta} _{1}^{\mathrm{boot}}) & =\mathbb{E} _{F _{0}}[\widehat{\theta} _{1}^{\mathrm{boot}}(z)-\theta _{0}]\\
 & =\mathbb{E} _{F _{0}}\big[2\widehat{\theta}(z)-\mathbb{E} _{\widehat{F}(z)}[\widehat{\theta}(z^{\*})]-\theta _{0}\big].
\end{align*}
$$

Applying the bootstrap principle again, we replace $F _{0}$ with $\widehat{F}\left(z\right)$, $\widehat{F}\left(z\right)$ with $\widehat{F}\left(z^{\*}\right)$, $z$ with $z^{\*}\sim\widehat{F}\left(z\right)$, and $z^{\*}$ with $z^{\*\*}\sim\widehat{F}\left(z^{\*}\right)$ to obtain the ideal second-order bootstrap estimate of the bias as 

$$
\begin{align*}
\widehat{\mathrm{bias}}(\widehat{\theta} _{1}^{\mathrm{boot}}) & =\mathbb{E} _{\widehat{F}(z)}\big[2\widehat{\theta}(z^{\*})-\mathbb{E} _{\widehat{F}(z^{\*})}[\widehat{\theta}(z^{\*\*})]-\widehat{\theta}(z)\big]\\
 & =\mathbb{E} _{\widehat{F}(z)}\big[2\widehat{\theta}(z^{\*})-\mathbb{E} _{\widehat{F}(z^{\*})}[\widehat{\theta}(z^{\*\*})]\big]-\widehat{\theta}.
\end{align*}
$$

Hence 

$$
\begin{align*}
\widehat{\theta} _{2}^{\mathrm{boot}} & =\widehat{\theta} _{1}^{\mathrm{boot}}-\widehat{\mathrm{bias}}(\widehat{\theta} _{1}^{\mathrm{boot}})\\
 & =3\widehat{\theta}-3\mathbb{E} _{\widehat{F}(z)}[\widehat{\theta}(z^{\*})]+\mathbb{E} _{\widehat{F}(z)}\mathbb{E} _{\widehat{F}(z^{\*})}[\widehat{\theta}(z^{\*\*})]
\end{align*}
$$

is an ideal second-order bootstrap bias-corrected estimator of $\theta _{0}$.

For higher order corrections, let $z^{\*(0)}=z$ and, for $b=0,1,2,\dots$, define $z^{\*(b+1)}\sim\widehat{F}(z^{\*(b)})$ recursively, so that $z^{\*(1)}=z^{\*}$, $z^{\*(2)}=z^{\*\*}$, and so on. It is easy to show that the $K$-th order ideal bootstrap bias corrections are  

$$
\begin{align*}
\widehat{\theta} _{K}^{\mathrm{boot}} & =(K+1)\widehat{\theta}-\sum _{k=2}^{K+1}(-1)^{k}\binom{K+1}{k}\mathbb{E} _{\widehat{F}(z)}\mathbb{E} _{\widehat{F}(z^{\*})}\cdots\mathbb{E} _{\widehat{F}(z^{\*(k-1)})}[\widehat{\theta}(z^{\*(k-1)})].
\end{align*}
$$

## Bootstrap simulations

Unlike the jackknife, for which there are closed-form expressions (see also the expression below for the second-order jackknife), the bootstrap bias corrections are usually not available in closed form. Except in very simple cases, the expectations have to be computed by simulation. This involves the following steps to draw nested bootstrap samples, all of which are of the same dimension, $T\times n$, as the original data set $z$. 

- Draw $B$ bootstrap samples from $\widehat{F}(z)$ as described in the beginning of this section and denote these as $z _{b _{1}}^{\*}$, $b _{1}=1,\dots,B$.
- For each $z _{b _{1}}^{\*}$, draw $B$ bootstrap samples from $\widehat{F}(z _{b _{1}}^{\*})$ and denote these as $z _{b _{1}b _{2}}^{\*(2)}$, $b _{2}=1,\dots,B$.
- Iterate this procedure to form further bootstrap samples $z _{b _{1}b _{2}b _{3}}^{\*(3)},\ldots,z _{b _{1}b _{2}b _{3}\cdots b _{K}}^{\*(K)}$, where $b _{j}=1,\ldots,B$ for each $j=1,\ldots,K$.


Now we can approximate $\widehat{\theta} _{K}^{\mathrm{boot}}$ using 

$$
\begin{align*}
\widehat{\theta} _{K}^{\mathrm{boot}} & \approx(K+1)\widehat{\theta}(z)-\sum _{k=2}^{K+1}\binom{K+1}{k}(-1)^{k}\frac{1}{B^{k-1}}\sum _{b _{1}=1}^{B}\sum _{b _{2}=1}^{B}\cdots\sum _{b _{k-1}=1}^{B}\widehat{\theta}(z _{b _{1}b _{2}\cdots b _{k-1}}^{\*(k-1)}).
\end{align*}
$$

##  Discussion

For the first order correction, because the bootstrap estimates the bias with relative error $O(T^{-1})$, the bootstrap can be thought of as eliminating $B _{1}/T$. At the same time, the bootstrap also modifies the higher-order terms in the expansion, but without modifying their order of magnitude. Thus, the bootstrap reduces the bias from $O(T^{-1})$ (the bias of $\widehat{\theta}$) down to $O(T^{-2})$ (the bias of $\widehat{\theta} _{1}^{\mathrm{boot}}$). Similar reasoning applies to higher-order bootstrap correction.

# Examples

Here we discuss some models where an incidental parameter problem occurs, and what the jackknife and nonparametric bootstrap deliver in terms of bias reduction.

## Logit model

The fixed-effect panel logit model is a well-known example of a parametric model that suffers from an incidental parameter problem. Let $z _{it}=(y _{it},x _{it})$, where $y _{it}$ is a binary outcome variable and $x _{it}$ is a vector
of exogenous regressors. The model specifies  

$$
\begin{align*}
\Pr[y _{it}=1|x _{it}]=\Lambda(\alpha _{i0}+x' _{it}\beta _{0})
\end{align*}
$$

where $\Lambda(w)=(1-e^{-w})^{-1}$ is the standard logistic distribution function. Hence the probability mass function is 

$$
\begin{align*}
f(y _{it}|x _{it};\beta,\alpha _{i})=[\Lambda(\alpha _{i}+x' _{it}\beta)]^{y _{it}}[1-\Lambda(\alpha _{i}+x' _{it}\beta)]^{1-y _{it}},\qquad y _{it}\in\{0,1\},
\end{align*}
$$

where $\theta=\beta$ is the parameter of interest and the log-likelihood function, normalized by the number of observations, is 

$$
\begin{align*}
l(\beta,\alpha;z) & =\frac{1}{nT}\sum _{i=1}^{n}\sum _{t=1}^{T}\log f(y _{it}|x _{it};\beta,\alpha _{i})\\
 & =\frac{1}{nT}\sum _{i=1}^{n}\sum _{t=1}^{T}\big[y _{it}\log\Lambda(\alpha _{i}+x _{it}\beta)+(1-y _{it})\log(1-\Lambda(\alpha _{i}+x _{it}\beta))\big].
\end{align*}
$$

The maximum likelihood estimator of $\beta _{0}$ is inconsistent as $n\to\infty$ with $T$ fixed.(In the logit model, there is a sufficient statistic for $\alpha _{i}$ and conditioning on it yields a conditional likelihood that resolves the incidental parameter problem. However, the logit model is the only binary-choice model where such a solution exists. Our interest here is in using the logit model as a test case for bootstrap corrections, which have a much wider scope of application than conditional likelihood.

The profile score function is 

$$
\begin{align*}
s(\beta,\widehat{\alpha}(\beta;z);z) & =\frac{1}{nT}\sum _{i=1}^{n}\sum _{t=1}^{T}\big(y _{it}-\Lambda(\widehat{\alpha} _{i}(\beta;z _{i})+x' _{it}\beta)\big)x _{it},
\end{align*}
$$

where $\widehat{\alpha} _{i}(\beta;z _{i})$ solves 

$$
\begin{align*}
\frac{1}{T}\sum _{t=1}^{T}\big(y _{it}-\Lambda(\alpha _{i}+x' _{it}\beta)\big)=0
\end{align*}
$$

for $\alpha _{i}$. Unfortunately, the solution $\widehat{\alpha} _{i}(\beta;z _{i})$ cannot be obtained in closed form from the latter equation. We shall, therefore, focus on the special case where $x _{it}\in\{0,1\}$ (i.e., $x _{it}$ is a binary scalar) because then $\widehat{\alpha} _{i}(\beta;z _{i})$ can be obtained in closed form. With binary $x _{it}$, the log-likelihood function simplifies to 

$$
\begin{align*}
l(\beta,\alpha;z) & =\frac{1}{nT}\sum _{i=1}^{n}(z _{i,10}\log\Lambda(\alpha _{i})+z _{i,00}\log\Lambda(-\alpha _{i})+z _{i,11}\log\Lambda(\alpha _{i}+\beta)\\
 & \qquad\qquad\quad+z _{i,01}\log\Lambda(-\alpha _{i}-\beta))+c,
\end{align*}
$$

where $c$ is an inessential constant and 

$$
\begin{align*}
z _{i,yx}=\sum _{t=1}^{T} 1(Y _{it}=y,X _{it}=x)\qquad{\textrm{ for }}y,x\in\{0,1\}.
\end{align*}
$$

Now $\widehat{\alpha} _{i}(\beta;z _{i})$ is given by the following lemma.

**Lemmma 1.** Let $x _{it}\in\{0,1\}$ for all $i$ and $t$. Then 

$$
\begin{align*}
\widehat{\alpha} _{i}(\beta;z _{i})=\log\frac{-b+\sqrt{b^{2}-4ac}}{2a}
\end{align*}
$$

where 

$$
\begin{align*}
a & =e^{\beta}(z _{i,00}+z _{i,01}),\\
b & =e^{\beta}z _{i,01}-e^{\beta}z _{i,10}+z _{i,00}-z _{i,11},\\
c & =-z _{i,10}-z _{i,11}.
\end{align*}
$$

The proof of the lemma is given in the appendix.

Now the profile score function simplifies to 

$$
\begin{align*}
s(\beta,\widehat{\alpha}(\beta;z);z) & =\frac{1}{nT}\sum _{i=1}^{n}\left(\frac{z _{i,11}-z _{i,01}e^{\widehat{\alpha} _{i}(\beta;z _{i})+\beta}}{1+e^{\widehat{\alpha} _{i}(\beta;z _{i})+\beta}}\right)
\end{align*}
$$

and it is easy to solve $s(\beta,\widehat{\alpha}(\beta;z);z)=0$ numerically for $\beta$, and similarly for solving the bias-corrected profile score equations.

We evaluated the performance of the bootstrap numerically in a small-scale setup with $n=10,000$ (so as to get close to the probability limits as $n\to\infty$), $T\in\{3,4,5,10\}$, and $B=10$. We generated the binary variables $x _{it}$ and $y _{it}$ according to 

$$
\begin{align*}
 & \Pr[x _{it}=1]=1/2,\\
 & \Pr[y _{it}=1|x _{it}]=\Lambda(\alpha _{i0}+x _{it}\beta _{0}),
\end{align*}
$$

with $\beta _{0}=1$ and $\alpha _{i0}\sim\mathcal{N}\left(-0.5,1\right)$. Table 1 reports the means and the standard deviations of the estimators, estimated from $1,000$ Monte Carlo replications. Clearly, the bootstrap performs reasonably well as a bias correction method, and is on par with the jackknife. The bias corrections have the effect of shrinking the maximum likelihood estimate and, as a result, also tend to reduce the standard deviation of the estimator, most prominently for the first-order bias corrections and when $T$ is very small.

<div align="center">
 
| $\beta _{0}=1$ | | $T=3$ | $T=4$ | $T=5$ | $T=10$ |
| --- | --- | --- | --- | --- | --- |
| $\widehat{\beta}$                    | mean<br>std | $1.5272$<br>$0.0534$ | $1.3509$<br>$0.0371$ | $1.2608$<br>$0.0294$ | $1.1152$<br>$0.0176$ |
| $\widehat{\beta} _{1}^{\text{jack}}$ | mean<br>std | $0.5810$<br>$0.0322$ | $0.8178$<br>$0.0219$ | $0.9049$<br>$0.0203$ | $0.9845$<br>$0.0154$ |
| $\widehat{\beta} _{2}^{\text{jack}}$ | mean<br>std | NA<br>NA             | $1.0534$<br>$0.0245$ | $1.0362$<br>$0.0236$ | $1.0027$<br>$0.0157$ |
| $\widehat{\beta} _{1}^{\text{boot}}$ | mean<br>std | $1.0998$<br>$0.0440$ | $0.9916$<br>$0.0307$ | $0.9758$<br>$0.0254$ | $0.9890$<br>$0.0165$ |
| $\widehat{\beta} _{2}^{\text{boot}}$ | mean<br>std | $0.8706$<br>$0.0504$ | $0.8925$<br>$0.0357$ | $0.9396$<br>$0.0304$ | $0.9935$<br>$0.0191$ |
| $\widehat{\beta} _{3}^{\text{boot}}$ | mean<br>std | $0.7640$<br>$0.0671$ | $0.8957$<br>$0.0468$ | $0.9630$<br>$0.0391$ | $0.9988$<br>$0.0232$ |

Table 1: Means and standard deviations of maximum likelihood and bias-corrected estimators in the panel logit model with fixed effects and a binary regressor.
</div>

## Probit model

The fixed-effect panel probit model is very similar to the logit model, but the incidental parameter problem is far more challenging. The probit model specifies 

$$
\begin{align*}
\Pr[y _{it}=1|x _{it}]=\Phi(\alpha _{i0}+x' _{it}\beta _{0})
\end{align*}
$$

where $\Phi(\cdot)$ is the standard normal distribution function. The probability mass function and the log-likelihood function, $f(y _{it}|x _{it};\beta,\alpha _{i})$ and $l(\beta,\alpha;z)$, are as in logit model, but with $\Phi(\cdot)$ instead of $\Lambda(\cdot)$. As in the logit model, the maximum likelihood estimator of $\beta _{0}$ is inconsistent as $n\to\infty$ with $T$ fixed. Unlike the logit model, there is no conditional likelihood that is free of incidental parameters. Moreover, \textcite{chamberlain2010} showed that, when $T=2$ and $x _{it}$ is binary, $\beta _{0}$ is not point-identified. One may conjecture that this holds for any fixed $T$ and for any $x _{it}$ with bounded support. Therefore, a complete resolution of the incidental parameter problem would require calculating (and estimating) the identified set for $\beta _{0}$. Here, instead, we explore the performance of the bootstrap as an approximate (and simpler) approach toward the incidental parameter problem.

The profile score function is 

$$
\begin{align*}
s(\beta,\widehat{\alpha}(\beta;z);z) & =\frac{1}{nT}\sum _{i=1}^{n}\sum _{t=1}^{T}\big(y _{it}-\Phi(\widehat{\alpha} _{i}(\beta;z _{i})+x' _{it}\beta)\big)g(\widehat{\alpha} _{i}(\beta;z _{i})+x' _{it}\beta)x _{it},
\end{align*}
$$

where $\widehat{\alpha} _{i}(\beta;z _{i})$ solves 

$$
\begin{align*}
\frac{1}{T}\sum _{t=1}^{T}\big(y _{it}-\Phi(\alpha _{i}+x' _{it}\beta)\big)g(\alpha _{i}+x' _{it}\beta)=0
\end{align*}
$$

for $\alpha _{i}$. Again, $\widehat{\alpha} _{i}(\beta;z _{i})$ cannot be obtained in closed form. As in the logit model, we consider the case where $x _{it}$ is binary to simplify the analysis. Then $\widehat{\alpha} _{i}(\beta;z _{i})$ solves the simpler equation 

$$
\begin{align*}
\frac{1}{T}\left(z _{i,00}\frac{-\phi(\alpha _{i})}{1-\Phi(\alpha _{i})}+z _{i,01}\frac{-\phi(\alpha _{i}+\beta)}{1-\Phi(\alpha _{i}+\beta)}+z _{i,10}\frac{\phi(\alpha _{i})}{1-\Phi(\alpha _{i})}+z _{i,11}\frac{\phi(\alpha _{i}+\beta)}{1-\Phi(\alpha _{i}+\beta)}\right)=0
\end{align*}
$$

for $\alpha _{i}$, where $z _{i,yx}$ (for $y,x\in\{0,1\}$) is defined as in the logit model. We solve this equation numerically to obtain $\widehat{\alpha} _{i}(\beta;z _{i})$ and subsequently obtain $\widehat{\beta}$ from solving $s(\beta,\widehat{\alpha}(\beta;z);z)=0$ for $\beta$, and similarly for the bias-corrected estimators. Thus, our numerical procedure to obtain any estimate of $\beta$ is a nested one, where the inner loop computes $\widehat{\alpha} _{i}(\beta;z _{i})$.

Our setup to evaluate the performance of the bootstrap is similar to that for the logit model: $n=10,000$, $T\in\{3,4,5,10\}$, $B=10$, (The probit model is computationally much more demanding than the logit model, so we set $B=3$ instead of $B=10$. The results, however, do not depend much on $B$. In fact, $B=1$ already suffices for bias reduction, although there is some variance inflation when $T$ is small.) and 

$$
\begin{align*}
\Pr[x _{it}=1]=1/2,\qquad\Pr[y _{it}=1|x _{it}]=\Phi(\alpha _{i0}+x _{it}\beta _{0}),
\end{align*}
$$

with $\beta _{0}=1$ and $\alpha _{i0}\sim\mathcal{N}\left(-0.5,1\right)$. Table 2 reports the means and standard deviations of the estimators, estimated from $1,000$ Monte Carlo replications. Again, the bootstrap performs well, broadly similar to the jackknife.

<div align="center">
 
| $\beta _{0}=1$ | | $T=3$ | $T=4$ | $T=5$ | $T=10$ |
| --- | --- | --- | --- | --- | --- |
| $\widehat{\beta}$                    | mean<br>std | $1.5962$<br>$0.0393$ | $1.3964$<br>$0.0272$ | $1.2946$<br>$0.0204$ | $1.1253$<br>$0.0114$ |
| $\widehat{\beta} _{1}^{\text{jack}}$ | mean<br>std | $0.6873$<br>$0.0311$ | $0.8014$<br>$0.0161$ | $0.8844$<br>$0.0133$ | $0.9767$<br>$0.0097$ | 
| $\widehat{\beta} _{2}^{\text{jack}}$ | mean<br>std | NA<br>NA             | $0.9160$<br>$0.0128$ | $1.0084$<br>$0.0149$ | $1.0014$<br>$0.0101$ | 
| $\widehat{\beta} _{1}^{\text{boot}}$ | mean<br>std | $1.2971$<br>$0.0374$ | $1.1030$<br>$0.0239$ | $1.0446$<br>$0.0173$ | $1.0032$<br>$0.0107$ | 
| $\widehat{\beta} _{2}^{\text{boot}}$ | mean<br>std | $1.1155$<br>$0.0417$ | $0.9923$<br>$0.0267$ | $0.9856$<br>$0.0196$ | $0.9994$<br>$0.0125$ | 
| $\widehat{\beta} _{3}^{\text{boot}}$ | mean<br>std | $1.0096$<br>$0.0510$ | $0.9633$<br>$0.0340$ | $0.9865$<br>$0.0248$ | $1.0031$<br>$0.0153$ | 

Table 2: Means and standard deviations of maximum likelihood and bias-corrected estimators in the panel probite model with fixed effects and a binary regressor.
</div>

# Appendix: Proof of Lemma 1

The score function for $\alpha _{i}$ is 

$$
\begin{align*}
\partial _{\alpha _{i}}l(\beta,\alpha;z) & =\frac{z _{i,10}}{1+e^{\alpha _{i}}}-\frac{z _{i,00}e^{\alpha _{i}}}{1+e^{\alpha _{i}}}+\frac{z _{i11}}{1+e^{\alpha _{i}+\beta}}-\frac{z _{i01}e^{\alpha _{i}+\beta}}{1+e^{\alpha _{i}+\beta}}.
\end{align*}
$$

Solving $\partial _{\alpha _{i}}l(\beta,\alpha;z)=0$ for $\alpha _{i}$ gives $\widehat{\alpha} _{i}(\beta;z _{i})$. We can rearrange $\partial _{\alpha _{i}}l(\beta,\alpha;z)=0$ as

$$
\begin{align*}
ae^{2\alpha _{i}}+be^{\alpha _{i}}+c & =0.
\end{align*}
$$

This is a quadratic equation in $e^{\alpha _{i}}$ with unique solution

$$
\begin{align*}
e^{\alpha _{i}}=\frac{-b+\sqrt{b^{2}-4ac}}{2a}
\end{align*}
$$

since $e^{\alpha _{i}}>0$ rules out the solution with the negative root. Taking the logarithm gives $\widehat{\alpha} _{i}(\beta;z _{i})$.

# Reference
Andersen, Erling Bernhard (1970). “Asymptotic Properties of Conditional Maximum-Likelihood Estimators”. In: Journal of the Royal Statistical Society. Series B (Methodological) 32.2, pp. 283–301.

Chamberlain, Gary (2010). “Binary Response Models for Panel Data: Identification and Information”. In: Econometrica 78.1, pp. 159–168.

Fernández-Val, Iván and Martin Weidner (2018). “Fixed Effects Estimation of Large-T Panel Data Models”. In: Annual Review of Economics 10.1, pp. 109–138.

Hahn, Jinyong and Whitney Newey (2004). “Jackknife and Analytical Bias Reduction for Nonlinear Panel Models”. In: Econometrica 72.4, pp. 1295–1319.

Kim, Min Seong and Yixiao Sun (2016). “Bootstrap and K-Step Bootstrap Bias Correction for Fixed Effects Estimators in Nonlinear Panel Data Models”. In: Econometric Theory 32.6, pp. 1523–1568.

Neyman, J. and Elizabeth L. Scott (1948). “Consistent Estimates Based on Partially Consistent Observations”. In: Econometrica 16.1, pp. 1–32.

Pace, Luigi and Alessandra Salvan (2006). “Adjustments of the Profile Likelihood from a New Perspective”. In: Journal of Statistical Planning and Inference 136.10, pp. 3554–3564.
