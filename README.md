# Introduction
Since its initial proposal by Neyman and Scott in 1948, the Incidental Parameter Problem (IPP) has gained recognition and extensive discussion among statisticians and econometricians. Despite various proposed solutions, many practitioners tend to sidestep the IPP by introducing additional assumptions into their models, as seen in the widely used BLP model for demand estimation. There are two primary reasons for this preference. First, numerous solutions, particularly analytical ones, are specific to partic- ular models. For example, Andersen’s (1970) conditional likelihood offers a solution for Logit but not for Probit models. In some cases, no suitable solution may exist for a practitioner’s model. Second, even when appropriate methods are available, practitioners might find it time-consuming to adapt these meth- ods and integrate them into their programming packages. To address this issue, this paper introduces a Python module, numipp.py, which offers two numerical methods for correcting the incidental parameter bias in large-T panel models: the delete-one panel jackknife, as proposed by Hahn and Newey (2004) , and nonparametric bootstrap, a part of the author’s Ph.D. work. These methods are not model-specific and are easy to implement.

# Incidental parameter problem
Suppose we observe a random variable $z_{it}$ (possibly a vector) for cross-sectional units $i=1,\dots,n$ and time periods $t=1,\dots,T$. Let the finite-dimensional parameters $\theta\in\Theta$ be the parameter of interest, $\alpha_{i}$ be a scalar individual effect, i.e., the incidental parameters, and $f\left(z\mid\theta,\alpha\right)$ be the density (or probability mass) function of $z_{it}$ (possibly given strictly exogenous covariates, embedded in $z_{it}$).

Assuming that $z_{it}$ are independent across $i$ and $t$, the fixed effect maximum likelihood estimator (MLE) of $\theta$ is defined as
$$ \widehat{\theta} =\arg\max_{\theta}\sum_{i=1}^{n}\sum_{t=1}^{T}\log f\left(z_{it}\mid\theta,\widehat{\alpha}_{i}\left(\theta\right)\right) $$

$$\begin{aligned}
\widehat{\theta} & =\arg\max_{\theta}\sum_{i=1}^{n}\sum_{t=1}^{T}\log f\left(z_{it}\mid\theta,\widehat{\alpha}_{i}\left(\theta\right)\right)\\
\widehat{\alpha}_{i} & =\arg\max_{\alpha}\sum_{t=1}^{T}\log f\left(z_{it}\mid\theta,\alpha_{i}\right).
\end{aligned}$$ 
That is, we 
\begin{itemize}
\item first use the $T$ observation for $i$ to compute $\widehat{\alpha}_{i}\left(\theta\right)$,
the MLE of $\alpha_{i}$ as a function of $\theta$;
\item plug $\widehat{\alpha}_{i}\left(\theta\right)$ into the likelihood
function $f\left(z_{it}\mid\theta,\alpha_{i}\right)$ and get $f\left(z_{it}\mid\theta,\widehat{\alpha}_{i}\left(\theta\right)\right)$,
known as the profiled likelihood;
\item compute $\widehat{\theta}$ as the maximizer of $\sum_{i=1}^{n}\sum_{t=1}^{T}\log f\left(z_{it}\mid\theta,\widehat{\alpha}_{i}\left(\theta\right)\right)$.
\end{itemize}
%
In short panels where $T$ is small, $\widehat{\theta}_{\text{MLE}}$
is often inconsistent, i.e., as $n\to\infty$ with $T$ fixed, we
have $\widehat{\theta}_{\text{MLE}}\overset{p}{\to}\theta_{T}\neq\theta_{0}$,
where $\theta_{0}$ is the true value of $\theta$. Intuitionally,
when the panel is short, we only have $T$ observations to estimated
each $\alpha_{i}$, leading the estimate of $\alpha_{i}$, $\widehat{\alpha}_{i}\left(\theta\right)$,
inaccurate even if $n$ goes to infinity. Further, this inaccuracy
transfers to the profiled likelihood function and as a consequence,
maximizing it gives an inconsistent $\widehat{\theta}$. For a formal
explanation, see \textcite{hahn2004}. This issue is known as the
incidental parameter problem (IPP), a statistical problem that was
identified by \textcite{neyman1948}. 
