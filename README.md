# composite-birth-death

Code for likelihood-based inference in finite-state birth-death processes with composite birth mechanisms.

The general problem is a continuous-time Markov population process with an absorbing state and additive birth rates,

$$
\lambda_k(\beta)=\sum_{i=1}^K \beta_i f_i(k),
$$

where each observed upward jump is produced by one of several latent mechanisms. When birth events are unmarked, inference must separate the mechanism-specific rates from a single aggregate trajectory. When observations are conditioned on non-extinction, inference must also correct the survival-selection bias caused by the absorbing state.

The documentation of functions and reproducibility workflow is available at:

[markolalovic.github.io/composite-birth-death](https://markolalovic.github.io/composite-birth-death/)


## Background

We consider likelihood-based inference for continuous-time birth-death processes, or Markov population processes, with composite birth rates. Several distinct mechanisms contribute additively to the total birth intensity. When events are unmarked, the observed path records an upward jump but not which mechanism produced it. Estimating the individual mechanism intensities from a single aggregate trajectory is therefore a deconvolution problem in event space: the state is one-dimensional, but its upward increments arise from a mixture of latent event types.

The second issue is survival conditioning. The state 0 is absorbing, so long observed trajectories are selected by non-extinction. Treating such paths as unconditional samples induces survival-selection bias. The methods implemented here use the Doob $h$-transform and the associated Q-process, an ergodic time-homogeneous surrogate for the law of the original process conditioned on long survival. This leads to survival-conditioned MLEs, QMLEs based on working scores, covariance estimates from Fisher and Godambe information, and one-sided boundary tests for the presence of a given birth mechanism.

The framework here addresses both. The code includes Gillespie simulation, sufficient-statistic extraction, unconditional and conditional likelihood estimators, working-score fixed-point algorithms, and Fisher and Godambe information calculations. The running example is an SIS process on a complete hypergraph, with pairwise and triadic transmission.

Beyond epidemic processes, similar issues of latent event types and survival conditioning can arise in population and extinction biology, chemical reaction networks, and models of social or information diffusion.


## Citation

If you find this useful in your work, please cite the corresponding [paper](https://arxiv.org/abs/2604.20422):

```bibtex
@misc{lalovic2026likelihoodbasedinferencebirthdeathprocesses,
      title={Likelihood-based inference for birth-death processes with composite birth mechanisms}, 
      author={Marko Lalovic and Nicos Georgiou and Istvan Z. Kiss},
      year={2026},
      eprint={2604.20422},
      archivePrefix={arXiv},
      primaryClass={math.ST},
      url={https://arxiv.org/abs/2604.20422}, 
}
```
