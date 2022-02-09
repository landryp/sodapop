# sodapop
Tools for neutron star population modeling and inference with gravitational-wave observations.

### Scripts

* reweight-masses /path/to/samples.csv -c m1_source m2_source dL z -o /path/to/samples_reweighted.csv -p flat_m1m2det_quad_dL -r flat_m1m2 -v

Reweight posterior samples according to flat-in-source-frame-component-mass prior to make them likelihood samples.

* INFER-NS-MASS

Specify observations, population models and hyperparameter priors, then infer the neutron-star mass distribution. Relies on sample-pop-params, infer-pop-params and build-ppd.
