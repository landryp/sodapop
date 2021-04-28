# sodapop
Tools for gravitational-wave population modeling and inference.

### Posterior reweighting

* reweight-masses /home/phil/Research/nspop/GW170817.dat -d ' ' -c m1_source m2_source distance -p flat_mcetadet -r flat_m1m2 -v

Recompute posterior weights with respect to a different mass prior.

* split-m1m2 /home/philippe.landry/nspop/GW170817_reweighted.csv -c m1_source m2_source -v

Make separate posteriors for component masses via marginalization.

### Population inference

* sample-pop-params bimodcut_mass -p flat,1.,3. flat,0.01,2. flat,1.,3. flat,0.01,2. flat,0.,1. flat,0.999,1.001 flat,1.5,3. -n 1e4 -v

Generate prior samples in population parameters for a given population model.

* infer-pop-params /home/phil/Research/nspop/bimodcut_mass_prior.csv /home/phil/Research/nspop/GW170817_reweighted.csv /home/phil/Research/nspop/GW190425_reweighted.csv /home/phil/Research/nspop/GW200105_reweighted.csv /home/phil/Research/nspop/GW200115_reweighted.csv -c m1_source m2_source likelihood -p bimodcut_m1m2 bimodcut_m1m2 unif_m1_bimodcut_m2 unif_m1_bimodcut_m2 -s snrcut -P 10000 -w 50 -b 1000 -m 1000 -S 1000 -o bimodcut_mass.csv -v

Infer population parameter posterior from observations by MCMC sampling from a given population parameter prior with emcee.

* plot-pop-model /home/phil/Research/nspop/bimod_mass.csv -p bimod_mass -b -n 1e4 -r /home/phil/Research/nspop/GW170817_reweighted-m1.csv,m /home/phil/Research/nspop/GW170817_reweighted-m2.csv,m /home/phil/Research/nspop/GW190425_reweighted-m1.csv,m /home/phil/Research/nspop/GW190425_reweighted-m2.csv,m /home/phil/Research/nspop/GW200105_reweighted-m2.csv,m /home/phil/Research/nspop/GW200115_reweighted-m2.csv,m -f 500 -M -v

* plot-pop-model /home/phil/Research/nspop/bimod_mass.csv -p bimod_mass -b -n 1e4 -r /home/phil/Research/nspop/GW170817_reweighted-m1.csv,m /home/phil/Research/nspop/GW170817_reweighted-m2.csv,m /home/phil/Research/nspop/GW190425_reweighted-m1.csv,m /home/phil/Research/nspop/GW190425_reweighted-m2.csv,m /home/phil/Research/nspop/GW200105_reweighted-m2.csv,m /home/phil/Research/nspop/GW200115_reweighted-m2.csv,m -q 0.05 0.5 0.95 -v

* extract-pop-param-quantiles /home/phil/Research/nspop/bimod_mass.csv mu sigma -w weight -v

* plothist /home/phil/Research/nspop/unif_mass.csv /home/phil/Research/nspop/unif_mass.csv -x mmax -W weight -o unif_mass-mmax.png -l "M_max [M_sun]" "" -L "post" "prior" -u True False -v
