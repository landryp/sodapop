# sodapop
Tools for gravitational-wave population modeling and inference.

*** this is a legacy branch to retain backwards compatibility for daniellmarc/reweighted-dns ***

### Posterior reweighting

* reweight-masses /home/phil/Research/nspop/GW170817.dat -d ' ' -c m1_source m2_source distance -p flat_mcetadet -r flat_m1m2 -v

Recompute posterior weights with respect to a different mass prior.

* split-m1m2 /home/philippe.landry/nspop/GW170817_reweighted.csv -c m1_source m2_source -v

Make separate posteriors for component masses via marginalization.

### Population inference

* POP-FIT-CONDOR /home/philippe.landry/nspop/dat/LR-bns_peakm1m2+nsbh_unifm1peakm2q2/ /home/philippe.landry/sodapop/ /home/philippe.landry/nspop/etc/GW170817_reweighted.csv /home/philippe.landry/nspop/etc/GW190425_reweighted.csv /home/philippe.landry/nspop/etc/GW200105_reweighted.csv /home/philippe.landry/nspop/etc/GW200115_reweighted.csv -c m1_source m2_source dL likelihood -l 100 -p peakcut_m1m2 peakcut_m1m2 unif_m1_peakcut_m2_qpair unif_m1_peakcut_m2_qpair -P mmin,mmax,mu+flat123,0.5,1.5,1.5,3.,1.,3. sigma+flat,0.01,2. -D quad,0.,1000. -n 108 -B 2. 3. 30. -f snrcut -S 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. -s 100000 -t 20000 -w 9 -b 1000 -o /home/philippe.landry/nspop/dat/LR-bns_peakm1m2+nsbh_unifm1peakm2q2/ -v

Launch condor jobs to infer population parameter posterior from observations by MCMC sampling with emcee relative to a specified population parameter prior.

* sample-pop-params bimodcut_mass -p flat,1.,3. flat,0.01,2. flat,1.,3. flat,0.01,2. flat,0.,1. flat,0.999,1.001 flat,1.5,3. -n 1e4 -v

Generate prior samples in population parameters for a given population model.

* infer-pop-params /home/phil/Research/emcee/test/bimodcut_mass_prior.csv /home/phil/Research/nspop/GW170817_reweighted.csv /home/phil/Research/nspop/GW190425_reweighted.csv /home/phil/Research/nspop/GW200105_reweighted.csv /home/phil/Research/nspop/GW200115_reweighted.csv -c m1_source m2_source likelihood -p bimodcut_m1m2 bimodcut_m1m2 unif_m1_bimodcut_m2 unif_m1_bimodcut_m2 -s snrcut -n 1000 -w 15 -b 100 -m 5000 -N 5000 -o bimodcut_mass.csv -B 3. 30. -S flat_m1m2_quad_dL,0.5,3.,0.5,3.,1.,1000. flat_m1m2_quad_dL,0.5,3.,0.5,3.,1.,1000. flat_m1m2_quad_dL,0.5,30.,0.5,3.,1.,1000. flat_m1m2_quad_dL,0.5,30.,0.5,3.,1.,1000. -P mmin,mmax+flat12,0.999,1.001,1.5,3. mu1,mu2+flat12,1.,3.,1.,3. sigma1+flat,0.01,2.0 sigma2+flat,0.01,2.0 alpha+flat,0.,1. -v

Infer population parameter posterior from observations by MCMC sampling from a given population parameter prior with emcee.

* plot-pop-model /home/phil/Research/nspop/bimod_mass.csv -p bimod_mass -b -n 1e4 -r /home/phil/Research/nspop/GW170817_reweighted-m1.csv,m /home/phil/Research/nspop/GW170817_reweighted-m2.csv,m /home/phil/Research/nspop/GW190425_reweighted-m1.csv,m /home/phil/Research/nspop/GW190425_reweighted-m2.csv,m /home/phil/Research/nspop/GW200105_reweighted-m2.csv,m /home/phil/Research/nspop/GW200115_reweighted-m2.csv,m -f 500 -M -v

* plot-pop-model /home/phil/Research/nspop/bimod_mass.csv -p bimod_mass -b -n 1e4 -r /home/phil/Research/nspop/GW170817_reweighted-m1.csv,m /home/phil/Research/nspop/GW170817_reweighted-m2.csv,m /home/phil/Research/nspop/GW190425_reweighted-m1.csv,m /home/phil/Research/nspop/GW190425_reweighted-m2.csv,m /home/phil/Research/nspop/GW200105_reweighted-m2.csv,m /home/phil/Research/nspop/GW200115_reweighted-m2.csv,m -q 0.05 0.5 0.95 -v

* extract-pop-param-quantiles /home/phil/Research/nspop/bimod_mass.csv mu sigma -w weight -v

* plothist /home/phil/Research/nspop/unif_mass.csv /home/phil/Research/nspop/unif_mass.csv -x mmax -W weight -o unif_mass-mmax.png -l "M_max [M_sun]" "" -L "post" "prior" -u True False -v

### Analyses

POP-FIT-CONDOR /home/philippe.landry/nspop/dat/LR-bns_peakm1m2+nsbh_unifm1peakm2q2/ /home/philippe.landry/sodapop/ /home/philippe.landry/nspop/etc/GW170817_reweighted.csv /home/philippe.landry/nspop/etc/GW190425_reweighted.csv /home/philippe.landry/nspop/etc/GW200105_reweighted.csv /home/philippe.landry/nspop/etc/GW200115_reweighted.csv -c m1_source m2_source dL likelihood -l 100 -p peakcut_m1m2 peakcut_m1m2 unif_m1_peakcut_m2_qpair unif_m1_peakcut_m2_qpair -P mmin,mmax,mu+flat123,0.5,1.5,1.5,3.,1.,3. sigma+flat,0.01,2. -D quad,0.,1000. -n 108 -B 2. 3. 30. -f snrcut -S 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. -s 100000 -t 20000 -w 9 -b 1000 -o /home/philippe.landry/nspop/dat/LR-bns_peakm1m2+nsbh_unifm1peakm2q2/ -v

POP-FIT-CONDOR /home/philippe.landry/nspop/dat/LR-bns_unifm1m2+nsbh_unifm1unifm2q2/ /home/philippe.landry/sodapop/ /home/philippe.landry/nspop/etc/GW170817_reweighted.csv /home/philippe.landry/nspop/etc/GW190425_reweighted.csv /home/philippe.landry/nspop/etc/GW200105_reweighted.csv /home/philippe.landry/nspop/etc/GW200115_reweighted.csv -c m1_source m2_source dL likelihood -l 100 -p unif_m1m2 unif_m1m2 unif_m1_unif_m2_qpair unif_m1_unif_m2_qpair -P mmin,mmax+flat12,0.5,1.5,1.5,3. -D quad,0.,1000. -n 100 -B 2. 3. 30. -f snrcut -S 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. -s 100000 -t 20000 -w 5 -b 1000 -o /home/philippe.landry/nspop/dat/LR-bns_unifm1m2+nsbh_unifm1unifm2q2/ -v

POP-FIT-CONDOR /home/philippe.landry/nspop/dat/LR-bns_bimodm1m2+nsbh_unifm1bimodm2q2/ /home/philippe.landry/sodapop/ /home/philippe.landry/nspop/etc/GW170817_reweighted.csv /home/philippe.landry/nspop/etc/GW190425_reweighted.csv /home/philippe.landry/nspop/etc/GW200105_reweighted.csv /home/philippe.landry/nspop/etc/GW200115_reweighted.csv -c m1_source m2_source dL likelihood -l 100 -p bimodcut_m1m2 bimodcut_m1m2 unif_m1_bimodcut_m2_qpair unif_m1_bimodcut_m2_qpair -P mmin,mmax,mu1,mu2+flat1234,0.5,1.5,1.5,3.,1.,3.,1.,3. sigma1+flat,0.01,2. sigma2+flat,0.01,2. alpha+flat,0.,1. -D quad,0.,1000. -n 105 -B 2. 3. 30. -f snrcut -S 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. 0.5,30.,0.5,3.,0.,1000. -s 100000 -t 20000 -w 15 -b 1000 -o /home/philippe.landry/nspop/dat/LR-bns_bimodm1m2+nsbh_unifm1bimodm2q2/ -v
