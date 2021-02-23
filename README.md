# sodapop
Tools for gravitational-wave population modeling and inference.

### Population inference

1. reweight-masses /home/phil/Research/nspop/GW170817.dat -d ' ' -c m1_source m2_source distance -p flat_mcetadet -r flat_m1m2 -v

2. sample-pop-params doublegaussian_m1m2 -p flat,1.2,1.6 flat,0.01,0.5 flat,1.7,2.1 flat,0.01,0.5 flat,0.,1. -n 1e3 -v

3. infer-pop-params /home/phil/Research/nspop/doublegaussian_m1m2_params.csv /home/phil/Research/nspop/GW170817_reweighted.csv -c m1_source m2_source likelihood -p doublegaussian_m1m2 -v

4. plot-pop-model /home/phil/Research/nspop/doublegaussian_m1m2.csv -p doublegaussian_m1m2 -v
