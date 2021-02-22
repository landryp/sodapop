# sodapop
Tools for gravitational-wave population modeling and inference.

### Population inference

1. reweight-masses /home/phil/Research/nspop/GW170817.dat -d ' ' -c m1_source m2_source distance -p flat_mcetadet -r flat_m1m2 -v

2. sample-pop-params flat_m1m2 -p flat,0.7,1.2 flat,2.,3. -v

3. infer-pop-params /home/phil/Research/nspop/flat_m1m2_params.csv /home/phil/Research/nspop/GW170817_reweighted.csv -c m1_source m2_source dL likelihood -p flat_m1m2 -v

4. plot-pop-model /home/phil/Research/nspop/flat_m1m2.csv -p flat_m1m2 -v
