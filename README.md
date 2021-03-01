# sodapop
Tools for gravitational-wave population modeling and inference.

### Population inference

1. reweight-masses /home/phil/Research/nspop/GW170817.dat -d ' ' -c m1_source m2_source distance -p flat_mcetadet -r flat_m1m2 -v

2. split-m1m2 /home/philippe.landry/nspop/GW170817_reweighted.csv -c m1_source m2_source -v

3. sample-pop-params bimod_mass -p flat,0.7,2.5 flat,0.01,2. flat,0.7,2.5 flat,0.01,2. flat,0.,1. -n 1e4 -v

4. infer-pop-params /home/phil/Research/nspop/bimod_mass_params.csv /home/phil/Research/nspop/GW170817_reweighted-m1.csv /home/phil/Research/nspop/GW170817_reweighted-m2.csv /home/phil/Research/nspop/GW190425_reweighted-m1.csv /home/phil/Research/nspop/GW190425_reweighted-m2.csv /home/phil/Research/nspop/GW200105_reweighted-m2.csv /home/phil/Research/nspop/GW200115_reweighted-m2.csv -c m1_source m2_source likelihood -p bimod_mass -v

5. plot-pop-model /home/phil/Research/nspop/bimod_mass.csv -p bimod_mass -b -n 1e4 -r /home/phil/Research/nspop/GW170817_reweighted-m1.csv,m /home/phil/Research/nspop/GW170817_reweighted-m2.csv,m /home/phil/Research/nspop/GW190425_reweighted-m1.csv,m /home/phil/Research/nspop/GW190425_reweighted-m2.csv,m /home/phil/Research/nspop/GW200105_reweighted-m2.csv,m /home/phil/Research/nspop/GW200115_reweighted-m2.csv,m -f 500 -M -v

6. plot-pop-model /home/phil/Research/nspop/bimod_mass.csv -p bimod_mass -b -n 1e4 -r /home/phil/Research/nspop/GW170817_reweighted-m1.csv,m /home/phil/Research/nspop/GW170817_reweighted-m2.csv,m /home/phil/Research/nspop/GW190425_reweighted-m1.csv,m /home/phil/Research/nspop/GW190425_reweighted-m2.csv,m /home/phil/Research/nspop/GW200105_reweighted-m2.csv,m /home/phil/Research/nspop/GW200115_reweighted-m2.csv,m -q 0.05 0.5 0.95 -v

7. extract-pop-param-quantiles /home/phil/Research/nspop/bimod_mass.csv mu sigma -w weight -v

8. plothist /home/phil/Research/nspop/unif_mass.csv /home/phil/Research/nspop/unif_mass.csv -x mmax -W weight -o unif_mass-mmax.png -l "M_max [M_sun]" "" -L "post" "prior" -u True False -v
