# README

Offsite is an offline autotuning framework, currently specialized on explicit ODE methods. Offsite can automatically derive implementation variants from user-provided abstract YAML description formats and then rates these variants based on the analytic ECM performance model. Currently supports x86 architectures: Intel and AMD. Two different backends are currently supported to derived the required ECM performance model contributions: (1) [kerncraft](https://github.com/RRZE-HPC/kerncraft) and (2) [YaskSite](https://github.com/seasite-project/YaskSite).


### How to install? ###

See [install.sh](https://github.com/seasite-project/Offsite/blob/master/install.sh) for a sample script to install Offsite. In general, the following steps are required to setup Offsite:

* `pip install --user wheel` (Make sure wheel is installed)
* `git clone https://github.com/RRZE-HPC/kerncraft && cd kerncraft && git checkout v0.8.5` 
* `python3 setup.py bdist_wheel && pip install --user dist/kerncraft*.whl` (Replace * with the used kerncraft version; current Offsite version requires kerncraft version 0.8.5)
* `iaca_get --I-accept-the-Intel-What-If-Pre-Release-License-Agreement-and-please-take-my-soul` (IACA is required by kerncraft)
* `git clone https://github.com/seasite-project/Offsite && cd Offsite`
* `python3 setup.py bdist_wheel && pip install --user dist/offsite*.whl` (Replace * with the used Offsite version)


### External dependencies ###

* Intel [IACA](https://software.intel.com/en-us/articles/intel-architecture-code-analyzer)
* [kerncraft](https://github.com/RRZE-HPC/kerncraft)


### How to run? ###

* `offsite_tune -h (to get help and options)`
* `Some example tuning scenarios can be found in the ` [examples](https://github.com/seasite-project/Offsite/tree/master/examples) folder.
* `Let's try for instance IVP InverterChain on a Cascade Lake CPU with backend kerncraft for a fixed ODE size of n=50,000:`
	
	`offsite_tune --machine examples/machines/CascadelakeSP_Gold-6248.yml --compiler icc --impl examples/impls/pirk/ --kernel examples/kernels/pirk/ --method examples/methods/implicit/radauIIA7.ode --ivp examples/ivps/InverterChain.ivp --bench examples/bench/OMP_BARRIER_CascadelakeSP_Gold-6248_icc19.0.2.187.bench --tool kerncraft -n 50000 --verbose`
* `For an example on using YaskSite as backend follow this ` [link](https://github.com/seasite-project/SC20_YASKSITE_AD/blob/master/README.md).
