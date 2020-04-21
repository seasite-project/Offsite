# README

Offsite is an offline autotuning framework, currently specialized on explicit ODE methods. Offsite can automatically derive implementation variants from user-provided abstract YAML description formats and then rates these variants based on the analytic ECM performance model. Currently supports x86 architectures: Intel and AMD. Two different backends are currently supported to derived the required ECM performance model contributions: (1) [kerncraft](https://github.com/RRZE-HPC/kerncraft) (2) [YaskSite](https://github.com/seasite-project/YaskSite)


### Features ###

* T
* O
* D
* O


### How to install? ###

See [install.sh](https://github.com/seasite-project/Offsite/blob/master/install.sh) for a sample script to install Offsite. In general, the following steps are required to setup Offsite:

* `pip install --user wheel` (Make sure wheel is installed)
* `git clone https://github.com/RRZE-HPC/kerncraft && cd kerncraft` 
* `python3 setup.py bdist_wheel && pip install --user dist/kerncraft*.whl` (Replace * with the used kerncraft version; newest Offsite version requires kerncraft version 0.8.4)
* `iaca_get --I-accept-the-Intel-What-If-Pre-Release-License-Agreement-and-please-take-my-soul` (IACA is required by kerncraft)
* `git clone https://github.com/seasite-project/Offsite && cd Offsite`
* `python3 setup.py bdist_wheel && pip install --user dist/offsite*.whl` (Replace * with the used Offsite version)


### External dependencies ###

* Intel [IACA](https://software.intel.com/en-us/articles/intel-architecture-code-analyzer)
* [kerncraft](https://github.com/RRZE-HPC/kerncraft)


### How to run? ###

* T
* O
* D
* O
