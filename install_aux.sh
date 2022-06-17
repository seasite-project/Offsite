#!/bin/bash

cur_path=${PWD}

# Default directories
install_dir="$HOME/installkit"

# Set PATH to include installed Python3 scripts.
export PATH="${PATH}:~/.local/bin/"

# Install build package required by setup.py.
pip install --user build

# Install wheel package required by setup.py.
pip install --user wheel

# Switch to installation directory.
mkdir -p ${install_dir} 
cd ${install_dir}

# Install kerncraft.
git clone https://github.com/RRZE-HPC/kerncraft && cd kerncraft && git checkout v0.8.13 && python -m build -nw && pip install --user dist/kerncraft*.whl

# Install IACA (required by kerncraft).
iaca_get --I-accept-the-Intel-What-If-Pre-Release-License-Agreement-and-please-take-my-soul

# Install Offsite.
cd ${install_dir}
git clone https://github.com/seasite-project/Offsite && cd Offsite && python -m build -nw && pip install --user dist/offsite*.whl

# Install Offsite auxiliary apps.
cd offsite_aux && python -m build -nw && pip install --user dist/offsite*.whl

# Attention!
echo ""
echo "Note: Add ~/.local/bin to your PATH variable to run the installed Python3 scripts."

cd ${cur_path}
