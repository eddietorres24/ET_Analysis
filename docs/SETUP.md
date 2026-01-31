# Setup (WSL Ubuntu + Miniforge + ET_Analysis)
This repo is designed to run on a local machine using *WSL Ubuntu* and *conda environments via Miniforge*.  
Large files (raw data, outputs) are not tracked by git.

## Assumptions
- OS: WSL Ubuntu (local machine)
- Conda distribution: Miniforge installed at `~/miniforge3`
- Repo location (recommended): `~/repos/ET_Analysis`
- Data is stored on a Windows drive (e.g., D:) and accessed from WSL via a symlink
- Outputs are written to `results/` (gitignored)

#Clone the repo (WSL)
From WSL:
```bash
mkdir -p ~/repos
cd ~/repos
git clone git@github.com:eddietorres24/ET_Analysis.git
cd ET_Analysis

#How to export env settings
conda activate env_name
conda env export --from-history > envs/env_name.from_history.yml

