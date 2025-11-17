#!/bin/bash

conda create -n flash_dock python=3.11 -y

conda activate flash_dock

pip install torch==2.9.1

pip install --no-build-isolation 'unicore @ git+ssh://git@github.com/dptech-corp/Uni-Core.git@ace6fae1c8479a9751f2bb1e1d6e4047427bc134'

pip install -r requirements.txt
