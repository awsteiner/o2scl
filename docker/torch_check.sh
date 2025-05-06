#!/bin/bash
python3 --version
nvidia-smi
nvcc --version
python3 -c "import torch"
python3 -c "import torch; print(torch.cuda.is_available())"
