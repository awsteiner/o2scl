#!/bin/bash
python3 --version
nvidia-smi
nvcc --version
python3 -c "import torch; print('CUDA avail:',torch.cuda.is_available()); print('CUDA version:',torch.version.cuda); print('CUDA built:',torch.backends.cuda.is_built()); print('CUDNN version:',torch.backends.cudnn.version()); print('CUDA device name:',torch.cuda.get_device_name(0))"
