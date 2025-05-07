#!/bin/bash
echo "PATH:"
echo $PATH
echo ""
echo "LD_LIBRARY_PATH:"
echo $LD_LIBRARY_PATH
echo ""
echo "Python version:"
python3 --version
echo ""
echo "Output of 'nvidia-smi':"
nvidia-smi
echo ""
echo "NVCC version:"
nvcc --version
echo ""
echo "Torch:"
python3 -c "import torch; print('CUDA avail:',torch.cuda.is_available()); print('CUDA version:',torch.version.cuda); print('CUDA built:',torch.backends.cuda.is_built()); print('CUDNN version:',torch.backends.cudnn.version()); print('CUDA device name:',torch.cuda.get_device_name(0))"
