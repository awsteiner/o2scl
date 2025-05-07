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
echo "TensorFlow:"
python3 -c "import tensorflow as tf; print(' '); print('tensorflow config:'); print('TF GPU devices:',tf.config.list_physical_devices('GPU')); print('TF build info:',tf.sysconfig.get_build_info())"
