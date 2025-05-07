#!/bin/bash
python3 --version
nvidia-smi
nvcc --version
python3 -c "import tensorflow as tf; print('TF GPU devices:',tf.config.list_physical_devices('GPU')); print('TF build info:',tf.sysconfig.get_build_info())"
