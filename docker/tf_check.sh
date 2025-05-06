#!/bin/bash
python3 --version
nvidia-smi
nvcc --version
python3 -c "import tensorflow"
python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
python -c 'import tensorflow as tf; print(tf.sysconfig.get_build_info())'
