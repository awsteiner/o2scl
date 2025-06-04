help:
	@echo "Docker targets (not part of autoconf package):"
	@echo "---------------------------------------------------"
	@echo "docker_clean:"
	@echo "docker_clean2:"
	@echo "docker_show:"
	@echo "docker_stop:"

docker_clean:
	-sudo docker rm \
	`sudo docker ps --all | grep -i -v container | awk '{print $$1}'`
	-sudo docker rmi \
	`sudo docker images --all | grep -i -v container | awk '{print $$3}'`

docker_show:
	- sudo docker ps --all
	- sudo docker images --all
	- sudo docker system df

docker_make_cache:
	sudo docker create -v /mnt/ccache:/ccache --name ccache debian

docker_mount_cache:
	sudo docker run -e CCACHE_DIR=/ccache --volumes-from ccache -it debian

DOCKER_PS_Q = $(shell sudo docker ps -q)

docker_stop:
	- sudo docker stop $(DOCKER_PS_Q)

docker_clean2:
	- sudo docker image prune -a
	- sudo docker stop $(sudo docker ps -q)
	- sudo docker container prune
	- sudo docker buildx prune -f
	- sudo docker system df

# ----------------------------------------------

VERSIONS := 0.930 0.931a1 0.931a2

OSS := u24.04 arch ot

BTYPES := min py

TEST_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).test)))

CUDA1_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).cuda1)))

CUDA2_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).cuda2)))

CPU_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).cpu)))

MP_CUDA_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).mp.cuda)))

CUDA_TEST1_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).cuda.test1)))

DOCKERFILES := $(shell ls v* | grep -v req)

dfl:
	@echo $(DOCKERFILES)

list:
	@echo $(TEST_TARGETS)

$(TEST_TARGETS) : %.test :
	@echo "test"
	@echo $@
	@echo $*

$(CPU_TARGETS) : %.cpu :
	sudo docker buildx build . \
		--build-arg BASE=awsteiner/foundation:tf_2.18_torch_2.7 \
		--no-cache -f $* -t awsteiner/o2scl:$* \
		--target working \
		> $*.cpu.out 2>&1 

$(CUDA1_TARGETS) : %.cuda1 :
	sudo docker buildx build . \
		--build-arg \
		BASE=awsteiner/foundation:cuda_12.6_tf_2.18_torch_2.7_m1 \
		--build-arg CONFIGX=--enable-cuda \
		-f $* -t awsteiner/o2scl:$*_cuda1 \
		--target working \
		> $*.cuda1.out 2>&1 

$(CUDA2_TARGETS) : %.cuda2 :
	sudo docker buildx build . \
		--build-arg \
		BASE=awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7_m2 \
		--build-arg CONFIGX=--enable-cuda \
		-f $* -t awsteiner/o2scl:$*_cuda2 \
		--no-cache \
		--target working \
		> $*.cuda2.out 2>&1 

$(MP_CUDA_TARGETS) : %.mp.cuda :
	sudo docker buildx build . \
		--build-arg BASE=nvidia/cuda:12.6.1-devel-ubuntu24.04 \
		--build-arg CONFIG=--enable-multip \
		-f $* -t awsteiner/o2scl:$*_cuda \
		--target working \
		> $*.cuda.out 2>&1 

$(CUDA_TEST1_TARGETS) : %.cuda.test1 :
	sudo docker run awsteiner/o2scl:$*_cuda --gpus all \
		sh -c "cd /opt/o2scl && make check && make o2scl-examples" \
		> test1_$*.out 2>&1 

tt2.build:
	sudo docker buildx build . \
		-f tf_2.18_torch_2.7 \
		-t awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7 \
		--target working

tt2.run:
	sudo docker run --gpus all \
		-t awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7 \
		sh -c "cd /opt; ./tf_check.sh; ./torch_check.sh"

tt2.runb:
	sudo docker run -it --gpus all \
		-t awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7 

tt.build:
	sudo docker buildx build . \
		-f cuda_12.6_tf_2.18_torch_2.7 \
		-t awsteiner/foundation:cuda_12.6_tf_2.18_torch_2.7 \
		--target working

tt.run:
	sudo docker run --gpus all \
		-t awsteiner/foundation:cuda_12.6_tf_2.18_torch_2.7 \
		sh -c "cd /opt; ./tf_check.sh; ./torch_check.sh"

torch_2.6.out:
	sudo docker buildx build . \
		-f torch_2.6 -t awsteiner/foundation:torch_2.6 \
		--target working

# $(TARGETS) : %.build.cuda : 
# 	@echo $@
# 	@echo $*

# $(TARGETS) : %.build : 
# 	@echo $@
# 	@echo $*

# %.build:
# 	sudo docker buildx build . \
# 		-f $* -t awsteiner/o2scl:$* \
# 		> $*.out 2>&1 

# %.test1:
# 	sudo docker run awsteiner/o2scl:$* sh -c \
# 		"cd /opt/o2scl && make check && make o2scl-examples" \
# 		> test1_$*.out 2>&1 

# %.test2:
# 	sudo docker run awsteiner/o2scl:$* sh -c \
# 		"cd /opt/o2sclpy_mp && make testq" \
# 		> test2_$*.out 2>&1 

# %.push:
# 	sudo docker push awsteiner/o2scl:$*

# %.pull:
# 	sudo docker pull awsteiner/o2scl:$*

# %.all: %.build %.test1 %.test2

# ----------------------------------------------
# cuda_test3

build_ct3:
	sudo docker buildx build -t cuda_test3 -f cuda_test3 .

test_ct3:
	sudo docker run --rm --gpus all cuda_test3

run_ct3:
	sudo docker run --rm -it --gpus all cuda_test3 /bin/bash

