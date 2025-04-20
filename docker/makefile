help:
	@echo "Docker targets (not part of autoconf package):"
	@echo "---------------------------------------------------"
	@echo "docker_clean:"
	@echo "docker_show:"
	@echo "build_opensuse_leap:"
	@echo "build_opensuse_leap_notest:"
	@echo "build_opensuse_leap_dev:"
	@echo "build_ubuntu:"
	@echo "build_ubuntu_notest:"
	@echo "build_ubuntu_dev_full:"
	@echo "build_ubuntu_dev_min:"

docker_clean:
	-sudo docker rm \
	`sudo docker ps --all | grep -i -v container | awk '{print $$1}'`
	-sudo docker rmi \
	`sudo docker images --all | grep -i -v container | awk '{print $$3}'`

docker_clear_cache:
	-sudo docker 

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

testq:
	sudo docker run --rm -it \
		--env="DISPLAY" \
		--volume="/etc/group:/etc/group:ro" \
		--volume="/etc/passwd:/etc/passwd:ro" \
		--volume="/etc/shadow:/etc/shadow:ro" \
		--volume="/etc/sudoers.d:/etc/sudoers.d:ro" \
		--volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" \
		awsteiner/o2scl:v0.927_u22.04_min

#		--user="awsteiner" 

# ----------------------------------------------

VERSIONS := 0.930 0.931a1 0.931a2

OSS := u24.04 arch ot

BTYPES := min py py_mp

TEST_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).test)))

CUDA_TARGETS := $(foreach vvar, $(VERSIONS), \
$(foreach ovar, $(OSS), \
$(foreach bvar, $(BTYPES), \
v$(vvar)_$(ovar)_$(bvar).cuda)))

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

$(CUDA_TARGETS) : %.cuda :
	sudo docker buildx build . \
		--build-arg BASE=nvidia/cuda:12.6.1-devel-ubuntu24.04 \
		-f $* -t awsteiner/o2scl:$*_cuda \
		--target working \
		> $*.cuda.out 2>&1 

$(CUDA_TEST1_TARGETS) : %.cuda.test1 :
	sudo docker run awsteiner/o2scl:$*_cuda --gpus all \
		sh -c "cd /opt/o2scl && make check && make o2scl-examples" \
		> test1_$*.out 2>&1 

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
