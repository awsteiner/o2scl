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
# v0.930a6_ot_min

build_v0.930a6_ot_min:
	sudo docker buildx build . -f v0.930a6_ot_min \
		-t awsteiner/o2scl:v0.930a6_ot_min \
		> v0.930a6_ot_min.out 2>&1

test_v0.930a6_ot_min:
	sudo docker run awsteiner/o2scl:v0.930a6_ot_min sh -c \
		"cd /opt/o2scl && make check && make o2scl-examples" \
		> test_v0.930a6_ot_min.out 2>&1

push_v0.930a6_ot_min:
	sudo docker push awsteiner/o2scl:v0.930a6_ot_min

pull_v0.930a6_ot_min:
	sudo docker pull awsteiner/o2scl:v0.930a6_ot_min

all_v0.930a6_ot_min: build_v0.930a6_ot_min test_v0.930a6_ot_min

# ----------------------------------------------
# v0.930a6_ot_plus

build_v0.930a6_ot_plus:
	sudo docker buildx build . -f v0.930a6_ot_plus \
		-t awsteiner/o2scl:v0.930a6_ot_plus \
		> v0.930a6_ot_plus.out 2>&1

test_v0.930a6_ot_plus:
	sudo docker run awsteiner/o2scl:v0.930a6_ot_plus sh -c \
		"cd /opt/o2scl && make check && make o2scl-examples" \
		> test_v0.930a6_ot_plus.out 2>&1

push_v0.930a6_ot_plus:
	sudo docker push awsteiner/o2scl:v0.930a6_ot_plus

pull_v0.930a6_ot_plus:
	sudo docker pull awsteiner/o2scl:v0.930a6_ot_plus

all_v0.930a6_ot_plus: build_v0.930a6_ot_plus test_v0.930a6_ot_plus

# ----------------------------------------------
# v0.930a6_ot_py

build_v0.930a6_ot_py:
	sudo docker buildx build . -f v0.930a6_ot_py \
		-t awsteiner/o2scl:v0.930a6_ot_py \
		> v0.930a6_ot_py.out 2>&1

test_v0.930a6_ot_py:
	sudo docker run awsteiner/o2scl:v0.930a6_ot_py sh -c \
		"cd /opt/o2scl && make check && make o2scl-examples" \
		> test_v0.930a6_ot_py.out 2>&1

test2_v0.930a6_ot_py:
	sudo docker run awsteiner/o2scl:v0.930a6_ot_py sh -c \
		"cd /opt/o2sclpy && make testq" \
		>> test2_v0.930a6_ot_py.out 2>&1

push_v0.930a6_ot_py:
	sudo docker push awsteiner/o2scl:v0.930a6_ot_py

pull_v0.930a6_ot_py:
	sudo docker pull awsteiner/o2scl:v0.930a6_ot_py

all_v0.930a6_ot_py: build_v0.930a6_ot_py test_v0.930a6_ot_py

# ----------------------------------------------
# v0.930a6_u24.04_min

#		--target working \
#		--no-cache \

build_v0.930a6_u24.04_min:
	sudo docker buildx build . \
		-f v0.930a6_u24.04_min \
		-t awsteiner/o2scl:v0.930a6_u24.04_min \
		> v0.930a6_u24.04_min.out 2>&1

#	sudo docker buildx build . \
#		--target working \
#		--no-cache \

#test1_v0.930a6_u24.04_min:
#	sudo docker run awsteiner/o2scl:v0.930a6_u24.04_min sh -c \
#	"cd /opt/o2scl && make check && make o2scl-examples" \
#		> test1_v0.930a6_u24.04_min.out 2>&1

test2_v0.930a6_u24.04_min:
	sudo docker run awsteiner/o2scl:v0.930a6_u24.04_min sh -c \
		"cd /opt/o2sclpy && make testq" \
		> test2_v0.930a6_u24.04_min.out 2>&1

push_v0.930a6_u24.04_min:
	sudo docker push awsteiner/o2scl:v0.930a6_u24.04_min

pull_v0.930a6_u24.04_min:
	sudo docker pull awsteiner/o2scl:v0.930a6_u24.04_min

all_v0.930a6_u24.04_min: build_v0.930a6_u24.04_min \
	test2_v0.930a6_u24.04_min push_v0.930a6_u24.04_min 

# ----------------------------------------------
# v0.930a6_u24.04_py

build_v0.930a6_u24.04_py:
	sudo docker buildx build . \
		-f v0.930a6_u24.04_py \
		-t awsteiner/o2scl:v0.930a6_u24.04_py \
		> v0.930a6_u24.04_py.out 2>&1

#	sudo docker buildx build . \
#		--no-cache \
#		--target working \

#test1_v0.930a6_u24.04_py:
#	sudo docker run awsteiner/o2scl:v0.930a6_u24.04_py sh -c \
#	"cd /opt/o2scl && make check && make o2scl-examples" \
#		> test1_v0.930a6_u24.04_py.out 2>&1

test2_v0.930a6_u24.04_py:
	sudo docker run awsteiner/o2scl:v0.930a6_u24.04_py sh -c \
		"cd /opt/o2sclpy && make testq" \
		> test2_v0.930a6_u24.04_py.out 2>&1

push_v0.930a6_u24.04_py:
	sudo docker push awsteiner/o2scl:v0.930a6_u24.04_py

pull_v0.930a6_u24.04_py:
	sudo docker pull awsteiner/o2scl:v0.930a6_u24.04_py

all_v0.930a6_u24.04_py: build_v0.930a6_u24.04_py \
	test2_v0.930a6_u24.04_py push_v0.930a6_u24.04_py 

# ----------------------------------------------
# v0.930a6_u24.04_py_mp

build_v0.930a6_u24.04_py_mp:
	sudo docker buildx build . \
		-f v0.930a6_u24.04_py_mp \
		-t awsteiner/o2scl:v0.930a6_u24.04_py_mp \
		> v0.930a6_u24.04_py_mp.out 2>&1

#	sudo docker buildx build . \
#		--no-cache \
#		--target working \

#test1_v0.930a6_u24.04_py_mp:
#	sudo docker run awsteiner/o2scl:v0.930a6_u24.04_py_mp sh -c \
#	"cd /opt/o2scl && make check && make o2scl-examples" \
#		> test1_v0.930a6_u24.04_py_mp.out 2>&1

test2_v0.930a6_u24.04_py_mp:
	sudo docker run awsteiner/o2scl:v0.930a6_u24.04_py_mp sh -c \
		"cd /opt/o2sclpy_mp && make testq" \
		> test2_v0.930a6_u24.04_py_mp.out 2>&1

push_v0.930a6_u24.04_py_mp:
	sudo docker push awsteiner/o2scl:v0.930a6_u24.04_py_mp

pull_v0.930a6_u24.04_py_mp:
	sudo docker pull awsteiner/o2scl:v0.930a6_u24.04_py_mp

all_v0.930a6_u24.04_py_mp: build_v0.930a6_u24.04_py_mp \
	test2_v0.930a6_u24.04_py_mp push_v0.930a6_u24.04_py_mp 

