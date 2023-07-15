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

docker_show:
	- sudo docker ps --all
	- sudo docker images --all

# ----------------------------------------------
# Standard builds

build_v0.927_ot:
	sudo docker build . -f v0.927_ot \
		> v0.927_ot.out 2>&1 &

build_v0.928_ot:
	sudo docker build . -f v0.928_ot \
		> v0.928_ot.out 2>&1 &

build_v0.927_u22.04_min:
	sudo docker build . -f v0.927_u22.04_min \
		> v0.927_u22.04_min.out 2>&1 &

build_v0.927_u22.04:
	sudo docker build . -f v0.927_u22.04 \
		> v0.927_u22.04.out 2>&1 &

build_v0.928_u23.04_min:
	sudo docker build . -f v0.928_u23.04_min \
		> v0.928_u23.04_min.out 2>&1 &

build_v0.928_u23.04:
	sudo docker build . -f v0.928_u23.04 \
		> v0.928_u23.04.out 2>&1 &

build_v0.928_u23.04_py:
	sudo docker build . -f v0.928_u23.04_py \
		> v0.928_u23.04_py.out 2>&1 &

build_all:
	sudo docker build . -f v0.927_ot \
		> v0.927_ot.out 2>&1 
	sudo docker build . -f v0.928_ot \
		> v0.928_ot.out 2>&1 
	sudo docker build . -f v0.927_u22.04_min \
		> v0.927_u22.04_min.out 2>&1 
	sudo docker build . -f v0.927_u22.04 \
		> v0.927_u22.04.out 2>&1 
	sudo docker build . -f v0.928_u23.04_min \
		> v0.928_u23.04_min.out 2>&1 
	sudo docker build . -f v0.928_u23.04 \
		> v0.928_u23.04.out 2>&1 
	sudo docker build . -f v0.928_u23.04_py \
		> v0.928_u23.04_py.out 2>&1 

# ----------------------------------------------
# Old builds

build_v0.928_ol15.5_min:
	sudo docker build . -f v0.928_ol15.5_min \
		> v0.928_ol15.5_min.out 2>&1 &

build_v0.928_ol15.5:
	sudo docker build . -f v0.928_ol15.5 \
		> v0.928_ol15.5.out 2>&1 &
