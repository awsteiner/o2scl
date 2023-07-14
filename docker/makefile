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

build_opensuse_leap:
	sudo docker build . -f opensuse_leap \
		> opensuse_leap.out 2>&1 &

build_opensuse_tumble:
	sudo docker build . -f opensuse_tumble \
		> opensuse_tumble.out 2>&1 &

build_opensuse_leap_notest:
	sudo docker build . -f opensuse_leap_notest \
		> opensuse_leap_notest.out 2>&1 &

build_opensuse_leap_dev:
	sudo docker build . -f opensuse_leap_dev \
		> opensuse_leap_dev.out 2>&1 &

build_ubuntu:
	sudo docker build . -f ubuntu \
		> ubuntu.out 2>&1 &

build_ubuntu_notest:
	sudo docker build . -f ubuntu_notest \
		> ubuntu_notest.out 2>&1 &

build_ubuntu_dev_full:
	sudo docker build . -f ubuntu_dev_full \
		> ubuntu_dev_full.out 2>&1 &

build_ubuntu_dev_min:
	sudo docker build . -f ubuntu_dev_min \
		> ubuntu_dev_min.out 2>&1 &

build_v0.928_ol15.5_min:
	sudo docker build . -f v0.928_ol15.5_min \
		> v0.928_ol15.5_min.out 2>&1 &
