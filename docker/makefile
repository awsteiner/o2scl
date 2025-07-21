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

# ----------------------------------------------------------------------

# Nicknames
NICKS := ost_min ost_plus ost_py u24.04_py u24.04_min u24.04_plus arch_min \
u25.04_min u25.04_plus arch_plus arch_py
# Dockerfiles
DFILE := ost_min ost_plus ost_py u24.04_py u24.04_min u24.04_plus arch_min \
u25.04_min u25.04_plus arch_plus arch_py
# Architectures
ARCHS := amd64 arm64 i386
PREFIX := v0.931rc

# Rule template
define RULE_tlate

# CPU-only builds

ifeq ($(O2SCL_NPROCS), undefined)

$(1).nc: empty
	sudo docker buildx build . \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		--no-cache \
		--target working 2>&1 > \
		out/$(1).out 2>&1 

$(1).build: empty
	sudo docker buildx build . \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		--target working > out/$(1).out 2>&1 

$(1).check: empty
	sudo docker run --rm \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		sh -c "cd /opt/o2scl && make check && make o2scl-examples" \
		> out/$(1).check

else

$(1).nc: empty
	sudo docker buildx build . \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		--build-arg NPROCS=$(O2SCL_NPROCS) \
		--no-cache \
		--target working 2>&1 > \
		out/$(1).out 2>&1 

$(1).build: empty
	sudo docker buildx build . \
		-f $$(word $(2), $(DFILE)) \
		--build-arg NPROCS=$(O2SCL_NPROCS) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		--target working > out/$(1).out 2>&1 

$(1).check: empty
	sudo docker run --rm \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		sh -c "cd /opt/o2scl && make -j $(O2SCL_NPROCS) check && make -j $(O2SCL_NPROCS) o2scl-examples" \
		> out/$(1).check

endif

$(1).pycheck: empty
	sudo docker run --rm \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64 \
		sh -c "cd /opt/o2sclpy; make testq" \
		> out/$(1).pycheck

# You might need to run xhost +local: on your host beforehand to allow
# Docker access to your display.
$(1).run:
	sudo docker run -it --rm \
		-e DISPLAY=$(DISPLAY) \
		-v /tmp/.X11-unix:/tmp/.X11-unix \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64

$(1).pull:
	sudo docker pull \
		awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64

$(1).push:
	sudo docker push \
		awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_amd64

# i386 targets

ifeq ($(O2SCL_NPROCS), undefined)

$(1)_i386.nc: empty
	sudo docker buildx build . \
		--platform linux/386 \
		--no-cache \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
		--target working > out/$(1)_i386.out 2>&1 

$(1)_i386.build: empty
	sudo docker buildx build . \
		--platform linux/386 \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
		--target working > out/$(1)_i386.out 2>&1 

$(1)_i386.check: empty
	sudo docker run --rm \
		--platform linux/386 \
                -t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
                sh -c "cd /opt/o2scl && make check && make o2scl-examples" \
		> out/$(1)_i386.check
else

$(1)_i386.nc: empty
	sudo docker buildx build \
		--platform linux/386 \
		--no-cache \
		--build-arg NPROCS=$(O2SCL_NPROCS) \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
		--target working . > out/$(1)_i386.out 2>&1 

$(1)_i386.build: empty
	sudo docker buildx build \
		--platform linux/386 \
		--build-arg NPROCS=$(O2SCL_NPROCS) \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
		--target working . > out/$(1)_i386.out 2>&1 

$(1)_i386.check: empty
	sudo docker run --rm \
		--platform linux/386 \
                -t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
                sh -c "cd /opt/o2scl && make -j $(O2SCL_NPROCS) check && make -j $(O2SCL_NPROCS) o2scl-examples" \
		> out/$(1)_i386.check

endif

$(1)_i386.run:
	sudo docker run -it --rm \
		--platform linux/386 \
		-e DISPLAY=$(DISPLAY) \
		-v /tmp/.X11-unix:/tmp/.X11-unix \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386

$(1)_i386.pull:
	sudo docker pull \
		awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386

$(1)_i386.push:
	sudo docker push \
		awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386

# arm64 targets

ifeq ($(O2SCL_NPROCS), undefined)

$(1)_arm64.nc: empty
	sudo docker buildx build \
		--platform linux/arm64 \
		-f $$(word $(2), $(DFILE)) \
		--no-cache \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
		--target working . > out/$(1)_arm64.out 2>&1 

$(1)_arm64.build: empty
	sudo docker buildx build \
		--platform linux/arm64 \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
		--target working . > out/$(1)_arm64.out 2>&1 

$(1)_arm64.check: empty
	sudo docker run --rm \
		--platform linux/arm64 \
                -t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
                sh -c "cd /opt/o2scl && make check && make o2scl-examples" \
		> out/$(1)_arm64.check

else

$(1)_arm64.nc: empty
	sudo docker buildx build . \
		--platform linux/arm64 \
		-f $$(word $(2), $(DFILE)) \
		--build-arg NPROCS=$(O2SCL_NPROCS) \
		--no-cache \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
		--target working > out/$(1)_arm64.out 2>&1 

$(1)_arm64.build: empty
	sudo docker buildx build . \
		--platform linux/arm64 \
		-f $$(word $(2), $(DFILE)) \
		--build-arg NPROCS=$(O2SCL_NPROCS) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
		--target working > out/$(1)_arm64.out 2>&1 

$(1)_arm64.check: empty
	sudo docker run --rm \
		--platform linux/arm64 \
                -t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
                sh -c "cd /opt/o2scl && make -j $(O2SCL_NPROCS) check && make -j $(O2SCL_NPROCS) o2scl-examples" \
		> out/$(1)_arm64.check

endif

$(1)_arm64.run:
	sudo docker run -it --rm \
		--platform linux/arm64 \
		-e DISPLAY=$(DISPLAY) \
		-v /tmp/.X11-unix:/tmp/.X11-unix \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64

$(1)_arm64.pull:
	sudo docker pull \
		awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64

$(1)_arm64.push:
	sudo docker push \
		awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64

# CUDA builds

$(1).cuda1:
	sudo docker buildx build . \
		--build-arg \
	BASE=awsteiner/foundation:cuda_12.6_tf_2.18_torch_2.7_m1 \
		--build-arg CONFIG=--enable-cuda \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE)) \
		--target working 2>&1 > out/$(1).cuda1

$(1).cuda2:
	sudo docker buildx build . \
		--build-arg \
	BASE=awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7_m2 \
		--build-arg CONFIG=--enable-cuda \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE)) \
		--target working 2>&1 > out/$(1).cuda2

endef

# Loop over the indices of the targets
INDEX := $(shell seq 1 $(words $(NICKS)))

$(foreach i,$(INDEX),$(eval $(call RULE_tlate,$(word $(i),$(NICKS)),$(i))))

# ----------------------------------------------------------------------

# Note that the ubuntu:24.04 manifest does not include i386 support

all: u24.04_min.build u24.04_min_arm64.build \
	u24.04_plus.build u24.04_plus_arm64.build \
	ost_min.build ost_min_arm64.build ost_min_i386.build \
	arch_min.build arch_min_arm64.build arch_min_i386.build
	@echo Done

#	-make u24.04_min_arm64.nc
#	-make u24.04_plus_arm64.nc
#	-make ost_min_arm64.nc
#	-make arch_min_arm64.nc
#	-make arch_min_i386.nc

all.nc: empty
	-make u24.04_min.nc
	-make u24.04_plus.nc
	-make u24.04_py.nc
	-make ost_min.nc
	-make ost_min_i386.nc
	-make arch_min.nc
	-make arch_py.nc

empty:
