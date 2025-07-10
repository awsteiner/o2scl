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
NICKS := ost_min ost_plus ost_py u24.04_py u24.04_min u24.04_plus arch_min
#u25.04_min u25.04_plus
# Dockerfiles
DFILE := ost_min ost_plus ost_py u24.04_py u24.04_min u24.04_plus arch_min
#u25.04_min u25.04_plus
# Architectures
ARCHS := amd64 arm64 i386
PREFIX := v0.931rc

# Rule template
define RULE_tlate

# CPU-only builds

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

ifeq ($(O2SCL_NPROCS), undefined)
$(1).check: empty
	sudo docker run --gpus all --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		sh -c "cd /opt/o2scl && make check && make o2scl-examples" \
		> $(1).check
else
$(1).check: empty
	sudo docker run --gpus all --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		sh -c "cd /opt/o2scl && make -j $(O2SCL_NPROCS) check && make -j $(O2SCL_NPROCS) o2scl-examples" \
		> $(1).check
endif

$(1).pycheck: empty
	sudo docker run --gpus all --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE))_amd64 \
		sh -c "cd /opt/o2sclpy; make testq" \
		> out/$(1).check

$(1).run:
	sudo docker run --gpus all -it --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE))_amd64

$(1).push:
	sudo docker push \
		awsteiner/o2scl:$$(word $(2), $(DFILE))

# i386 targets

$(1)_i386.nc: empty
	sudo docker buildx build . \
		--platform linux/386 \
		--no-cache
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
		--target working > out/$(1)_i386.out 2>&1 

$(1)_i386.build: empty
	sudo docker buildx build . \
		--platform linux/386 \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_i386 \
		--target working > out/$(1)_i386.out 2>&1 

$(1)_i386.run:
	sudo docker run --gpus all -it --rm \
		--platform linux/386 \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE))_i386

$(1)_i386.push:
	sudo docker push \
		awsteiner/o2scl:$$(word $(2), $(DFILE))_i386

# arm64 targets

$(1)_arm64.nc: empty
	sudo docker buildx build . \
		--platform linux/arm64 \
		-f $$(word $(2), $(DFILE)) \
		--no-cache \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
		--target working > out/$(1)_arm64.out 2>&1 

$(1)_arm64.build: empty
	sudo docker buildx build . \
		--platform linux/arm64 \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$(PREFIX)_$$(word $(2), $(DFILE))_arm64 \
		--target working > out/$(1)_arm64.out 2>&1 

$(1)_arm64.run:
	sudo docker run --gpus all -it --rm \
		--platform linux/arm64 \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE))_arm64

$(1)_arm64.push:
	sudo docker push \
		awsteiner/o2scl:$$(word $(2), $(DFILE))_arm64

# CUDA builds

$(1).cuda1:
	sudo docker buildx build . \
		--build-arg \
	BASE=awsteiner/foundation:cuda_12.6_tf_2.18_torch_2.7_m1 \
		--build-arg CONFIG=--enable-cuda \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		--target working 2>&1 > out/$(1).cuda1

$(1).cuda2:
	sudo docker buildx build . \
		--build-arg \
	BASE=awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7_m2 \
		--build-arg CONFIG=--enable-cuda \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
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

all.nc: u24.04_min.nc u24.04_min_arm64.nc \
	u24.04_plus.nc u24.04_plus_arm64.nc \
	ost_min.nc ost_min_arm64.nc ost_min_i386.nc \
	arch_min.nc arch_min_arm64.nc arch_min_i386.nc
	@echo Done

empty:
