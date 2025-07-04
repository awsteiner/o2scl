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
NICKS := ost_min ost_plus ost_py u24.04_py u24.04_min u24.04_plus \
	arch_min
# Dockerfiles
DFILE := ost_min ost_plus ost_py u24.04_py u24.04_min u24.04_plus \
	arch_min

# Rule template
define RULE_tlate

# CPU-only builds

$(1).nc: empty
	sudo docker buildx build . \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		--no-cache \
		--target working 2>&1 | tee out/$(1).out &

$(1).build: empty
	sudo docker buildx build . \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		--target working 2>&1 | tee out/$(1).out &

$(1).check: empty
	sudo docker run --gpus all --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		sh -c "cd /opt/o2scl; make -j 10 check o2scl-examples" \
		| tee out/$(1).check &

$(1).pycheck: empty
	sudo docker run --gpus all --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		sh -c "cd /opt/o2sclpy; make testq" \
		| tee out/$(1).check &

$(1).run:
	sudo docker run --gpus all -it --rm \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE))

$(1).push:
	sudo docker push \
		awsteiner/o2scl:$$(word $(2), $(DFILE))

# CUDA builds

$(1).cuda1:
	sudo docker buildx build . \
		--build-arg \
	BASE=awsteiner/foundation:cuda_12.6_tf_2.18_torch_2.7_m1 \
		--build-arg CONFIG=--enable-cuda \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		--target working 2>&1 | tee out/$(1).cuda1 &

$(1).cuda2:
	sudo docker buildx build . \
		--build-arg \
	BASE=awsteiner/foundation:cuda_12.8_tf_2.18_torch_2.7_m2 \
		--build-arg CONFIG=--enable-cuda \
		-f $$(word $(2), $(DFILE)) \
		-t awsteiner/o2scl:$$(word $(2), $(DFILE)) \
		--target working 2>&1 | tee out/$(1).cuda2 &

endef

# Loop over the indices of the targets
INDEX := $(shell seq 1 $(words $(NICKS)))
$(foreach i,$(INDEX),$(eval $(call RULE_tlate,$(word $(i),$(NICKS)),$(i))))

empty:
