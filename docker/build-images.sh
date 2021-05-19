#!/bin/bash
set -e
CONTAINER_ENGINE="podman"
CONTAINER_ENGINE_OPTIONS="--cgroup-manager=cgroupfs"
#CONTAINER_ENGINE="docker"
#CONTAINER_ENGINE_OPTIONS=""

${CONTAINER_ENGINE} pull fenicsproject/test-env:mpich
${CONTAINER_ENGINE} pull dolfinx/dev-env
${CONTAINER_ENGINE} build ${CONTAINER_ENGINE_OPTIONS} -t fractional-afem/dolfinx:debug .
${CONTAINER_ENGINE} build ${CONTAINER_ENGINE_OPTIONS} \
   --build-arg CMAKE_BUILD_TYPE=RelWithDebug --build-arg PIP_EXTRA_FLAGS='' \
   --build-arg CXXFLAGS='-O2' --build-arg IMAGE=dolfinx/dev-env -t fractional-afem/dolfinx:release .