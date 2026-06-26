# podman build -f Dockerfile --format docker --build-arg CPPBMAD_REF=$(git rev-parse HEAD) -t pybmad-build .
# podman run --rm -v $PWD:/workspace pybmad-build
#
# Stage 1: Build pybmad from source (bmad + cppbmad + pybmad wheel)
FROM condaforge/miniforge3:latest AS builder

WORKDIR /build

# Clone Bmad
ARG BMAD_REF=main
RUN git clone --depth=1 --branch=${BMAD_REF} \
    https://github.com/bmad-sim/bmad-ecosystem.git bmad_src

RUN mamba env create -n bmad-build -f bmad_src/.github/bmad-build-env.yaml && \
    mamba install -n bmad-build -y libacl libgomp cmake pydantic clang-format

SHELL ["conda", "run", "-n", "bmad-build", "/bin/bash", "-c"]

ENV ACC_ROOT_DIR=/build/bmad_src
ENV USE_MPI=N
ENV SHARED=Y
ENV USE_CONDA=1

RUN echo "export ACC_SET_GMAKE_JOBS=8" >> bmad_src/util/dist_prefs
RUN cd bmad_src && .github/scripts/install_bmad.sh

# Clone cppbmad at a specific commit
ARG CPPBMAD_REF=master
RUN git clone https://github.com/ken-lauer/cppbmad.git cppbmad && \
    cd cppbmad && git checkout ${CPPBMAD_REF}

WORKDIR /build/cppbmad

RUN python -m codegen
RUN cmake -DCMAKE_BUILD_TYPE=release -B release . && cmake --build release
RUN CMAKE_BUILD_PARALLEL_LEVEL=1 python -m pip wheel --no-deps -w dist .
RUN python -m pip install dist/*.whl
RUN python -c "import pybmad; print('pybmad loaded OK')"
RUN cp dist/*.whl /tmp/pybmad.whl

# Stage 2: Docs build (mirrors CI build-docs action)
FROM condaforge/miniforge3:latest AS docs

WORKDIR /workspace

COPY --from=builder /build/bmad_src/.github/bmad-build-env.yaml /tmp/env.yaml
RUN mamba env create -n bmad-build -f /tmp/env.yaml && \
    mamba install -n bmad-build -y libacl libgomp

SHELL ["conda", "run", "-n", "bmad-build", "/bin/bash", "-c"]

COPY --from=builder /tmp/pybmad.whl /tmp/pybmad-0.0.0-py3-none-any.whl
COPY --from=builder /build/bmad_src/production/lib /opt/bmad-lib

ENV LD_LIBRARY_PATH=/opt/bmad-lib

RUN python -m pip install /tmp/pybmad-0.0.0-py3-none-any.whl

# Install docs dependencies (same as CI action)
RUN python -m pip install mkdocs mkdocstrings[python] mkdocs-autorefs pymdown-extensions nanobind

RUN python -c "import pybmad; print('pybmad OK')"

ARG CPPBMAD_REF=master
RUN git clone https://github.com/ken-lauer/cppbmad.git /workspace && \
    cd /workspace && git checkout ${CPPBMAD_REF}

CMD ["conda", "run", "-n", "bmad-build", "mkdocs", "build", "-f", "python/docs/mkdocs.yml", "-d", "/workspace/docs-site"]
