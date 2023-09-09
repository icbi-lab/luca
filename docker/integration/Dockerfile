FROM python:latest

# Install dependencies
RUN pip install \
    scanpy \
    nxfvars \
    seaborn \
    pandas \
    anndata \
    papermill

RUN git clone https://github.com/icbi-lab/luca.git && \
    cd luca/lib && \
    pip install -e ./scanpy_helper_submodule