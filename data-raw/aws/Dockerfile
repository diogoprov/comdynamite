# Environment for the simulation scripts of data-raw/ on an AWS EC2 spot instance
# (data-raw/aws/launch_spot.sh; recipe adapted from crywolf/dev/dyn-sim/aws). The checkout,
# from git HEAD (pack_code.sh), is mounted at /project when the container runs; the scripts load
# the package from it with pkgload::load_all(), so only the dependencies are installed here.
# Build from the package root:  docker build -f data-raw/aws/Dockerfile -t comdynamite-sim .
FROM rocker/r-ver:4.6.1

# One R process per core: multithreaded BLAS would only compete with the other workers.
ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1

RUN R -q -e 'install.packages(c("ape", "pkgload"))'

COPY DESCRIPTION /src/DESCRIPTION
RUN R -q -e 'install.packages("remotes"); remotes::install_deps("/src", dependencies = NA)' \
 && R -q -e 'sessionInfo()'

WORKDIR /project
