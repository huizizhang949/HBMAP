FROM rocker/binder:4.5.1
LABEL maintainer='eddy_a'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN wget https://github.com/huizizhang949/HBMAP/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2025-07-16/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
