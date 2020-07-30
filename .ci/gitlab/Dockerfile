ARG BASE
FROM dunecommunity/$BASE

COPY . ${SUPERDIR}/${MY_MODULE}/
RUN ls -l ${SUPERDIR}/${MY_MODULE}/
RUN cd ${SUPERDIR} && \
    echo ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} all && \
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} all
