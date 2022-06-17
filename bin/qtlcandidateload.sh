#!/bin/sh

#
# This script is a wrapper around the process that loads 
# QTL Candidate Gene relationships
#
#
#     qtlcandidateload.sh 
#

if [ -z ${MGICONFIG} ]
then
        MGICONFIG=/usr/local/mgi/live/mgiconfig
        export MGICONFIG
fi

. ${MGICONFIG}/master.config.sh

CONFIG=${TSSGENELOAD}/qtlcandidateload.config

#
# verify & source the configuration file
#

if [ ! -r ${CONFIG} ]
then
    echo "Cannot read configuration file: ${CONFIG}"
    exit 1
fi

. ${CONFIG}

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
#

preload ${OUTPUTDIR}

#
# rm all files/dirs from OUTPUTDIR
#

cleanDir ${OUTPUTDIR}

#
# run the load
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Run qtlcandidateload.py"  | tee -a ${LOG_DIAG}
${PYTHON} ${QTLCANDIDATELOAD}/bin/qtlcandidateload.py | tee -a ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "${QTLCANDIDATELOAD}/bin/qtlcandidateload.py"

# run postload cleanup and email logs

shutDown

