#  Purpose:
#
#      Find QTL Candidate Genes via Mapping Experiments
#       and load QTL/Candidate Gene relationships
#
#  Inputs:
#
#       QTLs in the database
#	Mapping experiments for QTL with type 'TEXT-QTL-Candidate-Genes'
#               in the database
#
#  Outputs:
#
#	1 BCP file:
#	A pipe-delimited file:
#       	MGI_Relationship.bcp
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#      2:  bcp fails
#
#  Assumes:
#
#	
#  Implementation:
#
#      1) Validate the arguments to the script.
#      2) Perform initialization steps.
#      3) Query DB for QTLs and their Mapping Experiments
#      4) Determine Candidate Genes
#      5) Write out to relationship bcp
#      6) Delete existing relationships
#      7) BCP in new relationships:
#
# History:
#
# sc	06/16/2022
#	- WTS2-896
#

import sys
import os
import string
import db
import mgi_utils

#db.setTrace()

CRT = '\n'
TAB = '\t'

# from configuration file
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
outputDir = os.environ['OUTPUTDIR']
reportDir = os.environ['RPTDIR']
bcpFile = 'MGI_Relationship.bcp'
relationshipFile = '%s/%s' % (outputDir, bcpFile)

cdate  = mgi_utils.date("%m/%d/%Y")
fpRelationship = ''

# The category key 'qtl_to_candidate_gene'
catKey = 1009

# the relationship term key 'Candidate Gene'
relKey = 105563920

# the qualifier key 'Not Specified'
qualKey = 11391898

# the evidence key 'Not Specified'
evidKey = 17396909

# qtl candidate geneload user key
userKey = 1631

# database primary keys, will be set to the next available from the db
nextRelationshipKey = 1000	# MGI_Relationship._Relationship_key

#
# Purpose: open files, create db connection, gets max keys from the db
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables, exits if a file can't be opened,
#
def init():

    global nextRelationshipKey

    #
    # Open input and output files
    #
    openFiles()

    #
    # create database connection
    #
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    #
    # get next MGI_Relationship key
    #
    results = db.sql('''select nextval('mgi_relationship_seq') as nextKey''', 'auto')
    if results[0]['nextKey'] is None:
        nextRelationshipKey = 1000
    else:
        nextRelationshipKey = results[0]['nextKey']

    return 0

# end init() -------------------------------

# Purpose: Open input/output files.
# Returns: exit 1 if cannot open input or output file
# Assumes: Nothing
# Effects: sets global variables, exits if a file can't be opened
#
def openFiles ():

    global fpRelationship

    try:
        fpRelationship = open(relationshipFile, 'w')
    except:
        print(('Cannot open relationships bcp file: %s' % relationshipFile))
        sys.exit(1)

    return 0

# end openFiles() -------------------------------

#
# Purpose: Close the files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
#
def closeFiles ():

    global fpRelationship

    fpRelationship.close()

    return 0

# end closeFiles() -------------------------------

# Purpose: Determines Candidate genes for QTLs, writes to relationship bcp
# Returns: Nothing
# Assumes: file descriptors have been initialized
# Effects: sets global variables, writes to the file system
#
def findRelationships( ): 
    global nextRelationshipKey

    db.sql('''-- get all the official qtl markers
        select distinct _marker_key as qtl_key, symbol as qtl_symbol
        into temporary table temp1
        from mrk_marker
        where _marker_type_key = 6
        and _marker_status_key = 1''', None)

    db.sql('''create index idx1 ON temp1 (qtl_key)''', None)

    db.sql('''-- get qtl marker with mapping type TEXT-QTL-Candidate Genes
    select distinct t1.qtl_key, t1.qtl_symbol, memv1._expt_key
    into temporary table temp2
    from temp1 t1, mld_expt_marker_view memv1
    where memv1.expttype = 'TEXT-QTL-Candidate Genes'
    and memv1._marker_key = t1.qtl_key''', None)

    db.sql('''create index idx2 ON temp2 (_expt_key)''', None)

    db.sql('''-- get the set of markers for each experiment
    select distinct t2.qtl_key, t2.qtl_symbol, t2._expt_key, 
        memv2._marker_key as gene_key, memv2.symbol as gene_symbol, 
        mev._refs_key, mev.jnumid
    into temporary table temp3
    from temp2 t2, mld_expt_marker_view memv2, mld_expt_view mev
    where memv2._expt_key = t2._expt_key
    and memv2._expt_key = mev._expt_key''', None)

    db.sql('''create index idx3 ON temp3 (qtl_key, gene_key)''', None)

    db.sql('''-- get all Gene, Pseudogene, Other Genome Feature or Complex/Cluster/Region
    select distinct m._marker_key as gene_key, m.symbol as gene_symbol, 
        m._marker_type_key, mt.name as non_qtl_marker_type
    into temporary table temp4
    from mrk_marker m, mrk_types mt
    where m._marker_type_key in (1, 7, 9, 10)
    and m._marker_status_key = 1
    and mt._marker_type_key = m._marker_type_key''', None)

    db.sql('''create index idx4 ON temp4 (gene_key, non_qtl_marker_type)''', None)

    results = db.sql('''-- get the set of qtl and non-qtl markers for Candidate gene mapping type
    select distinct t3.qtl_key, t3.qtl_symbol, t3._expt_key, t4.gene_key, 
        t4.gene_symbol, t4.non_qtl_marker_type, t3._refs_key, t3.jnumid
    from temp3 t3, temp4 t4
    where t3.gene_key = t4.gene_key''', 'auto')

    for r in results:        
        objKey1 = r['qtl_key']
        objKey2 = r['gene_key']
        refsKey = r['_refs_key']
        # MGI_Relationship
        fpRelationship.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' % \
            (nextRelationshipKey, catKey, objKey1, objKey2, relKey, qualKey, evidKey, refsKey, userKey, userKey, cdate, cdate))

        nextRelationshipKey += 1
    
    return 0

# end findRelationships() -------------------------------------

# Purpose: deletes existing relationships
# Returns: None
# Assumes: None
# Effects: None
#
def writeReports():

    return 0

# Purpose: deletes existing relationships
# Returns: None
# Assumes: None
# Effects: None
#
def doDeletes():
    db.sql('''delete from MGI_Relationship where _CreatedBy_key = %s ''' % userKey, None)
    db.commit()
    db.useOneConnection(0)

    return 0

# end doDeletes() -------------------------------------

# Purpose: loads bcp file
# Returns: exist 2 if bcp fails
# Assumes: None
# Effects: None
#
def bcpFiles():

    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

    bcpCmd = '%s %s %s %s %s %s "|" "\\n" mgd' % \
            (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), 'MGI_Relationship', outputDir, bcpFile)
    print(bcpCmd)
    rc = os.system(bcpCmd)

    if rc != 0:
        closeFiles()
        sys.exit(2)
    # update mgi_relationship auto-sequence
    db.sql(''' select setval('mgi_relationship_seq', (select max(_Relationship_key) from MGI_Relationship)) ''', None)

    return 0

# end bcpFiles() -------------------------------------

#
# Main
#

print('%s' % mgi_utils.date())
print ('init()')
if init() != 0:
    exit(1, 'Error in  init \n' )

print('%s' % mgi_utils.date())
print('findRelationships()')
# determine qtl candidate genes and write to bcp
if findRelationships() != 0:
    exit(1, 'Error in  findRelationships \n' )

#print('%s' % mgi_utils.date())
#print('writeReports()')
## write out info for each bucket of relationships
#if writeReports() != 0:
#    exit(1, 'Error in  writeReports \n' )

print('%s' % mgi_utils.date())
print('doDeletes()')
# delete existing relationships
if doDeletes() != 0:
    exit(1, 'Error in  doDeletes \n' )

print('%s' % mgi_utils.date())
print('closeFiles()')
# close all output files
if closeFiles() != 0:
    exit(1, 'Error in  closeFiles \n' )

print('%s' % mgi_utils.date())
print('bcpFiles()')
# bcp the relationships
if bcpFiles() != 0:
    exit(1, 'Error in  bcpFiles \n' )

print('%s' % mgi_utils.date())
print('done')
sys.exit(0)
