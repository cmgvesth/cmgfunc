import sys
import MySQLdb
import subprocess
import os

mysql_settings = {
    "passwd":"<insetr password>",
    "user":"<insert user>",
    "host":"mysql",
    "port":3306,
    "db":"<insert db name>"
}

load_feat_script = "/home/projects/steve/ProtFun/load_feat3.pl"

conn = MySQLdb.connect(**mysql_settings)
conn.autocommit(False)
curs = conn.cursor()

should_run = True
while( should_run ):
    # row_count is either one or zero if finished
    try:
	curs.execute("LOCK TABLES feature_processing WRITE")
        row_count = curs.execute("SELECT feature_processing.protein_id, stage FROM feature_processing WHERE stage=0 LIMIT 1")
        if( row_count is 0 ):
            sys.exit(1)
        result = curs.fetchone();
        if( not result ):
            sys.exit(1)
        protein_id = result[0];
        stage       = result[1];
        result = curs.execute("UPDATE feature_processing SET stage=1 WHERE protein_id=%s AND stage=0", (protein_id,));
        conn.commit();
	curs.execute("UNLOCK TABLES")
        if( not result ):
            continue;
    except Exception, e:
        print "Error finding protein id for processing"
        print e.message
	continue;
    print "Processing protein: %i" % protein_id
    temp_file = "tempfeat_%i" % os.getpid();
    ret_val = subprocess.call([load_feat_script, "-project", "1600genomes", "-temp", temp_file, "-protein", "%i" % protein_id]);
    if( not ret_val ): # Success
        curs.execute("UPDATE feature_processing SET stage=2 WHERE protein_id=%s", (protein_id,));
        conn.commit()
    else: # Failure?
        curs.execute("UPDATE feature_processing SET stage=-1 WHERE protein_id=%s", (protein_id,));
        conn.commit()
