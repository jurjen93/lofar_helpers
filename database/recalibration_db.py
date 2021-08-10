"""
This code is based on https://github.com/mhardcastle/lotss-query/blob/master/surveys_db.py
"""

from __future__ import print_function
from builtins import object
import sshtunnel
import socket
import os
from time import sleep

try:
    import MySQLdb as mdb
    import MySQLdb.cursors as mdbcursors
except ImportError:
    import pymysql as mdb
    import pymysql.cursors as mdbcursors


class SurveysDB(object):
    '''Database for recalibration'''

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        # if 'closed' doesn't exist, then we are most likely being called through __del__ due to a failure in the init call. So skip the rest.
        if hasattr(self, 'closed'):
            if not self.closed:
                if not self.readonly:
                    self.cur.execute('unlock tables')
                self.con.close()
                if self.usetunnel:
                    self.tunnel.stop()
                self.closed = True  # prevent del from trying again

    def __init__(self, readonly=False, verbose=False, survey=None):

        if survey is None:
            survey = 'hba'  # preserve old default behaviour
        # get the config file -- this must exist
        home = os.getenv("HOME")
        mysql_host = os.getenv('DDF_PIPELINE_MYSQLHOST')
        if not mysql_host:
            mysql_host = 'lofar-server.data'
        if verbose:
            print('MySQL host is', mysql_host)
        cfg = [l.rstrip() for l in open(home + '/.surveys').readlines()]
        self.password = cfg[0]
        try:
            self.ssh_user = cfg[1]
        except:
            self.ssh_user = None

        try:
            self.ssh_key = cfg[2]
        except:
            self.ssh_key = "id_rsa"

        self.readonly = readonly
        self.verbose = verbose
        self.survey = survey
        if self.survey == 'hba':
            self.database = 'surveys'
        elif self.survey == 'lba':
            self.database = 'lba'
        else:
            raise NotImplementedError('Survey "%s" not known' % self.survey)

        # set up an ssh tunnel if not running locally
        self.usetunnel = False
        self.hostname = socket.gethostname()
        if self.hostname == 'lofar-server':
            if verbose:
                print('Using direct connection to localhost')
            self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, self.database,
                                   cursorclass=mdbcursors.DictCursor)
        else:
            try:
                socket.gethostbyname(mysql_host)
            except socket.gaierror:
                if verbose:
                    print('Cannot find host', mysql_host, 'will use tunnel')
                self.usetunnel = True

            if self.usetunnel:
                self.tunnel = sshtunnel.SSHTunnelForwarder('lofar.herts.ac.uk',
                                                           ssh_username=self.ssh_user,
                                                           ssh_pkey=home + '/.ssh/%s' % self.ssh_key,
                                                           remote_bind_address=('127.0.0.1', 3306),
                                                           local_bind_address=('127.0.0.1',))

                self.tunnel.start()
                localport = self.tunnel.local_bind_port
                self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, self.database, port=localport,
                                       cursorclass=mdbcursors.DictCursor)
            else:
                connected = False
                retry = 0
                while not connected and retry < 10:
                    try:
                        self.con = mdb.connect(mysql_host, 'survey_user', self.password, self.database,
                                               cursorclass=mdbcursors.DictCursor)
                        connected = True
                    except mdb.OperationalError as e:
                        time = 60
                        print('Database temporary error! Sleep %i seconds to retry\n' % time, e)
                        retry += 1
                        sleep(time)
                if not connected:
                    raise RuntimeError("Cannot connect to database server after repeated retry")
        self.cur = self.con.cursor()

        # get the tables list for locking
        self.cur.execute('show tables')
        result = self.cur.fetchall()
        self.tables = [d.items()[0][1] for d in result]

        if self.readonly:
            pass
            # can't use this feature on lofar's version of MariaDB
            # self.cur.execute('set session transaction read only')
        else:
            command = 'lock table '
            for table in self.tables:
                command += table + ' write, '
            command = command[:-2]
            self.cur.execute(command)
        self.closed = False

    def execute(self, *args):
        if self.verbose:
            print(args)
        self.cur.execute(*args)

    def db_get(self, table, id):
        table = self.check_table(table)
        self.execute('select * from ' + table + ' where id=%s', (id,))
        result = self.cur.fetchall()
        if len(result) == 0:
            return None
        else:
            return result[0]

    def check_table(self, table):
        if table not in self.tables:
            table += 's'
            if table not in self.tables:
                raise RuntimeError('Unknown table %s requested' % table)
        return table
