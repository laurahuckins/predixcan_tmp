#!/usr/bin/env python

import sqlite3
import os

SOURCE_DIR =  'Where_resfiles_are/'

conn = sqlite3.connect('DLPFC.tmp.db')

c = conn.cursor()

#c.execute('''CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)''')

c.execute("CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
c.execute("CREATE INDEX weights_rsid ON weights (rsid)")
c.execute("CREATE INDEX weights_gene ON weights (gene)")
c.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
c.execute("CREATE TABLE extra (ENS TEXT, gene TEXT, R2 DOUBLE,  `n.snps` INTEGER)")
c.execute("CREATE INDEX extra_gene ON extra (gene)")

#Gene	SNP	eff.allele	ref.allele	beta

def weights_per_file(filename):
    file=open(filename)
    for line in file:
        raw=line.split()
        if raw[0] == "genename":
            print "Ignoring Header Row!"
        else:
            c.execute("INSERT INTO weights VALUES (?,?,?,?,?, NULL, NULL, NULL)",
                      (raw[1], raw[0], raw[4], raw[3], raw[2]))

def extras_per_file(filename):
    file=open(filename)
    for line in file:
        raw=line.split()
        if raw[0] == "ENSID":
            print "Ignoring Header Row!"
        else:
            c.execute("INSERT INTO extra VALUES(?, ?, ?, ?)", (raw[0], raw[1], raw[2], raw[3]))



# Get all the file lists that we need

source_files = os.popen("ls NEW_EU_ONLY/noSVA_ancestry_1mb/betas* | grep -v txt ")
smartlist=source_files.read().split()

#extras_files = os.popen("ls /sc/orga/scratch/huckil01/predictors_CMC_Best/*extra*") # Honestly you don't need this step IMO
#extras_smartlist=extras_files.read().split()

for i in xrange(0, len(smartlist)):
    weights_per_file(smartlist[i])

conn.commit()

#for i in xrange(0, len(extras_smartlist)):
#    extras_per_file(extras_smartlist[i])

#conn.commit()
conn.close()
