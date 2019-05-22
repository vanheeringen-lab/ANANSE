#!/usr/bin/env python
import sys
from pita.io import read_gff_transcripts
from BCBio import GFF
from pita.util import get_overlapping_models

if not len(sys.argv) == 2:
    sys.stderr.write("Usage: gff3tobed12 <gff3> > <bed12>\n")
    sys.exit(1)

infile = sys.argv[1]
#infile = "test.gff"
fobj = open(infile)

def _gff_type_iterator(feature, ftypes):
    if not hasattr(feature, "sub_features") or len(feature.sub_features) == 0: #feature.type in ftypes:
        #print feature.type, len(feature.sub_features)
        yield feature
    else:
        for feature in feature.sub_features:
            for f in _gff_type_iterator(feature, ftypes):
                yield f
bedline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}"
smap = {"1":"+",1:"+","-1":"-",-1:"-", None:"+"}
for rec in GFF.parse(fobj):
    chrom = rec.id
    for feature in rec.features:

        gstart = int(feature.location.start.position)
        gend = int(feature.location.end.position)    
        gstrand = smap[feature.strand]
        gname = feature.qualifiers["Name"][0]
     
        allexons=[]

        bla = feature.sub_features
        while bla != [] and bla[0].type not in ["mRNA", "ncRNA", "transcript"]:
            tmp = []
            for f in bla:
                tmp += f.sub_features
            bla = tmp
        #if bla == []:
        #    bla = [feature]
        #print bla
        #for f in _gff_type_iterator(feature, []):
        #    print f
        for gene in bla: 
            exons = []
            for f in _gff_type_iterator(gene, []):
                start = int(f.location.start.position)
                end = int(f.location.end.position)    
                strand = smap[f.strand]
                name = f.qualifiers['Parent'][0]
                typ = f.type
                exons.append([start, end, name, strand, typ])
            exons = sorted(exons, lambda x,y: cmp(x[0], y[0]))
            exons = [e for e in exons if e[-1] == "exon"]
            
            if len(exons) == 0:
                sys.stderr.write(str(gene))
                sys.stderr.write("\n")
                continue
            cds_start = exons[0][0] 
            cds_end = exons[-1][1]
            if exons[0][4].upper().find("UTR") > -1:
                c = 0
                while exons[c + 1][4].upper().find("UTR") > -1:
                    c += 1
                cds_start = exons[c][1] 
                exons[c + 1][0] = exons[c][0]
                del exons[c]
            if exons[-1][4].upper().find("UTR") > -1:
                c = -1
                while exons[c - 1][4].upper().find("UTR") > -1:
                    c -= 1
                cds_end = exons[c][0]
                exons[c - 1][1] = exons[c][1]
                del exons[c]
           
            start = exons[0][0] 
            end = exons[-1][1]
            starts = [e[0] - start for e in exons]
            sizes = [e[1] - e[0] for e in exons]
            
            for eee in exons:
                allexons.append(eee)
        #allexons=list(set(allexons))
        nallexons=[]
        for j in allexons:
            ex=[j[0],j[1],j[3]]
            if ex not in nallexons:
                nallexons.append(ex)

        nstarts=[]
        nsizes=[]       
        for ee in nallexons:
            nstart = ee[0] 
            nend = ee[1]
            nstarts.append(ee[0] - gstart)
            nsizes.append(ee[1] - ee[0])

        # print nstarts
        # print nsizes
        if len(nallexons)>0:
            print(bedline.format(
                                 chrom,
                                 gstart,
                                 gend,
                                 gname,
                                 0,
                                 gstrand,
                                 gstart,
                                 gend,
                                 "0,0,0",
                                 len(nallexons),
                                 ",".join([str(x) for x in nsizes] + [""]),
                                 ",".join([str(x) for x in nstarts] + [""]),
                                 ))                                 
