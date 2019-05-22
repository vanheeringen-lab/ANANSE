

fa2name={}
fa2=open("/home/qxu/projects/regulatoryNetwork/run20180719/scripts/data/gene2name.txt","r")
for i in fa2:
    a=i.split()
    if a[0].startswith("gene"):
        fa2name[a[1]]=a[0]

fl1=open(r"gimme.vertebrate.v3.3.xt.motif2factors.txt","r")
fl2=open(r"gimme.vertebrate.v3.3.1.xt.motif2factors.txt","wb")

for i in fl1:
    a=i.split()
    fl2.write(a[0]+"\t")
    if len(a)==2:
        b=a[1].split(",")
        for j in b:
            if j in fa2name:
                fl2.write(fa2name[j]+"\t")
    fl2.write("\n")
