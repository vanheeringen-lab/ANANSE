
fl2=open("gimme.vertebrate.v5.1.motif2factors.txt","w")

motifs={}

with open("gimme.vertebrate.v5.0.motif2factors.txt") as m2f:
	next(m2f)
	for line in m2f:
		a=line.split()
		if a[0] in motifs:
			motifs[a[0]]=motifs[a[0]]+','+a[1].upper()
		else:
			motifs[a[0]]=a[1].upper()
for motif in motifs:	
	fl2.write(motif+"\t")
	m=set()
	for i in motifs[motif].split(","):
		m.add(i)
	ll=""
	for j in m:
		ll+=","+j
	fl2.write(ll[1:])
	fl2.write("\n")


fl1=open("gimme.vertebrate.v5.1.pfm","w")

with open ("gimme.vertebrate.v5.0.pfm") as pfm:
	for i in pfm:
		if i.startswith("#"):
			pass
		else:
			fl1.write(i)
			