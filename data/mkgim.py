fl2=open("gimme.vertebrate.v3.3.1.xt.motif2factors.txt","w")
with open("gimme.vertebrate.v3.3.xt.motif2factors.txt","r") as fl1:
    for i in fl1:
        a=i.split()
        if len(a)==2:
            fl2.write('{0}\t{1}\n'.format(a[0],a[1].upper()))
        else:
            fl2.write(i)
