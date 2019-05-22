
fl2=open("gimme.vertebrate.v3.3.1.xt.factortable.txt","w")
fl2.write("motif\tfactor\n")

with open("gimme.vertebrate.v3.3.1.xt.motif2factors.txt") as m2f:
    for line in m2f:
        # print(line)
        if len(line.split()) >1:
            for i in line.split()[1].split(","):
                fl2.write(line.split()[0]+"\t"+i+"\n")



