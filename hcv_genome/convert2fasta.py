info = open('HCV.info','r').read().split('\n')
fasta = open('HCV.fasta','w')
fasta.write('>hcv\n')
genome = ''
for line in info:
   try:
      num = int(line.split()[0])
      fasta.write(''.join(line.split()[1:]).upper()+'\n')
      genome = genome + ''.join(line.split()[1:]).upper()
   except:
      pass

print(len(genome))

#fasta.write(info[0]+'\n'+genome)

