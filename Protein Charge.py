from collections import Counter


class Residue:
    def __init__(self,name,charge,count):
        self.name = name
        self.charge = charge
        self.count = count

    def printout(self):
        print(self.name, self.charge, self.count)

    def totalcharge(self):
        return float(self.charge * self.count)


class Termini:
    def __init__(self, name,pKa, charge, nc):
        self.name = name
        self.pKa = pKa
        self.charge = charge
        self.nc = nc

    def printout(self):
        print(self.name, self.pKa, self.charge, self.nc)

    def totalcharge(self):
        return float(self.charge)


resiSeq = input("Enter Sequence:  ")
pH = float(input("Enter pH:   "))
AA = ['H', 'R', 'K', 'D', 'E', 'C', 'Y']
pka_dict = {}
charge_dict = {}
summation_dict = []
count_dict = {}
# charge_dict = {}
counter = Counter(resiSeq)
for resi in AA:
    count_dict[resi] = counter[resi]

termN = resiSeq[0]
termC = resiSeq[len(resiSeq) - 1]
termNpKa = 0
termCpKa = 0
val_Pka = open('AA_pKa', 'r')
for line in val_Pka:
    line = line.split('-')

    if line[0] == termN:
        # pKa of N terminal Residue
        termNpKa = line[2]
    if line[0] == termC:
        # pKa of C terminal Residue
        termCpKa = line[1]
    if line[0] in AA:
        # pKa of charged Amino acids
        pka_dict[line[0]] = line[3].strip()
val_Pka.close()
nTerminalCharge = 10 ** -pH / (10 ** -pH + 10 ** (-1 * float(termNpKa)))
cTerminalCharge = -1 + 10**-pH / (10**-pH + 10**(-1 * float(termCpKa)))
for value in pka_dict:
    if value not in ["H", "R", "K"]:
        charge_dict[value] = -1 + 10**-pH / (10**-pH + 10**(-1 * float(pka_dict[value])))
    else:
        charge_dict[value] = 10 ** -pH / (10 ** -pH + 10 ** (-1 * float(pka_dict[value])))

nTerminus = Termini(termN, termNpKa, nTerminalCharge, 'n')
cTerminus = Termini(termC, termCpKa, cTerminalCharge, 'c')
nTerminus.printout()
cTerminus.printout()

aaList = []
for i in AA:
    aaList.append(Residue(i, charge_dict[i], count_dict[i]))

overallCharge = nTerminus.totalcharge() + cTerminus.totalcharge()
for val in aaList:
    overallCharge = val.totalcharge() + overallCharge
    val.printout()
print(overallCharge)

##Don't forget about N and C terminus!!