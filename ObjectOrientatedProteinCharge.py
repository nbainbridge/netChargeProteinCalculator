from collections import Counter
import csv


class Residue:
    def __init__(self, name, count, pKa, charge):
        self.name = name
        self.count = count
        self.pKa = pKa
        self.charge = None

    def printout(self):
        """ Printout all values stored within the class """
        print(self.name, self.count, self.pKa, self.charge, )

    def totalcharge(self):
        """Return value of charge of 1 residue * number of that residues in the sequence"""
        return float(self.charge * self.count)

    def findcharge(self):
        """Determine the charge of specific group: HRK and NH are basic side chains, others are acidic """
        if self.name not in ["H", "R", "K", "NH"]:
            self.charge = -1 + 10 ** -pH / (10 ** -pH + 10 ** (-1 * float(self.pKa)))

        else:
            self.charge = 10 ** -pH / (10 ** -pH + 10 ** (-1 * float(self.pKa)))


pH = 0.606
totalCharge = 0
aaList = []
""" All charged side chains:"""
AA = ['H', 'R', 'K', 'D', 'E', 'C', 'Y']
seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

counter = Counter(seq)
val_Pka = open('E:\Konnerman Group\Protein Charge\SigmaAldrich_pKaValues', 'r')

for line in val_Pka:
    line = line.split('-')

    """pKa for charge of N terminal"""
    if line[0] == seq[0]:
        aaList.append(Residue('NH', 1, line[2], None))
    """pKa for charge of C terminal"""
    if line[0] == seq[len(seq) - 1]:
        aaList.append(Residue('COOH', 1, line[1], None))
    """pKa for charge of all charged R group residues"""
    if line[0] in AA and counter[line[0]] != 0:
        aaList.append(Residue(line[0], counter[line[0]], line[3].strip(), None))
val_Pka.close()


header = ['pH', 'Charge']
with open('charge.csv', 'a', encoding='UTF8', newline='') as f:
    f.truncate(0)
    writer = csv.writer(f)
    writer.writerow(['pH', 'Charge'])

    while pH <= 14:
        totalCharge = 0
        for val in aaList:
            val.findcharge()
            totalCharge += val.totalcharge()
        writer.writerow([round(pH, 2), totalCharge])
        print(round(pH, 2), ',', totalCharge)
        pH = pH + 0.1

        totalCharge = 0




