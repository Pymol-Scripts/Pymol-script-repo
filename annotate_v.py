"""
annotate_v
A simple-to-use Python script annotating VH or VL sequences of an antibody.

Description of usage, see: https://pymolwiki.org/index.php/Annotate_v

Author: Xin Yu, 2020

"""

import requests
from pymol import cmd

class annotate():
    """
    class `annotate`.

    Initiator `__init__` has 2 parameters:

    :param aaseq: STRING: A SINGLE-LETTER, amino acid sequence corresponding to the complete VH or VL chain. Both uppercase and lowercase are accepted.

    :param scheme: STRING: "kabat", "chothia", "contact", or "imgt". Must be in LOWERCASE

    Class has 3 methods. `retrieve()`: retrieves numbered seqs from Abnum website, then sends it to method `analyze` to determine the FR and CDR regions, and to `output() ` to print the result and return a list of 2 dictionaries, the first of which contains to region:seq pairs, the second of which contains number:residue pairs.

    """

    def __init__(self, aaseq, scheme):

        self.aaseq = aaseq
        self.scheme = scheme

    def __repr__(self):
        return "Annotation of VH or VL sequence using Kabat, Chothia, Contact, or IMGT scheme"

    def output(self, regionlst):
        """
        Prints the FR and CDR regions and their corresponding seq. It returns a `list` of 2 `dict`.

        :param chain: STRING, either "H" or "L" in uppercase
        :param lst:  LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :param regionlst: LIST, a list of peptides, each corresponds to a FR or CDR region
        :return: LIST, a list of 2 `dict`, The first dict consists of region: seq pairs. The second dict consists of number:residue pairs.

        """

        self.regionlst = regionlst

        self.regiondict= {}


        if self.scheme == "kabat":
            print("Annotation in Kabat scheme:")
        elif self.scheme == "chothia":
            print("Annotation in Chothia scheme:")
        elif self.scheme == "contact":
            print("Annotation in Contact scheme:")
        else:
            print("Annotation in IMGT scheme:")

        if self.chain == "L":
            print("L-FR1:  ", self.regionlst[0])
            print("L-CDR1: ", self.regionlst[1])
            print("L-FR2:  ", self.regionlst[2])
            print("L-CDR2: ", self.regionlst[3])
            print("L-FR3:  ", self.regionlst[4])
            print("L-CDR3: ", self.regionlst[5])
            print("L-FR4:  ", self.regionlst[6])

            for region, seq in zip(["L-FR1", "L-CDR1", "L-FR2", "L-CDR2", "L-FR3", "L-CDR3", "L-FR4"], self.regionlst):
                self.regiondict[region] = seq

            return self.regiondict

        else:
            print("H-FR1:  ", self.regionlst[0])
            print("H-CDR1: ", self.regionlst[1])
            print("H-FR2:  ", self.regionlst[2])
            print("H-CDR2: ", self.regionlst[3])
            print("H-FR3:  ", self.regionlst[4])
            print("H-CDR3: ", self.regionlst[5])
            print("H-FR4:  ", self.regionlst[6])

            for region, seq in zip(["H-FR1", "H-CDR1", "H-FR2", "H-CDR2", "H-FR3", "H-CDR3", "H-FR4"], self.regionlst):
                self.regiondict[region] = seq

            return self.regiondict

    def analyze(self, chain, lst):
        """
        Define CDR and FR regions based on the numbered sequence returned from website

        :param chain: STRING, "H" or "L" in uppercase
        :param lst: LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :return: LIST, a list of strings, where each string is a peptide corresponding to the a region, in the order of: FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4

        :raises: `ValueError` if any of the FR or CDR region is missing

        """

        self.chain = chain
        self.lst = lst
        if self.chain == "L":
            self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4 = ["" for i in
                                                                                                     range(0, 7)]

            try:
                if self.scheme in ["kabat", "chothia"]:
                    self.L_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("L24"), 2)])
                    self.L_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L24"), self.lst.index("L35"), 2)])
                    self.L_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L35"), self.lst.index("L50"), 2)])
                    self.L_CDR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L50"), self.lst.index("L57"), 2)])
                    self.L_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L57"), self.lst.index("L89"), 2)])
                    self.L_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("L98"), len(self.lst), 2)])

                elif self.scheme == "contact":
                    self.L_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("L30"), 2)])
                    self.L_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L30"), self.lst.index("L37"), 2)])
                    self.L_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L37"), self.lst.index("L46"), 2)])
                    self.L_CDR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L46"), self.lst.index("L56"), 2)])
                    self.L_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L56"), self.lst.index("L89"), 2)])
                    self.L_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L89"), self.lst.index("L97"), 2)])
                    self.L_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("L97"), len(self.lst), 2)])

                else:  # IMGT scheme
                    self.L_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("L27"), 2)])
                    self.L_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L27"), self.lst.index("L33"), 2)])
                    self.L_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L33"), self.lst.index("L50"), 2)])
                    self.L_CDR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L50"), self.lst.index("L52"), 2)])
                    self.L_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L52"), self.lst.index("L89"), 2)])
                    self.L_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("L98"), len(self.lst), 2)])

                return [self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4]

            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured")
        else:
            self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4 = ["" for i in
                                                                                                     range(0, 7)]
            try:
                if self.scheme == "kabat":
                    self.H_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("H31"), 2)])
                    self.H_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H31"), self.lst.index("H36"), 2)])
                    self.H_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H36"), self.lst.index("H50"), 2)])
                    self.H_CDR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H50"), self.lst.index("H66"), 2)])
                    self.H_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H66"), self.lst.index("H95"), 2)])
                    self.H_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("H103"), len(self.lst), 2)])

                elif self.scheme == "chothia":
                    self.H_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("H26"), 2)])
                    self.H_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H26"), self.lst.index("H33"), 2)])
                    self.H_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H33"), self.lst.index("H52"), 2)])
                    self.H_CDR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H52"), self.lst.index("H57"), 2)])
                    self.H_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H57"), self.lst.index("H95"), 2)])
                    self.H_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("H103"), len(self.lst), 2)])

                elif self.scheme == "contact":
                    self.H_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("H30"), 2)])
                    self.H_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H30"), self.lst.index("H36"), 2)])
                    self.H_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H36"), self.lst.index("H47"), 2)])
                    self.H_CDR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H47"), self.lst.index("H59"), 2)])
                    self.H_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H59"), self.lst.index("H93"), 2)])
                    self.H_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H93"), self.lst.index("H102"), 2)])
                    self.H_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("H102"), len(self.lst), 2)])

                else:  # IMGT scheme
                    self.H_FR1 = "".join([self.lst[i + 1] for i in range(0, self.lst.index("H26"), 2)])
                    self.H_CDR1 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H26"), self.lst.index("H34"), 2)])
                    self.H_FR2 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H34"), self.lst.index("H51"), 2)])
                    self.H_CDR2 = "".join([self.lst[i + 1] for i in range(self.lst.index("H51"), self.lst.index("H58"),
                                                                          2)])  # 51>57 (instead of 56)
                    self.H_FR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H58"), self.lst.index("H93"), 2)])
                    self.H_CDR3 = "".join(
                        [self.lst[i + 1] for i in range(self.lst.index("H93"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i + 1] for i in range(self.lst.index("H103"), len(self.lst), 2)])

                return [self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4]

            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured in the `analyze()` method")

    def retrieve(self):
        """
        Retrieve numbered residues from Abnum website

        :return: returns same object from the `output()` method.

        :raises: `ValueError` if input scheme is not among "kabat", "chothia", "contact", and "imgt"

        """

        self.url = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi"

        try:
            if self.scheme not in ["kabat", "chothia", "contact", "imgt"]:
                raise Exception

        except ValueError:
            print("Incorrect scheme mode. Must be one of the following (lowercase): kabat, chothia, contact, imgt")

        else:
            if self.scheme == "kabat":
                self.sche = "-k"
            else:
                self.sche = "-c"

        try:
            self.d = {"plain": 1, "scheme": self.sche, "aaseq": self.aaseq}
            self.myPage = requests.get(self.url, params=self.d)
            self.text = self.myPage.text
            self.lst = self.text.split()

            if len(self.lst) > 1:
                self.chain = self.lst[0][0]
                self.regiondict = self.output(self.analyze(self.chain, self.lst))
                return self.regiondict
            else:
                print("No annotation retrieved. Did you enter the complete VH or VL sequence?")
        except:
            print("An error occured in the `retrieve()` method")


def annotate_v(selection, scheme):
    aaseq="".join(cmd.get_fastastr(selection).split("\n")[1:])
    obj=annotate(aaseq, scheme)
    result=obj.retrieve()
    for i in result.keys():
        cmd.select(i, "pepseq %s" % result[i])


cmd.extend("annotate_v", annotate_v)

