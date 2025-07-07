from pymol import cmd
from io import StringIO

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio import SeqIO
    from Bio.Seq import Seq
except ModuleNotFoundError:
    # Note that Bio package might be missing from Pymol 2 installation!
    print("Oops! Protparam: Biopython is missing!\n If you want to install it, run protparam_dependencies_install command")


@cmd.extend
def protparam(selection='enabled', bychain=0):
    '''
    DESCRIPTION:
    Given selection, calculates common protein properties, like Mw, pI, length and aminoacid content.
    By default, combines all chains of each object into the single sequence.
    
    USAGE:
    protparam selection, [bychain]

    DEPENDENCIES:
    biopython
    '''
    #TODO: add pretty output suitable for copy-pasting
    for entry in cmd.get_object_list(selection):
        sequence_obj = cmd.get_fastastr(f"({selection}) and {entry}")
        fasta_io = StringIO(sequence_obj)
        sequences = list(SeqIO.parse(fasta_io, "fasta"))
        sequences = [s.seq for s in sequences]
        if not bychain:
            #by default combine all chains into single sequence
            sequences = [Seq('').join(sequences)]
        for sequence in sequences:
            sequence = str(sequence).replace('?','').strip()
            analysis = ProteinAnalysis(sequence)
            counts_aa = analysis.count_amino_acids() #Dict is useful when only specific residues should be reported
            print(f"Protein name: {entry}")
            print(f"Sequence: {sequence}")
            print(f"\nProtein length: {analysis.length} aa")
            print(f"Molecular Weight: {analysis.molecular_weight():.1f} Da")
            print(f"Isoelectric point: {analysis.isoelectric_point():.2f}")
            print(f"Count of aminoacids: {counts_aa}\n\n")

@cmd.extend
def protparam_dependencies_install():
    import sys
    import subprocess
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", 'biopython'])
        print(f"Successfully installed biopython! Reload Protparam plugin or restart PyMOL.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to install biopython: {e}")

def test_protparam(capsys):
    cmd.reinitialize()
    cmd.fab("A// ACD B// EFG", "m1")
    cmd.fab("HIKL", "m2")
    cmd.alter("resn CYS", "resn='UNK'")
    protparam()
    captured = capsys.readouterr()
    assert "Protein name: m1" in captured.out
    assert "Protein name: m2" in captured.out
    assert "Sequence: ADEFG\n" in captured.out
    assert "Sequence: HIKL\n" in captured.out
    assert "Protein length: 2 aa" not in captured.out
    assert "Protein length: 3 aa" not in captured.out
    assert "Protein length: 4 aa" in captured.out
    assert "Protein length: 5 aa" in captured.out
    assert "Count of aminoacids: {'A': 1," in captured.out
    protparam(bychain=1)
    captured = capsys.readouterr()
    assert "Protein name: m1" in captured.out
    assert "Protein name: m2" in captured.out
    assert "Sequence: AD\n" in captured.out
    assert "Sequence: EFG\n" in captured.out
    assert "Protein length: 2 aa" in captured.out
    assert "Protein length: 3 aa" in captured.out
    assert "Protein length: 4 aa" in captured.out
    assert "Protein length: 5 aa" not in captured.out
    assert "Molecular Weight: 204.2 Da" in captured.out
    protparam("resn LYS")
    captured = capsys.readouterr()
    assert "Protein name: m1" not in captured.out
    assert "Protein name: m2" in captured.out
    assert "Protein length: 1 aa" in captured.out
    assert "Isoelectric point: 8.75" in captured.out
    protparam("resn TRP")
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""