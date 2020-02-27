'''
http://pymolwiki.org/index.php/uniprot_features

(c) 2010-2012 Thomas Holder, MPI for Developmental Biology
(c) 2012 Troels Linnet, SBiNLab Copenhagen University

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException


class resid_mapper(dict):
    '''
    DESCRIPTION

    Residue identifier mapping
    '''
    @classmethod
    def from_seq_sel(cls, sequence, selection):
        '''
        Constructor with sequence and PyMOL selection.

        Requires psico, biopython and needle.
        '''
        from psico import one_letter
        from psico.seqalign import needle_alignment, alignment_mapping

        NL, VL = [], []
        cmd.iterate('(%s) and guide' % (selection),
                    'NL.append(resn);VL.append(resv)', space=locals())
        seq = ''.join(one_letter.get(resn, 'X') for resn in NL)
        align = needle_alignment(sequence, seq)

        return cls((i + 1, VL[j]) for (i, j) in alignment_mapping(*align))

    def __call__(self, k):
        if isinstance(k, tuple):
            a, b = map(int, k)
            return self.search(a), self.search(b, -1)
        return self[int(k)]

    def search(self, resv, d=1):
        assert d in (1, -1)
        stop = max(self) if d > 0 else min(self)
        for resv in range(resv, stop + d, d):
            if resv in self:
                return self.get(resv)
        raise KeyError


def uniprot_features(uniprot_id, selection='(all)', withss=0,
                     prefix='feature_', sm=None, quiet=1):
    '''
DESCRIPTION

    Fetch feature list from uniprot.org and create named selections.

    Requires residue numbering (resi) to match uniprot sequence!

ARGUMENTS

    uniprot_id = string: UniProtKB name or accession

    selection = string: atom selection {default: all}

    withss = 0/1: update secondary structure {default: 0}
    '''
    import xml.etree.ElementTree as etree
    try:
        from urllib import urlopen
    except ImportError:
        from urllib.request import urlopen

    withss, quiet = int(withss), int(quiet)

    url = 'http://www.uniprot.org/uniprot/%s.xml' % uniprot_id
    if not quiet:
        print('Downloading', url)
    doc = etree.parse(urlopen(url))

    ns = 'http://uniprot.org/uniprot'
    NS = {'u': ns}
    if not quiet:
        print('Parsing Features')
    features = doc.findall('{%s}entry/{%s}feature' % (ns, ns))
    sequence = doc.findtext('{%s}entry/{%s}sequence' % (ns, ns))
    sequence = ''.join(sequence.split())

    try:
        if sm is None:
            sm = resid_mapper.from_seq_sel(sequence, selection)
    except:
        print(' Warning: sequence mapping failed')
        sm = lambda x: x

    if withss == 1:
        cmd.alter(selection, 'ss="L"')
        ssmap = {'helix': 'H', 'strand': 'S', 'turn': 'L'}

    norange_types = ['disulfide bond']

    count = 0
    for feature in features:
        type_ = feature.get('type')
        begin = feature.find('{%s}location/{%s}begin' % (ns, ns))
        if begin is not None:
            end = feature.find('{%s}location/{%s}end' % (ns, ns))
            values = begin.get('position'), end.get('position')
            if type_ in norange_types:
                try:
                    values = sm(values[0]), sm(values[1])
                    sel = '(' + selection + ') and resi %s+%s' % values
                except KeyError:
                    print(' Warning: could not map', values)
                    sel = 'none'
            else:
                try:
                    values = sm(values)
                    sel = '(' + selection + ') and resi %s-%s' % values
                except KeyError:
                    print(' Warning: could not map', values)
                    sel = 'none'
        else:
            position = feature.find('{%s}location/{%s}position' % (ns, ns))
            try:
                value = sm(position.get('position'))
                sel = '(%s) and resi %s' % (selection, value)
            except KeyError:
                sel = 'none'
        if type_ in ['helix', 'strand', 'turn'] and withss < 2:
            if withss == 1:
                cmd.alter(sel, 'ss="%s"' % ssmap.get(type_, 'L'))
        else:
            count += 1
            name = cmd.get_legal_name('%s%03d_%s' % (prefix, count,
                                                     feature.get('description', '').replace('.', '')))
            groupname = cmd.get_legal_name('%s%s' % (prefix, type_))
            cmd.select(name, sel)
            cmd.group(groupname, name, 'add')

    if not quiet:
        print('Found %d feature records (without secondary structures)' % count)


def uniprot_auto(pdb_id, selection='', withss=0, quiet=1):
    '''
DESCRIPTION

    Like "uniprot_features" but with automatic fetching of UniProtKB accession
    and sequence mapping for given pdb_id from http://www.bioinf.org.uk/pdbsws/

ARGUMENTS

    pdb_id = string: PDB accession ID

    selection = string: atom selection {default: <pdb_id>, will be fetched if
    no such object is loaded}

    withss = 0/1: update secondary structure {default: 0}
    '''
    try:
        from urllib import urlopen
    except ImportError:
        from urllib.request import urlopen

    if len(pdb_id) != 4 or not pdb_id[0].isdigit():
        raise CmdException('invalid pdb_id: ' + pdb_id)

    if not selection:
        selection = pdb_id
        if pdb_id not in cmd.get_names('all'):
            cmd.fetch(pdb_id)

    sele_chains = cmd.get_chains(selection)
    mappings = {}

    pdb_id = pdb_id.lower()
    url = 'http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&all=yes&id=' + pdb_id

    try:
        for line in urlopen(url):
            if not isinstance(line, str):
                line = line.decode('utf-8')
            if not line.startswith(pdb_id):
                continue

            chain = line[5]
            resno = line[20:25].strip()
            acc = line[27:36].strip()
            number = line[40:50].strip()

            if not acc or not number:
                continue

            if chain not in mappings:
                mappings[chain] = acc, resid_mapper()

            if mappings[chain][0] != acc and chain in sele_chains:
                raise ValueError('multiple accessions per chain not supported')

            mappings[chain][1][int(number)] = resno
    except Exception as e:
        raise CmdException(str(e))

    for chain, (acc, sm) in mappings.items():
        uniprot_features(acc, '(%s) and chain %s' % (selection, chain),
                         withss, 'feature_' + chain + '_', sm, quiet)

cmd.extend('uniprot_features', uniprot_features)
cmd.extend('uniprot_auto', uniprot_auto)

# vi:expandtab:smarttab
