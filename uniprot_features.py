'''
http://pymolwiki.org/index.php/uniprot_features

(c) 2010-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd

def uniprot_features(uniprot_id, selection='(all)', withss=0,
        prefix='feature_', quiet=1):
    '''
DESCRIPTION

    Fetch feature list from uniprot.org and create named selections.

    Requires residue numbering (resi) to match uniprot sequence!

ARGUMENTS

    uniprot_id = string: UniProtKB name or accession

    selection = string: atom selection {default: all}

    withss = 0/1: update secondary structure {default: 0}
    '''
    from lxml import etree
    from urllib import urlopen

    withss, quiet = int(withss), int(quiet)

    url = 'http://www.uniprot.org/uniprot/%s.xml' % uniprot_id
    if not quiet:
        print 'Downloading', url
    doc = etree.parse(urlopen(url))

    ns = 'http://uniprot.org/uniprot'
    NS = {'u': ns}
    if not quiet:
        print 'Parsing Features'
    features = doc.findall('{%s}entry/{%s}feature' % (ns, ns))

    if withss == 1:
        cmd.alter(selection, 'ss="L"')
        ssmap = {'helix': 'H', 'strand': 'S', 'turn': 'L'}

    count = 0
    for feature in features:
        type = feature.get('type')
        begin = feature.xpath('u:location/u:begin/@position', namespaces=NS)
        if len(begin):
            end = feature.xpath('u:location/u:end/@position', namespaces=NS)
            sel = '(%s) and resi %s-%s' % (selection, begin[0], end[0])
        else:
            position = feature.xpath('u:location/u:position/@position',
                    namespaces=NS)
            sel = '(%s) and resi %s' % (selection, position[0])
        if type in ['helix', 'strand', 'turn'] and withss < 2:
            if withss == 1:
                cmd.alter(sel, 'ss="%s"' % ssmap.get(type, 'L'))
        else:
            count += 1
            name = cmd.get_legal_name('%s%03d_%s' % (prefix, count,
                feature.get('description', '').replace('.', '')))
            groupname = cmd.get_legal_name('%s%s' % (prefix, type))
            cmd.select(name, sel)
            cmd.group(groupname, name, 'add')

    if not quiet:
        print 'Found %d feature records (without secondary structures)' % count

cmd.extend('uniprot_features', uniprot_features)

# vi:expandtab:smarttab
