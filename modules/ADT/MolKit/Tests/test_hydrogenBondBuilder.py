import unittest
from string import split
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder



class HydrogenBondBuilderConstructorTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an HydrogenBuilder
        """
        hb_builder = HydrogenBondBuilder()
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)


    def test_constructorOptions(self):
        """
         test possible constructor options
            options = distCutoff, distCutoff2, d2min, d2max,
                      d2min, d3max, a2min, a2max, a3min, a3max,
                      donorTypes, acceptorTypes, distOnly
        """
        val = 1.0
        hb_builder = HydrogenBondBuilder(distCutoff=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['distCutoff'],val) 
        hb_builder = HydrogenBondBuilder(distCutoff2=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['distCutoff2'],val) 

        val = 85
        hb_builder = HydrogenBondBuilder(d2max=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['d2max'],val) 
        hb_builder = HydrogenBondBuilder(d2min=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['d2min'],val) 
        hb_builder = HydrogenBondBuilder(a2min=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['a2min'],val) 
        hb_builder = HydrogenBondBuilder(a2max=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['a2max'],val) 

        hb_builder = HydrogenBondBuilder(d3min=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['d3min'],val) 
        hb_builder = HydrogenBondBuilder(d3max=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['d3max'],val) 
        hb_builder = HydrogenBondBuilder(a3min=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['a3min'],val) 
        hb_builder = HydrogenBondBuilder(a3max=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['a3max'],val) 

        val = ['Nam','Ng+']
        hb_builder = HydrogenBondBuilder(donorTypes=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['donorTypes'],val) 
        val = ['S3','03']
        hb_builder = HydrogenBondBuilder(acceptorTypes=val)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['acceptorTypes'],val) 

        val = True
        hb_builder = HydrogenBondBuilder(distOnly=True)
        self.assertEquals(hb_builder.__class__, HydrogenBondBuilder)
        self.assertEquals(hb_builder.paramDict['distOnly'],val) 



class HydrogenBondBuilderTests(unittest.TestCase):

    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/1crn_Hs.pdbq')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    def test_buildHbonds_babelTypes(self):
        """
        test assigning babel_types 
        """
        hb_builder = HydrogenBondBuilder()
        hb_builder.check_babel_types(self.mol.allAtoms)
        d = {}
        for a in self.mol.allAtoms:
            d[a.babel_type] = 1
        all_types= d.keys()
        all_types.sort()
        #print "all_types=", all_types
        #note: these are polar only hydrogens so no type 'HC'
        canned_list = ['C+', 'C2', 'C3', 'Cac', 'H',  'N3+', 
                        'Nam', 'Ng+', 'O-', 'O2', 'O3', 'S3']
        for k, v in zip(all_types, canned_list):
            self.assertEqual(k,v)


    def test_buildD(self):
        """
        test buildD: donor3Ats, acceptor3Ats, acceptor2Ats, donor2Ats
        """
        hb_builder = HydrogenBondBuilder()
        hb_builder.check_babel_types(self.mol.allAtoms)
        result = hb_builder.buildD(self.mol.allAtoms, hb_builder.paramDict)
        #result_keys = ['donor3Ats', 'acceptor3Ats', 'acceptor2Ats', 'donor2Ats']
        result_keys = ['donor3Ats', 'acceptor3Ats', 'hAts', 'acceptorAts', 'acceptor2Ats', 'donor2Ats']
        d = {}
        d['donor3Ats'] = ['N', 'OG1', 'OG1', 'SG', 'SG', 'OG', 'OG', 'SG', 'OG1', 'SG', 'OG1', 'OH', 'OG1', 'SG', 'OG1', 'SG', 'OH']
        d['acceptor3Ats'] = ['OG1', 'OG1', 'SG', 'SG', 'OG', 'OG', 'SG', 'OG1', 'SG', 'OG1', 'OH', 'OG1', 'SG', 'OG1', 'SG', 'OH'] 
        d['hAts'] = ['HN1', 'HN2', 'HN3', 'HG1', 'HN', 'HG1', 'HN', 'HN', 'HN', 'HG', 'HN', 'HN', 'HN', 'HN', 'HE', 'HH11', 'HH12', 'HH21', 'HH22', 'HN', 'HG', 'HN', 'HD21', 'HD22', 'HN', 'HN', 'HD21', 'HD22', 'HN', 'HN', 'HN', 'HE', 'HH11', 'HH12', 'HH21', 'HH22', 'HN', 'HN', 'HN', 'HG1', 'HN', 'HN', 'HN', 'HN', 'HN', 'HN', 'HG1', 'HN', 'HH', 'HN', 'HG1', 'HN', 'HN', 'HN', 'HN', 'HN', 'HN', 'HN', 'HN', 'HG1', 'HN', 'HN', 'HN', 'HN', 'HH', 'HN', 'HN', 'HD21', 'HD22']
        d['donor2Ats'] = ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'NE', 'NH1', 'NH2', 'N', 'N', 'ND2', 'N', 'N', 'ND2', 'N', 'N', 'N', 'NE', 'NH1', 'NH2', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'ND2']
        d['acceptor2Ats'] = ['O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'OD1', 'O', 'O', 'OD1', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'OE1', 'OE2', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'OD1', 'OD2', 'O', 'O', 'O', 'OD1', 'OXT']
        #acceptorAts = acceptor2Ats + acceptor3Ats
        d['acceptorAts'] = ['O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'OD1', 'O', 'O', 'OD1', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'OE1', 'OE2', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'OD1', 'OD2', 'O', 'O', 'O', 'OD1', 'OXT', 'OG1', 'OG1', 'SG', 'SG', 'OG', 'OG', 'SG', 'OG1', 'SG', 'OG1', 'OH', 'OG1', 'SG', 'OG1', 'SG', 'OH'] 
        for k in result_keys:
            #print 'testing ', k, ': '
            if result[k] is not None and len(result[k]):
                self.assertEqual(result[k].name, d[k])
            else:
                self.assertEqual(result[k], d[k])



    def test_process_part1(self):
        """
        check hAts and acceptorAts from buildD
        """
        hb_builder = HydrogenBondBuilder()
        hb_builder.check_babel_types(self.mol.allAtoms)
        result = hb_builder.buildD(self.mol.allAtoms, hb_builder.paramDict)
        hAts = result['hAts']
        self.assertEqual(len(hAts), 69)
        acceptorAts = result['acceptorAts']
        self.assertEqual(len(acceptorAts), 70)
        #hb_builder.process(result, result, hb_builder.paramDict)
        

    def test_process_part2(self):
        """
        check results of calls to removeNeighbors, filterBasedOnAngs,
        removeBadAts...
        """
        hb_builder = HydrogenBondBuilder()
        hb_builder.check_babel_types(self.mol.allAtoms)
        result = hb_builder.buildD(self.mol.allAtoms, hb_builder.paramDict)
        hAts = result['hAts']
        acceptorAts = result['acceptorAts']
        dist = hb_builder.paramDict['distCutoff']
        atDict = hb_builder.distSelector.select(hAts, acceptorAts, dist)
        self.assertEqual(len(atDict), 40)
        atDict = hb_builder.removeNeighbors(atDict)
        self.assertEqual(len(atDict), 30)

        donor2Ats = result['donor2Ats']
        donor3Ats = result['donor3Ats']
        acceptor2Ats = result['acceptor2Ats']
        acceptor3Ats = result['acceptor3Ats']
        badAtDict = hb_builder.filterBasedOnAngs(atDict, donor2Ats, donor3Ats,
                        acceptor2Ats, acceptor3Ats, hb_builder.paramDict)
        self.assertEqual(len(badAtDict), 30)
        atDict = hb_builder.removeBadAts(atDict, badAtDict)
        self.assertEqual(len(atDict), 28)
        #hb_builder.process(result, result, hb_builder.paramDict)
        #self.assertEqual(len(self.mol.allAtoms.get(lambda x: hasattr(x, 'hbonds')), 28)
        


    def test_buildHbonds(self):
        """
        check results of single call to build
        """
        hb_builder = HydrogenBondBuilder()
        hb_builder.build(self.mol)
        hbond_ats = self.mol.allAtoms.get(lambda x: hasattr(x, 'hbonds'))
        self.assertEquals(len(hbond_ats), 82)
        d = {}
        for a in self.mol.allAtoms:
            if hasattr(a, 'hbonds'):
                for b in a.hbonds:
                    d[b] = 1
        self.assertEquals(len(d.keys()), 28)


        
if __name__ == '__main__':
    test_cases = [
        'HydrogenBondBuilderConstructorTests',
        'HydrogenBondBuilderTests',
    ]

    unittest.main(argv=([__name__, ] + test_cases))

    #unittest.main(argv=([__name__, '-v'] + test_cases))
    #unittest.main()

