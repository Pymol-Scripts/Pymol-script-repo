#########################################################################
#
# Date: Nov. 2002  Author: Daniel Stoffler
#
# Copyright: Daniel Stoffler and TSRI
#
#########################################################################

import sys
from mglutil.regression import testplus
from mglutil.web import HTMLParser
from time import sleep

def pause(sleepTime=0.4):
    ed.master.update()
    sleep(sleepTime)

def test_ParseCGIForms():
    f = open('testfile.html','r')
    txt = f.readlines()
    P = HTMLParser.ParseHTML(mode='forms')
    forms = P.parse(txt)
    f.close()    

    # we have 4 forms in this testfile (INPUT)
    assert len(forms) == 4
    
    # FORM 0:
    # form 0 has a text entry and a sumbit button (INPUT)
    assert len(forms[0].input.keys()) == 2
    # form 0 has no name, the parser should set the name to 'form_0'
    assert forms[0].name == 'form_0'
    # the method of form 0 is 'post'
    assert forms[0].method == 'post'
    # the action of form 0 is http://hoohoo.ncsa.uiuc.edu/cgi-bin/post-query
    assert forms[0].action == 'http://hoohoo.ncsa.uiuc.edu/cgi-bin/post-query' 
    # testing the inputs:

    assert forms[0].input['input_0']['type'] == 'submit'
    assert forms[0].input['input_0']['value'] == 'Submit Query'

    ##########################################################################
    # FORM 1
    # form[1] has 2 text entry, 3 checkbuttons, 7 radiobuttons, 1 submit button

    assert len(forms[1].input.keys()) == 3 # 2 text entries, 1 submit button
    assert len(forms[1].radiobutton.items()) == 2  # 2 categories
    assert len(forms[1].radiobutton['paymethod']) == 5 
    assert len(forms[1].radiobutton['callfirst']) == 2 # 5 + 2 = 7

    assert len(forms[1].checkbutton.items()) == 1 # 1 category
    assert len(forms[1].checkbutton['topping']) == 3 # 3 buttons

    ##########################################################################
    # FORM 2
    # form[2] has 2 comboboxes (SELECT)  and 2 buttons (INPUT)
    assert len(forms[2].select.keys()) == 2
    assert len(forms[2].input.keys())  == 2
    # testing comboboxes:
    assert len(forms[2].select['what-to-do']['options']) == 5
    assert forms[2].select['what-to-do']['options'][0] == 'Drink Coffee'

    assert len(forms[2].select['who-to-do-it-with']['options']) == 6
    assert forms[2].select['who-to-do-it-with']['options'][-1] == 'Chouck'

    ##########################################################################
    # FORM 3
    # form[3] has 3 textareas (TEXTAREA) and 2 buttons (INPUT)

    print forms[3].textarea
    assert len(forms[3].textarea.items()) == 3
    assert len(forms[3].input) == 2
    # testing textareas:
    assert forms[3].textarea['positive']['cols'] == '60'
    assert forms[3].textarea['positive']['rows'] == '20'
    assert forms[3].textarea['username']['text'][0] == 'Your Name Here'
    

harness = testplus.TestHarness( __name__,
                                funs = testplus.testcollect( globals()),
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))
