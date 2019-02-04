from bs4 import BeautifulSoup
import os
import json

webpages = {
            '3d': {'j0': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode5.html',
                   'j2': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode9.html',
                   'j4': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode14.html'},
            '4d': {'j0': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode6.html',
                   'j2': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode10.html',
                   'j4': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode15.html'},
            're': {'j0': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode7.html',
                   'j2': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode11.html',
                   'j4': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode16.html',
                   'j6': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode18.html'},
            'ac': {'j0': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode8.html',
                   'j2': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode12.html',
                   'j4': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode17.html',
                   'j6': 'https://www.ill.eu/sites/ccsl/ffacts/ffactnode19.html'}
            }
            
def _get_webpage_contents_as_string(filename):

    f = open(filename, 'r')
    data = f.read()
    f.close()
    
    return data
    
def _read_htmltable_to_list_from_string(string):
    """
    Adapted from
    https://www.reddit.com/r/Python/comments/2eoeji/how_to_turn_a_beautifulsoup_object_into_a_string/
    """

    _soup_object = BeautifulSoup(string, 'html.parser')
    _table_content = list(_soup_object.table.descendants)
    _table_content_aslist = [str(x.string) for x in _table_content]
    
    # Filter out None-type and newline elements in list
    _table_content_aslist = [x for x in _table_content_aslist if x!='None']
    _table_content_aslist = [x for x in _table_content_aslist if x!='\n']
    
    # Above list contains two of every value. These are filtered out below.
    _table_content_aslist = _table_content_aslist[::2]
    
    # Only return values after element 'D'. j2 for rare earths has a slightly different structure
    # than the others. This operations takes that difference into account.
    index0 = _table_content_aslist.index('D')
    
    return _table_content_aslist[index0+1:]
    
def _list_of_lists_from_tablelist(tablelist, N):

    _listoflists = [tablelist[n:n+N] for n in range(0,len(tablelist),N)]
    
    return _listoflists
    
def _dictionary_from_list_of_lists(lol):

    _newdict = {}
    for e in lol:
        ion = e[0]
        a = [float(x) for x in e[1:-2:2]]
        b = [float(x) for x in e[2:-1:2]]
        c = float(e[-1])
        _newdict[e[0]] = {'a': a, 'b': b, 'c': c}
    
    return _newdict

def _save_dict_to_file(_lol_dict, type):

    f = open('m{}.py'.format(type),'w')
    
    f.write('m{} = '.format(type) + '{\n')
    for e in _lol_dict[type]:
        e = [e[0]]+[float(x) for x in e[1:]]
        f.write("'{}': \n".format(e[0]) + '{\n')
        f.write("'a': [{:>8.4f},{:>8.4f},{:>8.4f}],\n".format(*e[1:-2:2]))
        f.write("'b': [{:>8.4f},{:>8.4f},{:>8.4f}],\n".format(*e[2:-1:2]))
        f.write("'c': {:>8.4f}\n".format(e[-1])+'},\n')
    f.write('}')
    f.close()
    
    
if __name__ == '__main__':
    
    N = 8
    
    _formfactor_pages = [f for f in os.listdir('htmlpages') if f.endswith('.html')]
    
    _lol_dict = {'j0': [], 'j2': [], 'j4': [], 'j6': []}
    
    for f in _formfactor_pages:
        
        Bessel = f.split()[0][1:3]
        
        _content = _get_webpage_contents_as_string('htmlpages\\'+f)
        _content_as_list = _read_htmltable_to_list_from_string(_content)
        _content_as_lol = _list_of_lists_from_tablelist(_content_as_list, N)
        
        _lol_dict[Bessel] += _content_as_lol
    
    for key in _lol_dict.keys():
        print(key)
        _save_dict_to_file(_lol_dict, key)