"""!Common utilities used by other parts of the produtil.testing package."""

import sys, re, StringIO, collections, os, datetime, logging

##@var __all__
# List of variables exported by "from produtil.testing.utilities import *"
__all__=[ 'module_logger', 'BASELINE', 'EXECUTION', 'elipses', 'splitkey', 
          'dqstring2bracestring', 'is_valid_workflow_name', 'unknown_file',
          'peekable','bashify_string', 'ListableSet' ] #, 'yell' ]

# def yell(arg):
#     """!Unused; needs to be removed.

#     @todo Remove this function."""
#     pass 

##@var module_logger
# The default logging.Logger to use if no logger is specified.
module_logger=logging.getLogger('produtil.testing')

##@var BASELINE
# A constant that indicates the suite is being run to generate a new baseline.
BASELINE=object()

##@var EXECUTION
# A constant that indicates the suite is being run to verify against
# an existing baseline.
EXECUTION=object()

##@var unknown_file
# A constant string used for logging purposes to indicate a filename
# was unspecified or unknown.
unknown_file='(**unknown**)'

def bashify_string(string):
        output=StringIO.StringIO()
        for m in re.finditer('''(?xs)
            (
                (?P<quotes>'+)
              | (?P<dquotes>"+)
              | (?P<printable>[!-&(-\[\]-~ ]+)
              | (?P<control>.)
            )''',string):
            if m.group('quotes'):
                output.write('"' + m.group('quotes') + '"')
            elif m.group('dquotes'):
                output.write("'" + m.group('dquotes') + "'")
            elif m.group('printable'):
                output.write("'"+m.group('printable')+"'")
            elif m.group('control'):
                output.write("$'\%03o'"%ord(m.group('control')))
        ret=output.getvalue()
        output.close()
        return ret

def elipses(long_string,max_length=20,elipses='...'):
    """!Returns a shortened version of long_string.

    If long_string is longer than max_length characters, returns a
    string of length max_length that starts with long_string and ends
    with elipses.  Hence, the number of characters of long_string that
    will be used is max_length-len(elipses)

    @param long_string a basestring or subclass thereof
    @param max_length maximum length string to return
    @param elipses the elipses string to append to the end"""
    strlen=len(long_string)
    if strlen<max_length:
        return long_string
    else:
        return long_string[0:max_length-len(elipses)]+elipses

def splitkey(key):
    """!Splits a string on "%" and returns the list, raising an
    exception if any components are empty.

    @returns a list of substrings of key, split on "%"
    @param key a string to split
    @raise ValueError if any substrings are empty"""
    names=key.split("%")
    if any([ not s  for  s in names ]):
        raise ValueError("Empty name component in \"%s\""%(key,))
    return names

def dqstring2bracestring(dq):
    """!Converts a bash-style double quote string to a tripple brace
    string.
    @param dq The bash-style double quote string, minus the 
      surrounding double quotes."""
    output=StringIO.StringIO()
    for m in re.finditer(r'''(?xs)
        (
            \\ (?P<backslashed>.)
          | (?P<braces> [\]\[]+ )
          | (?P<text> [^\\@\]\[]+)
          | (?P<atblock>
                @ \[ @ \]
              | @ \[ ' [^']+ ' \]
              | @ \[ [^\]]+ \]
            )
          | (?P<literal_at> @ (?!\[) )
          | (?P<error> . )
        ) ''',dq):
        if m.group('backslashed'):
            s=m.group('backslashed')
            if s=='@':
                output.write('@[@]')
            elif s in '[]':
                output.write("@['"+s+"']")
            else:
                output.write(s)
        elif m.group('literal_at'):
            output.write('@[@]')
        elif m.group('atblock'):
            output.write(m.group('atblock'))
        elif m.group('braces'):
            output.write("@['"+m.group('braces')+"']")
        elif m.group('text'):
            output.write(m.group('text'))
        else:
            raise ValueError(
                'Cannot convert double-quote string \"%s\" to brace string: '
                'parser error around character \"%s\"."'%(dq,m.group()))
    value=output.getvalue()
    output.close()
    return value

def is_valid_workflow_name(name):
    """!is this a valid name for a produtil.testing workflow?

    Workflow names have to fit within certain restrictions of workflow
    automation suites.  For this reason, we restrict names to begin
    with a letter and only contain letters, numbers and
    underscores.
    @param name the name to check
    @returns True if the name meets requirements and False otherwise"""
    return bool(re.match('(?s)^[a-zA-Z][a-zA-Z0-9_]*$',name))

##@var _HAVE_NOT_PEEKED
# Special constant used by peekable to indicate nothing has been peeked yet.
# @warning Terrible things will happen if you overwrite this.
# @private
_HAVE_NOT_PEEKED=object()

class peekable(object):
    def __init__(self,iterator):
        self.__child=iterator
        self.__iterator=iter(iterator)
        self.__peek=_HAVE_NOT_PEEKED
    @property
    def child(self):
        return self.__child
    def next(self):
        if self.__peek is not _HAVE_NOT_PEEKED:
            p,self.__peek = self.__peek,_HAVE_NOT_PEEKED
        else:
            p=self.__iterator.next()
        return p
    def peek(self):
        if self.__peek is _HAVE_NOT_PEEKED:
            self.__peek=self.__iterator.next()
        return self.__peek
    def at_end(self):
        if self.__peek is not _HAVE_NOT_PEEKED:
            return False
        try:
            self.__peek=self.__iterator.next()
        except StopIteration as se:
            return True
        return False
    def __iter__(self):
        p,self.__peek = self.__peek,_HAVE_NOT_PEEKED
        if p is not _HAVE_NOT_PEEKED:
            yield p
        for v in self.__iterator:
            yield v

class ListableSet(object):
    def __init__(self,contents=None):
        super(ListableSet,self).__init__()
        if contents is None:
            self._set=set()
            self._list=list()
        else:
            self._set=set(contents)
            self._list=list(contents)
        self._NOTHING=object()
    def __contains__(self,item):
        return item in self._set
    def __iter__(self):
        for s in self._list:
            yield s
    def __str__(self):
        return '{ '+', '.join([ repr(s) for s in self._list ])+' }'
    def __repr__(self):
        return 'ListableSet([ '+', '.join(
            [ repr(s) for s in self._list ])+' ])'
    def add(self,item):
        if item not in self:
            self._set.add(item)
            self._list.append(item)
    def minus(self,other):
        for s in other:
            if s in self:
                self._set.discard(s)
                self._list.remove(s)
    def inter(self,other):
        remove=set()
        for s in self:
            if s not in other:
                remove.add(s)
        for s in remove:
            self._set.discard(s)
            self._list.remove(s)
    def union(self,other):
        prior=self._NOTHING
        iprior=-1
        for s in other:
            if s in self:
                prior=s
                iprior=self._list.index(s)
            elif prior is not self._NOTHING:
                self._list.insert(iprior+1,s)
                self._set.add(s)
                prior=s
                iprior=iprior+1
            else:
                self._list.append(s)
                self._set.add(s)
                prior=s
                iprior=len(self._list)
