import os, re

from produtil.testing.utilities import peekable, ListableSet

__all__=[ 'ArithException', 'arithparse', 'ArithKeyError' ]

def arithparse(spec,sets,elements):
    assert(isinstance(spec,basestring))
    tokiter=peekable(_arithtoken(spec))
    return _arithparse_top(tokiter,sets,elements)

class ArithException(Exception):
    """!Raised for parser errors in the inputs to the set arithmetic
    parser (arithparse)."""

class ArithKeyError(ArithException,KeyError):
    """!Raised when an element or set is requested that does not exist"""

########################################################################

def _arithtoken(spec):
    """!Tokenizes a set arithmetic expression.
    @protected
    @see arithparse()"""
    assert(isinstance(spec,basestring))
    for m in re.finditer(r'''(?xs)
        (
            (?P<oper>[A-Za-z]+) \(
          | (?P<endoper> \) )
          | (?P<setname> [A-Za-z][A-Za-z%._0-9]* )
          | (?P<comma> , )
          | (?P<lset> \{ )
          | (?P<rset> \} )
          | (?P<whitespace>\s+)
          | (?P<error> . )
        )''',spec):
        if m.group('whitespace'):
            continue # ignore whitespace
        elif m.group('oper'):
            yield 'oper',m.group('oper')
        elif m.group('endoper'):
            yield ')',m.group('endoper')
        elif m.group('setname'):
            yield 'set',m.group('setname')
        elif m.group('comma'):
            yield ',',','
        elif m.group('lset'):
            yield '{','{'
        elif m.group('rset'):
            yield '}','}'
        elif m.group('error'):
            raise ArithException('Unexpected character %s'%(
                    repr(m.group('error')),))
    yield '',''

########################################################################

def each_not_all(sets):
    for s in sets.iterkeys():
        if s[0]!='*':
            yield s

# Internal implementation routines

def _arithparse_set(tokiter,elements):
    typ,data=tokiter.peek()
    result=ListableSet()
    while typ=='set':
        if data in elements:
            result.add(elements[data])
        else:
            raise ArithKeyError(
                'Unknown test %s. Please select one of: { %s }'%(
                    repr(data), ', '.join([
                            str(k) for k in elements.iterkeys()])))
        tokiter.next() # discard element name
        typ,data=tokiter.peek()
        if typ==',':
            typ,data=tokiter.next() # discard comma and then
            typ,data=tokiter.peek() # peek next token
    if not typ:
        raise ArithException(
            'Unexpected end of text when parsing an '
            'argument list to an operator.')
    if typ!='}':
        raise ArithException(
            'Unexpected %s when parsing an set {} definition '%(repr(data),))
    else:
        tokiter.next()
    return result

def _arithparse_list(tokiter,sets,elements):
    """!Iterates over the arguments to a set arithmetic function,
    yielding sets.  Calls _arithparse_expr() on each argument.

    @protected
    @see arithparse()
    @returns Nothing; this is an iterator."""
    yielded=False
    typ,data=tokiter.peek()
    while typ in ['oper','set','{']:
        result=_arithparse_expr(tokiter,sets,elements)
        assert(isinstance(result,ListableSet))
        yield result
        yielded=True
        typ,data=tokiter.peek()
        if typ==',':
            typ,data=tokiter.next() # discard comma and then
            typ,data=tokiter.peek() # peek next token
    if not typ:
        raise ArithException(
            'Unexpected end of text when parsing an '
            'argument list to an operator.')
    if typ!=')':
        raise ArithException(
            'Unexpected %s when parsing an argument '
            'list to an operator.'%(repr(data),))
    else:
        typ,data=tokiter.next() # discard )

def _arithparse_expr(tokiter,sets,elements):
    """!Evaluates one set arithmetic expression.

    @param tokiter a produtil.testing.utilities.peekable object that
    iterates over the expression.
    @param sets A mapping from set name to
    produtil.testing.utilities.ListableSet objects that represent each
    set."""
    result='invalid value that should be replaced by below code'
    typ,data=tokiter.next()
    if typ=='oper':
        if data=='union':
            # Union of no sets is the empty set:
            result=None
            for subset in _arithparse_list(tokiter,sets,elements):
                if result is None:
                    result=subset
                else:
                    result.union(subset)
            if result is None:
                result=ListableSet()
        elif data=='inter':
            result=None
            for subset in _arithparse_list(tokiter,sets,elements):
                if result is None:
                    result=subset
                else:
                    result.inter(subset)
            if result is None:
                # Intersection of no sets is the empty set:
                result=ListableSet()
        elif data=='minus':
            result=None
            for subset in _arithparse_list(tokiter,sets,elements):
                if result is None:
                    result=subset
                else:
                    result.minus(subset)
        else:
            raise ArithException('Invalid operator %s.  Should be union, '
                            'inter, or minus.'%(data,))
    elif typ=='set':
        if data not in sets:
            raise ArithKeyError('Unknown runset %s.  Known sets: { %s }'%(
                    repr(data), ', '.join([
                            str(k) for k in each_not_all(sets)])))
        result=ListableSet(sets[data])
    elif typ=='{':
        result=_arithparse_set(tokiter,elements)
    else:
        raise ArithException('Unexpected %s in set spec.'%(repr(data),))
    assert(isinstance(result,ListableSet))
    return result

def _arithparse_top(tokiter,sets,elements):
    typ,data = tokiter.peek()
    if not typ:
        return ListableSet()
    return _arithparse_expr(tokiter,sets,elements)
