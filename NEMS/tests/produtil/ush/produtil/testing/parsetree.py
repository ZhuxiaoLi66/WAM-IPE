import sys, re, StringIO, collections, os, datetime, logging, math
import produtil.run, produtil.log, produtil.setup

# This module really does use everything public from utilities and
# tokenize, hence the "import *"
from produtil.testing.utilities import *
from produtil.testing.tokenize import *

__all__=[ 'Context', 'BaseObject', 'null_value', 'TypelessObject', 'Scope',
          'make_params', 'make_scope', 'call_scope', 'Builtin', 'Copy', 
          'CopyDir', 'Link', 'AtParse', 'BitCmp', 'Criteria', 'Filters',
          'Rank', 'SpawnProcess', 'EmbedBash', 'Task', 'Build', 'Platform',
          'Test', 'AutoDetectPlatform', 'Numeric', 'String', 'Environ',
          'Md5Cmp']

class Context(object):
    def __init__(self,scopes,token,run_mode,logger,verbose=True):
        super(Context,self).__init__()
        self.run_mode=run_mode
        self.token=token
        self.scopes=scopes
        if logger is None:
            logger=module_logger
        self.logger=logger
        self.verbose=bool(verbose)
    @property
    def filename(self):
        return self.token.filename
    @property
    def lineno(self):
        return self.token.lineno
    def mpirunner(self,spawnProcess):
        raise NotImplementedError('Base Context class does not know how '
                                  'to run MPI processes.')
    def info(self,message):
        message="%s:%s: %s"%(
                str(self.token.filename),
                repr(self.token.lineno),
                message)
        self.logger.info(message)
        return message
    def warning(self,message):
        message="%s:%s: %s"%(
                str(self.token.filename),
                repr(self.token.lineno),
                message)
        self.logger.warning(message)
        return message
    def error(self,message):
        message="%s:%s: %s"%(
                str(self.token.filename),
                repr(self.token.lineno),
                message)
        self.logger.error(message)
        return message

def fileless_context(scopes=None,filename=None,lineno=None,
                     run_mode=None,logger=None,verbose=True):
    if scopes is None: scopes=[]
    if filename is None: filename=unknown_file
    if lineno is None: lineno=-1
    if logger is None: logger=module_logger
    if run_mode is None: run_mode=EXECUTION
    return Context(
        scopes,Token(end_of_line_type,end_of_line_type,filename,lineno),
        run_mode,logger,verbose)

class BaseObject(object):
    def __init__(self,defscopes):
        for s in defscopes:
            assert(isinstance(s,Scope))
        self.defscopes=defscopes
        self.is_scope=False
        self.is_filters=False
        self.is_criteria=False
        self.is_scalar=False
        self.can_be_used=False
        self.is_test=False
    def bash_context(self,con):
        raise Exception("Cannot express null_value in a bash string.")
    def is_valid_rvalue(self,con): return True
    def string_context(self,con): return "null"
    def logical_context(self,con): return False
    def numeric_context(self,con): return 0.0
    def _apply_rescope(self,scopemap=None,prepend=None):
        if not prepend: prepend=[]
        if not scopemap: scopemap={}
        self.defscopes=prepend+[ scopemap[s] if s in scopemap else s
                                 for s in self.defscopes ]
    def rescope(self,scopemap=None,prepend=None):
        return self
    def run(self,con): 
        con.logger.info(self.string_context(con))
    def iterdeps(self):
        return
        yield 'a' # Syntactic trick to ensure this is an iterator.
    def __repr__(self):
        return '<BaseObject@%s>'%(id(self),)
    def __str__(self):
        return '<BaseObject@%s>'%(id(self),)

##@var null_value
# A special constant that indicates a variable without a value. 
# @warning Terrible things will happen if you overwrite this.
null_value=BaseObject([])

########################################################################

class TypelessObject(BaseObject):
    """!Represents an object that cannot be evauated in any context.
    This is a convenience class intended to be used by subclasses to
    disable all but certain contexts."""
    def bash_context(self,con):
        raise TypeError('Cannot evaluate %s in a bash context.'%(
                type(self).__name__,))
    def string_context(self,con):
        raise TypeError('Cannot evaluate %s in a string context.'%(
                type(self).__name__,))
    def logical_context(self,con):
        raise TypeError('Cannot evaluate %s in a logical context.'%(
                type(self).__name__,))
    def numeric_context(self,con):
        raise TypeError('Cannot evaluate %s in a numeric context.'%(
                type(self).__name__,))
    def run(self,con):
        raise TypeError('Cannot run objects of type %s.'%(
                type(self).__name__,))

########################################################################

class Scope(BaseObject):
    def __init__(self,defscopes=None):
        if defscopes is None:
            defscopes=list()
        super(Scope,self).__init__(defscopes)
        self.__const=dict()
        self.__vars=dict()
        self.__parameters=dict()
        self.is_scope=True
        self.can_be_used=True
        self.__overrides=collections.defaultdict(list)

    def __eq__(self,other):
        # Scopes are equal iff they are the same object.
        return self is other

    def override_local(self,defscopes,key,val):
        names=splitkey(key)
        the_string=String(self.defscopes+[self],val,False)
        #print 'overrides %s%%%s = %s'%(repr(names[0]),repr(names[1:]),repr(the_string))

        self.__overrides[names[0]].append([ names[1:],the_string])

    def validate_parameter(self,name): pass

    def new_empty(self): return Scope(self.defscopes)

    def bash_context(self,con):
        raise Exception("Cannot express a hash in a bash string.")

    def string_context(self,con): 
        if self.haslocal('_as_string'):
            value=self.getlocal('_as_string')
            return value.string_context(con)
        else:
            return str(id(self))

    def no_nulls(self):
        for k,v in self.iterlocal():
            if v is null_value:
                return False
        return True

    def _set_parameters(self,update):
        self.__parameters.update(update)
        for p in self.__parameters.iterkeys():
            self.validate_parameter(p)

    def _set_constants(self,update):
        self.__const.update(update)
        for p in self.__const.iterkeys():
            self.validate_parameter(p)

    def numeric_context(self,con):
        return len(self.__vars) + len(self.__parameters) + len(self.__const)
    def logical_context(self,con): 
        return bool(self.__vars) or bool(self.__parameters) or bool(self.__const)

    def get_type(self,key):
        if key in self.__parameters:
            return 'parameter'
        elif key in self.__const:
            return 'constant'
        elif key in self.__vars:
            return 'var'
        else:
            return None

    def as_parameters(self,con):
        """!Changes all variables to parameters."""
        self.__parameters.update(self.__vars)
        self.__vars=dict()
        for p in self.__parameters.iterkeys():
            self.validate_parameter(p)
        return self

    def _apply_rescope(self,scopemap=None,prepend=None):
        super(Scope,self)._apply_rescope(scopemap,prepend)
        for d in [ self.__parameters, self.__vars, self.__const]:
            for k,v in d.iteritems():
                v._apply_rescope(scopemap,prepend)

    def rescope(self,scopemap=None,prepend=None):
        scope=self.new_empty()
        scope._apply_rescope(scopemap,prepend)
        for k,v in self.__parameters.iteritems():
            scope.__parameters[k]=v.rescope(scopemap,prepend)
        for k,v in self.__vars.iteritems():
            scope.__vars[k]=v.rescope(scopemap,prepend)
        for k,v in self.__const.iteritems():
            scope.__const[k]=v.rescope(scopemap,prepend)
        return scope

    def has_parameters(self):
        return bool(self.__parameters)

    def has_constants(self):
        return bool(self.__const)

    def use_from(self,used_scope,only_scalars=False):
        if used_scope.has_parameters():
            raise Exception('Cannot "use" a function.  Fixme: use better '
                            'exception here.')
        if not used_scope.is_scope:
            raise Exception('Target of "use" statement is not a scope.')
        if not used_scope.can_be_used:
            raise Exception('Target of "use" statement cannot be used.')
        prepend_me=[ self ]
        found_non_scalars=False
        for k,v in used_scope.iterlocal():
            if only_scalars and v.is_scope:
                found_non_scalars=True
            if k not in self.__const:
                my_v=v.rescope({used_scope:self})
                self.force_define(k,my_v)

        if used_scope.haslocal('RUNDIR'):
            mycon=fileless_context([self])
            usedcon=fileless_context([used_scope])
                
        return found_non_scalars
    def apply_parameters(self,scope,con):
        assert(not scope.__parameters)
        assert(self.__parameters)
        if not self.__parameters:
            return self
        s=self.new_empty()
        s._apply_rescope({self:s,scope:s})
        caller=con.scopes[0]
        #print('APPLY PARAMETERS from %s to %s type %s\n'%(
        #        repr(scope),repr(self),type(s).__name__))
        for k,v in self.__parameters.iteritems():
            if scope.haslocal(k):
                s.__vars[k]=scope.getlocal(k).rescope({self:s,scope:s})
            elif v is not null_value:
                s.__vars[k]=v.rescope({self:s,scope:s})
            else:
                raise Exception('%s: no argument sent for this parameter'%(
                        k,))
        for k,v in self.__vars.iteritems():
            if k not in s.__vars:
                s.__vars[k]=v.rescope({self:s,scope:s})
        for k,v in self.__const.iteritems():
            if k not in s.__const:
                s.__const[k]=v.rescope({self:s,scope:s})
        #print('APPLY RESULT IS %s %s\n'%(
        #        type(s).__name__,repr(s)))
        return s

    def __str__(self):
        return '{' + ','.join( [
                "%s=%s"%(str(s),repr(k)) for s,k in self.iterlocal()
                ] ) + '}'

    def __repr__(self):
        return '{' + ','.join( [
                "%s=%s"%(str(s),repr(k)) for s,k in self.iterlocal()
                ] ) + '}'

    def subscope(self,key):
        """!Returns a Scope within this Scope, with the given name.
        
        @param key a valid identifier within this scope
        @returns a Scope with the given name, within this Scope
        @raise ValueError if the key contains a "%"
        @raise TypeError if the key refers to something in this Scope
          that is not a Scope.  This is detected through the is_scope 
          attribute or property."""
        if "%" in key:
            raise ValueError("Key \"%s\" is not a valid identifier"%(key,))
        if key in self.__parameters:
            value=self.__parameters[key]
        elif key in self.__const:
            value=self.__const[key]
        else:
            value=self.__vars[key]
        try:
            if value.is_scope:
                return value
        except AttributeError as ae:
            pass # value does not define is_scope
        raise TypeError("Key \"%s\" refers to something that is not a Scope."
                        %(key,))

    def getlocal(self,key):
        """!Return the value of a key local to this scope without
        searching other scopes.  Will raise ValueError if the key
        contains a "%" 

        @param key a valid identifier within this scope
        @returns The value of the key from this scope.
        @raise ValueError if the key is syntactically not a valid identifier
           such as one that contains a "%"
        @raise KeyError if the key is a valid identifier but is not
        within this scope."""
        if "%" in key:
            raise ValueError("Key \"%s\" is not a valid identifier"%(key,))
        if key in self.__parameters:
            return self.__parameters[key]
        elif key in self.__const:
            return self.__const[key]
        elif key in self.__vars:
            return self.__vars[key]
        raise KeyError(key)

    def setlocal(self,key,value):
        """!Sets the value of a key within this scope.

        @param key a valid identifier to set within this scope
        @param valuel the value of the identifier
        @raise ValueError if the key is not a valid identifier, such as
          one that contains a "%" """
        if not isinstance(value,BaseObject):
            raise TypeError('The value argument to setlocal must be a '
                            'subclass of BaseObject, not a %s %s'%(
                    type(value).__name__,repr(value)))
        if '%' in key:
            raise ValueError("Key \"%s\" contains a \"%\""%(key,))
        if key in self.__parameters:
            raise Exception("Key \"%s\" is already a parameter.  FIXME: "
                            "Need better exception class."%(key,))
        elif key in self.__const:
            raise Exception("Key \"%s\" is already a constant.  FIXME: "
                            "Need better exception class."%(key,))
        self.__vars[key]=value
        global over
        try:
            if key in self.__overrides:
                for names,override in self.__overrides[key]:
                    #print 'Check override %s %s'%(names[-1],override)
                    subscope=self
                    for name in [key]+names[0:-1]:
                        subscope=subscope.resolve(name)
                        assert(subscope.is_scope)
                        assert(isinstance(subscope,produtil.testing.parsetree.Scope))
                        if not subscope.is_scope:
                            continue # raise KeyError(name)
                    assert(isinstance(override,BaseObject))
                    subscope.force_define(names[-1],override)
                return override
        except KeyError as ke:
            raise
        return value

    def haslocal(self,key):
        return key in self.__parameters or key in self.__vars or \
            key in self.__const

    def iterlocal(self):
        for k,v in self.__parameters.iteritems():
            yield k,v
        for k,v in self.__const.iteritems():
            yield k,v
        for k,v in self.__vars.iteritems():
            yield k,v

    def resolve(self,key,scopes=None):
        assert(isinstance(key,basestring))
        names=splitkey(key)
        con=fileless_context()
        if scopes is None:
            search=[self]+self.defscopes
        else:
            search=[self]+scopes
        #scopestack=list()
        # yell('search for %s = %s in %s\n'%(repr(key),repr(names),repr(search)))
        for name in names:
            found=None
            for scope in search:
                try:
                    found=scope.getlocal(name)
                    # if name=='TEST_NAME':
                    #     print 'Found %s in scope %s = "%s" name "%s"'%(
                    #         str(name),id(scope),elipses(scope.getlocal(
                    #                 'TEST_DESCR').string_context(con)),
                    #         scope.getlocal('TEST_NAME').string_context(con))
                    #scopestack.insert(0,scope)
                    break # Done searching scopes.
                except KeyError as ke:
                    # yell('Key %s not in scope %s from top %s\n'%(
                    #         repr(name),repr(scope),repr(self)))
                    continue # Check for name in next scope.
            if found is None:
                raise KeyError(key)
            search=[found]
        if found is None:
            raise KeyError(key)
        # if subscopes: 
        #     return ( found, scopestack )
        # else:
        return found

    def force_define(self,key,value,skip_constants=False):
        names=splitkey(key)
        lval=self
        for i in xrange(len(names)-1):
            lval=lval.getlocal(names[i])
        if skip_constants and lval.get_type(names[-1])=='constant':
            #module_logger.info('Do not redefine %s.'%(key,))
            return None
        assert(names[-1]!='TEST_NAME')
        lval.setlocal(names[-1],value)
        return value

    def check_define(self,key,value):
        names=splitkey(key)
        lval=self
        for i in xrange(len(names)-1):
            lval=lval.subscope(names[i])
        assert(names[-1]!='TEST_NAME')
        if lval.haslocal(value):
            raise Exception('Key %s already exists. (FIXME: Put better '
                            'exception class here.)'%(key,))
        #yell('setlocal %s = %s\n'%(names[-1],value))
        lval.setlocal(names[-1],value)
        return value

    def expand_string(self,string,con,scopes=None):
        stream=StringIO.StringIO()
        # if string.find('TEST_NAME')>-1:
        #     print 'Expand "%s"'%(elipses(string,max_length=80),)
        # yell('Expand %s in %s\n'%(repr(string),repr(self)))
        def streamwrite(s):
            #yell("Append \"%s\" to output string.\n"%(s,))
            stream.write(s)
        for m in re.finditer(r'''(?sx)
            (
                (?P<text>[^@]+)
              | @ \[ (?P<escaped_at>@ ) \]
              | @ \[ ' (?P<escaped_text> [^']+ ) ' \]
              | @ \[ (?P<varexpr>[^\]]+) \]
              | (?P<literal_at>@ ) (?! \[ )
              | (?P<error>.)
            )''',string):
            if m.group('text'):
                streamwrite(m.group())
            elif m.group('escaped_text'):
                streamwrite(m.group('escaped_text'))
            elif m.group('escaped_at'):
                streamwrite(m.group('escaped_at'))
            elif m.group('literal_at'):
                streamwrite(m.group('literal_at'))
            elif m.group('varexpr'):
                streamwrite(self.resolve(m.group('varexpr'),scopes) \
                                .string_context(con))
            else:
                raise ValueError("Parser error: invalid character \"%s\" in"
                                 " \"%s\"\n"%(m.group(0),string))
        val=stream.getvalue()
        # if string.find('TEST_NAME')>-1:
        #     print 'Result "%s"'%(elipses(val,max_length=120),)
            #assert(val.find('nmm_cntrl')==-1)
        stream.close()
        return val

def make_params(defscopes,**kwargs):
    s=Scope(defscopes)
    for k,v in kwargs.iteritems():
        s.setlocal(k,v)
    return s.as_parameters()

def make_scope(defscopes,**kwargs):
    s=Scope(defscopes)
    for k,v in kwargs.iteritems():
        s.setlocal(k,v)
    return s

def call_scope(scope,con,defscopes,**kwargs):
    parms=make_scope(defscopes,**kwargs)
    assert(parms.no_nulls())
    s=scope.apply_parameters(parms,con)
    assert(s.no_nulls())
    return s
    
########################################################################

class Builtin(Scope):
    def __init__(self,defscopes,opname):
        super(Builtin,self).__init__(defscopes)
        self.__opname=opname
    def bash_context(self,con):
        raise TypeError('Cannot evaluate built-in operator %s in a '
                        'bash context.'%(self.__opname,))
    def string_context(self,con):
        raise TypeError('Cannot evaluate built-in operator %s in a '
                        'string context.'%(self.__opname,))
    def logical_context(self,con):
        raise TypeError('Cannot evaluate built-in operator %s in a '
                        'logical context.'%(self.__opname,))
    def numeric_context(self,con):
        raise TypeError('Cannot evaluate built-in operator %s in a '
                        'numeric context.'%(self.__opname,))

    def getcom(self,con):
        """!Returns the COM directory for the test that contains this
        operation.

        Walks up the nested scopes from inner to outer, searching for
        the innermost Test (detected via BaseObject.is_test).  When
        found, asks the test to resolve the "COM" symbol, and returns
        the result of that objects string_context function.

        @param con a Context from which this function is called

        @returns The COM directory for the innermost test in the
        defscopes as a string."""
        for scope in self.defscopes:
            if scope.is_test:
                return scope.resolve('COM').string_context(con)
        raise KeyError('COM')

    def run(self,con):
        raise TypeError('Cannot run built-in operator %s.'%(self.__opname,))

    def new_empty(self):
        return Builtin(self.defscopes,self.__opname)

########################################################################

class Copy(Builtin):
    def __init__(self,defscopes,empty=False):
        super(Copy,self).__init__(defscopes,'.copy.')
        if not empty:
            self._set_parameters({'src':null_value, 'tgt':null_value})
    def new_empty(self):
        return Copy(self.defscopes,empty=True)
    def run(self,con):
        assert(self.no_nulls())
        src=self.resolve('src').string_context(con)
        tgt=self.resolve('tgt').string_context(con)
        produtil.fileop.deliver_file(src,tgt)
    def bash_context(self,con):
        assert(self.no_nulls())
        src=self.resolve('src').bash_context(con)
        tgt=self.resolve('tgt').bash_context(con)
        return 'deliver_file %s %s\n'%(src,tgt)

########################################################################

class CopyDir(Builtin):
    def __init__(self,defscopes,empty=False):
        super(CopyDir,self).__init__(defscopes,'.copydir.')
        if not empty:
            self._set_parameters({'src':null_value, 'tgt':null_value})
    def new_empty(self):
        return CopyDir(self.defscopes,empty=True)
    def run(self,con):
        raise NotImplementedError('CopyDir.run is not implemented yet.')
    def bash_context(self,con):
        assert(self.no_nulls())
        src=self.resolve('src').bash_context(con)
        tgt=self.resolve('tgt').string_context(con)
        return '''
(
  shopt -u failglob ;
  shopt -s nullglob ;
  for srcfile in %s/%s ; do
    deliver_file "$srcfile" . ;
  done
)
'''%(src,tgt)

########################################################################

class Link(Builtin):
    def __init__(self,defscopes,empty=False):
        super(Link,self).__init__(defscopes,'.link.')
        if not empty:
            self._set_parameters({'src':null_value, 'tgt':null_value})
    def new_empty(self):
        return Link(self.defscopes,empty=True)
    def run(self,con):
        assert(self.no_nulls())
        src=self.resolve('src').string_context(con)
        tgt=self.resolve('tgt').string_context(con)
        produtil.fileop.make_symlink(src,tgt)
    def bash_context(self,con):
        assert(self.no_nulls())
        src=self.resolve('src').bash_context(con)
        tgt=self.resolve('tgt').bash_context(con)
        return 'rm -f %s\nln -sf %s %s\n'%(tgt,src,tgt)

########################################################################

class AtParse(Builtin):
    def __init__(self,defscopes,empty=False):
        super(AtParse,self).__init__(defscopes,'.atparse.')
        if not empty:
            self._set_parameters({'src':null_value, 'tgt':null_value})
    def new_empty(self):
        return AtParse(self.defscopes,empty=True)
    def run(self,con):
        raise NotImplementedError("FIXME: Sam has not implemented AtParse.run")
    def bash_context(self,con):
        out=StringIO.StringIO()
        src=self.resolve('src').bash_context(con)
        tgt=self.resolve('tgt').bash_context(con)
        out.write("echo input to atparse from %s:\ncat %s\n"%(src,src))
        out.write("echo send to %s\n"%(tgt,))
        out.write("cat %s | atparse \\\n"%(src,))
        seen=set()
        for scope in self.defscopes:
            for k,v in scope.iterlocal():
                if k in seen: continue
                seen.add(k)
                if '%' in k or '.' in k:
                    pass#out.write('# $%s: skip; invalid shell variable name\n'%(k,))
                elif k[0:2] == '__':
                    pass#out.write("# $%s: skip; name begins with __\n"%(k,))
                elif v is null_value:
                    pass#out.write('# $%s: skip; has no value\n'%(k,))
                elif not v.is_scalar:
                    pass#out.write('# $%s: skip; value is not scalar\n'%(k,))
                else:
                    out.write('  %s=%s \\\n'%(k,v.bash_context(con)))
        out.write("  > %s\n"%(tgt,))
        if con.verbose:
            out.write("set -xe\n")
        else:
            out.write("set -e\n")
        out.write('cat %s\n\n'%(tgt,))
        ret=out.getvalue()
        out.close()
        return ret

########################################################################

class BitCmp(Builtin):
    def __init__(self,defscopes,empty=False):
        super(BitCmp,self).__init__(defscopes,'.bitcmp.')
        if not empty:
            self._set_parameters({'src':null_value, 'tgt':null_value})
    def new_empty(self):
        return BitCmp(self.defscopes,empty=True)
    def run(self,con):
        src=self.resolve('src').string_context(con)
        tgt=self.resolve('tgt').string_context(con)
        if con.run_mode==BASELINE:
            produtil.fileop.deliver_file(src,tgt)
            return
        if os.path.samefile(src,tgt):
            # Same file object in filesystems.
            return True
        with open(src,'rt') as srcf:
            with open(tgt,'rt') as tgtf:
                srcstat=os.fstat(src.fileno())
                tgtstat=os.fstat(tgt.fileno())
                if not srcstat: return False # file stopped existing
                if not tgtstat: return False # file stopped existing
                if srcstat.st_size!=tgtstat.st_size:
                    # Different size according to stat
                    return False
                eof=False
                while not eof:
                    srcdat=src.read(1048576)
                    tgtdat=tgt.read(1048576)
                    if len(srcdat)!=len(tgtdat):
                        return False # Lengths differ
                    if srcdat!=tgtdat:
                        return False # Contents differ
                    eof=not len(srcdat) or not len(tgtdat)
                return True
    def bash_context(self,con):
        src=self.resolve('src').bash_context(con)
        tgt=self.resolve('tgt').bash_context(con)
        compath=os.path.join(self.getcom(con),
                             os.path.basename(tgt))
        if con.run_mode==BASELINE:
            return 'deliver_file %s %s\n'%(tgt,src) + \
                   'deliver_file %s %s\n'%(tgt,compath)
        else:
            return '''deliver_file %s %s
set +e # bitcmp failures are reported by report_failed\nbitcmp %s %s\nset -e\n'''%(
                tgt,compath, # deliver_file $tgt $compath
                tgt,src)     # bitcmp $tgt $src

########################################################################

class Md5Cmp(Builtin):
    def __init__(self,defscopes,empty=False):
        super(Md5Cmp,self).__init__(defscopes,'.md5cmp.')
        if not empty:
            self._set_parameters({'src':null_value, 'tgt':null_value})
    def new_empty(self):
        return Md5Cmp(self.defscopes,empty=True)
    def run(self,con):
        raise NotImplementedError("Md5Cmp.run is not implemented.")
    def bash_context(self,con):
        md5ref=self.resolve('src').string_context(con) # reference md5sum
        exe=self.resolve('tgt').string_context(con) # executable
        md5sum=os.path.join(self.getcom(con),
                            os.path.basename(md5ref))
        md5ref=bashify_string(md5ref)
        exe=bashify_string(exe)
        md5sum=bashify_string(md5sum)
        # Behaves the same way in baseline and verification mode
        # because this is checking the executable used during the run
        return '''md5sum {exe} > {md5sum}
report_line md5sum: $( cat {md5sum} )
report_line md5sum "local="{md5sum}
report_line md5sum "reference="{md5ref}
'''.format(
            md5ref=md5ref,md5sum=md5sum,exe=exe)

########################################################################

class Criteria(TypelessObject):
    def __init__(self,defscopes):
        super(Criteria,self).__init__(defscopes)
        self.__opmap=collections.defaultdict(list)
        self.__tgtlist=list()
        self.is_criteria=True
    def add_binary_operator(self,tgt,op,src,con):
        if tgt not in self.__opmap:
            self.__tgtlist.append(tgt)
        callme=call_scope(op,con,self.defscopes,tgt=tgt,src=src)
        for mycall in self.__opmap[tgt]:
            assert mycall is not null_value
            if mycall==callme:
                return
        callme=self.__opmap[tgt].append(callme)
        assert callme is not null_value
    def rescope(self,scopemap=None,prepend=None):
        if prepend is None: prepend=[]
        if scopemap is None: scopemap={}
        f=Criteria(prepend+[ scopemap[s] if s in scopemap else s
                        for s in self.defscopes ])
        for tgt in self.__tgtlist:
            rtgt=tgt.rescope(scopemap,prepend)
            f.__tgtlist.append(rtgt)
            f.__opmap[rtgt]=[ 
                callme.rescope(scopemap,prepend)
                for callme in self.__opmap[rtgt] ]
        return f
    def _apply_rescope(self,scopemap=None,prepend=None):
        super(Criteria,self)._apply_rescope(scopemap,prepend)
        for tgt in self.__tgtlist:
            tgt._apply_rescope(scopemap,prepend)
            for callme in self.__opmap[tgt]:
                callme._apply_rescope(scopemap,prepend)
    def use_from(self,criteria,only_scalars=False):
        if only_scalars:
            raise ValueError('In Criteria.use_from, only_scalars must '
                             'be False.')
        if not criteria.is_criteria:
            raise TypeError('Criteria blocks can only use criteria blocks.')
        for tgt,callme in criteria.itercriteria():
            found=False
            for mycall in self.__opmap[tgt]:
                if callme==mycall:
                    found=True
                    break
            if not found: 
                self.__opmap[tgt].append(callme)
    def itercriteria(self):
        for tgt in self.__tgtlist:
            for callme in self.__opmap[tgt]:
                yield tgt,callme
    def bash_context(self,con):
        out=StringIO.StringIO()
        if con.run_mode==BASELINE:
            out.write('\n########################################################################\necho BASELINE GENERATION:\n\n')
        else:
            out.write('\n########################################################################\necho OUTPUT VALIDATION:\n\n')
        for tgt in self.__tgtlist:
            # out.write('echo criteria for target %s:\n'%(
            #         tgt.bash_context(con),))
            for callme in self.__opmap[tgt]:
                out.write(callme.bash_context(con))
                if con.run_mode==BASELINE: break
        if con.run_mode==BASELINE:
            out.write('\necho END OF BASELINE GENERATION\n########################################################################\n\n')
        else:
            out.write('\necho END OF OUTPUT VALIDATION\n########################################################################\n\n')
        ret=out.getvalue()
        out.close()
        return ret

########################################################################

class Filters(TypelessObject):
    def __init__(self,defscopes):
        super(Filters,self).__init__(defscopes)
        self.__opmap=dict()
        self.__tgtlist=list()
        self.is_filters=True
    def add_binary_operator(self,tgt,op,src,con):
        if tgt not in self.__opmap:
            self.__tgtlist.append(tgt)
        self.__opmap[tgt]=call_scope(op,con,self.defscopes,
                                     tgt=tgt,src=src)
        assert self.__opmap[tgt] is not null_value
    def rescope(self,scopemap=None,prepend=None):
        if prepend is None: prepend=[]
        if scopemap is None: scopemap={}
        f=Filters(prepend+[ scopemap[s] if s in scopemap else s
                        for s in self.defscopes ])
        for tgt in self.__tgtlist:
            rtgt=tgt.rescope(scopemap,prepend)
            f.__tgtlist.append(rtgt)
            f.__opmap[rtgt]=self.__opmap[tgt].rescope(scopemap,prepend)
        return f
    def _apply_rescope(self,scopemap=None,prepend=None):
        super(Filters,self)._apply_rescope(scopemap,prepend)
        for tgt in self.__tgtlist:
            tgt._apply_rescope(scopemap,prepend)
    def use_from(self,filters,only_scalars=False):
        if only_scalars:
            raise ValueError('In Filters.use_from, only_scalars must '
                             'be False.')
        if not filters.is_filters:
            raise TypeError('Filters blocks can only use filters blocks.')
        for tgt,callme in filters.iterfilters():
            have_tgt=tgt in self.__opmap
            if have_tgt and self.__opmap[tgt]==callme:
                continue
            self.__opmap[tgt]=callme
            assert self.__opmap[tgt] is not null_value
            if not have_tgt:
                self.__tgtlist.append(tgt)
    def iterfilters(self):
        for tgt in self.__tgtlist:
            yield tgt,self.__opmap[tgt]
    def bash_context(self,con):
        out=StringIO.StringIO()
        out.write('\n########################################################################\necho INPUT FILTERS:\n\n')
        for tgt in self.__tgtlist:
            # out.write('echo Filter for target %s:\n'%(
            #         tgt.bash_context(con),))
            out.write(self.__opmap[tgt].bash_context(con))
        out.write('\necho END OF INPUT FILTERS\n########################################################################\n\n')
        ret=out.getvalue()
        out.close()
        return ret

########################################################################

class Rank(TypelessObject):
    def __init__(self,args,opts):
        self.__args=args
        self.__opts=opts
    def __repr__(self):
        return 'Rank(args=%s,opts=%s)'%(
            repr(self.__args),repr(self.__opts))
    @property
    def args(self):
        return self.__args
    def argiter(self):
        for arg in self.__args:
            yield arg
    def ranks(self,con):
        if self.__opts is None or not self.__opts.haslocal('ranks'):
            return 0
        return int(self.__opts.getlocal('ranks').numeric_context(con))
    def ppn(self,con):
        if self.__opts is None or not self.__opts.haslocal('ppn'):
            return 0
        return int(self.__opts.getlocal('ppn').numeric_context(con))
    def threads(self,con):
        if self.__opts is None or not self.__opts.haslocal('threads'):
            return 0
        return int(self.__opts.getlocal('threads').numeric_context(con))

########################################################################

def pack_ranks(nodesize,count):
    out=list()
    n=count
    if n<nodesize:  # special case: smaller than one node
        out.append([1,n])
    elif nodesize*(n//nodesize)==n: #exact number of nodes
        out.append([n//nodesize,nodesize])
    else:
        need=math.ceil(n/float(nodesize))
        averagef=n/math.ceil(need)
        af,ai = math.modf(averagef)
        n1=need-round(af*need)
        n2=round(af*need)
        if n1: out.append([n1, ai])
        if n2: out.append([n2, ai+1])
    return out

########################################################################

class SpawnProcess(TypelessObject):
    def __init__(self,defscopes):
        super(SpawnProcess,self).__init__(defscopes)
        self.__ranks=list()
    def add_rank(self,args,opts):
        self.__ranks.append(Rank(args,opts))
    def iterrank(self):
        for rank in self.__ranks:
            yield rank
    def mpi_comm_size(self,con):
        size=0
        for rank in self.__ranks:
            size+=rank.ranks(con)
        return size

    def rocoto_resources(self,con):
        nodes=list()
        accumranks=0
        accumthreads=0
        accumppn=0
        nodesize=int(self.defscopes[-1].resolve('plat%cores_per_node')
                     .numeric_context(con))
        if nodesize<1:
            raise Exception(con.error(
                    'Value of plat%cores_per_node must be an integer >= 1'))

        def accumulate():
            con.info('enter accumulate with ppn=%d ranks=%d threads=%d nodesize=%d'%(
                    accumppn,accumranks,accumthreads,nodesize))
            if not accumranks:
                con.info('No ranks.')
                return
            if accumppn:
                con.info('accumulate ppn=%d ranks=%d'%(accumppn,accumranks))
                ranks=pack_ranks(accumppn,accumranks)
            elif accumthreads:
                con.info('accumulate nodesize=%d threads=%d ranks=%d'%(
                        nodesize,accumthreads,accumranks))
                ranks=pack_ranks(nodesize//accumthreads,accumranks)
            else:
                con.info('accumulate nodesize=%d ranks=%d'%(
                        nodesize,accumranks))
                ranks=pack_ranks(nodesize,accumranks)
            nodes.extend(ranks)

        missing_ranks=False
        for rank in self.__ranks:
            nranks=max(0,rank.ranks(con))
            nthreads=max(0,rank.threads(con))
            ppn=max(0,rank.ppn(con))
            if nthreads==2 and ppn!=8:
                con.info(repr(rank))
            if not nranks:
                missing_ranks=True
            elif accumranks:
                if accumthreads==nthreads and accumppn==ppn:
                    # New block of MPI ranks not the same distribution
                    # as previous:
                    accumulate()
                    accumranks,accumthreads,accumppn = nranks,nthreads,ppn
                else:
                    # New block of MPI ranks has the same distribution
                    # as the previous, so we can just increase the
                    # number of ranks:
                    accumranks+=nranks
                    # FIXME: update parser to give a way to override this.
            else: # First block of MPI ranks seen so far:
                accumranks,accumthreads,accumppn = nranks,nthreads,ppn

        if accumranks: accumulate()

        MPI=con.scopes[-1].resolve('plat%MPI').string_context(con)
        con.info('MPI=%s'%(repr(MPI),))
        if MPI.upper().find('LSF')>=0:
            bfppn=int(con.scopes[-1].resolve('plat%cores_per_node').
                      numeric_context(con))//max(1,nthreads)
            if ppn>0: bfppn=min(ppn,bfppn)
            before='<native>-R span[ptile=%d]</native>'%(bfppn,)
        else:
            before=''

        if missing_ranks:
            #if world:
            raise Exception(con.error('Missing a ranks= value in an '
                                          'MPI program.'))
            # Serial program.  Request 2 cores to ensure exclusive access.
            return before+'<cores>2</cores>'

        if len(nodes)==1 and nodes[0][1]==nodesize:
            return before+'<cores>%d</cores>'%(nodes[0][1]*nodes[0][0])
        else:
            return '<nodes>%s</nodes>'%(
                '+'.join(['%d:ppn=%d'%(node[0],node[1]) for node in nodes]))
        
    def bash_context(self,con):
        out=StringIO.StringIO()
        out.write('# Embedded process execution:\n')
        need_ranks=len(self.__ranks)>1
        have_ranks=False
        ranks=list()
        threads=list()
        for rank in self.__ranks:
            nranks=rank.ranks(con)
            if nranks>0: have_ranks=True
            if nranks<1 and need_ranks:
                nranks=1
            ranks.append(nranks)
            threads.append(rank.threads(con))
            con.info('Rank block: threads=%s ranks=%s'%(repr(threads[-1]),
                                                          repr(ranks[-1])))
        nthreads=max(0,max(threads))
        assert(nthreads>0)
        out.write('export OMP_NUM_THREADS=%d MKL_NUM_THREADS=0\n'%(
                nthreads,))
        if not have_ranks:
            # Serial or openmp program.
            out.write(' '.join([r.bash_context(con) 
                                for r in self.__ranks[0].args]))
            out.write('\n')
        elif len(self.__ranks)==1:
            out.write('%s\n'%(con.mpirunner(self),))
            #out.write('mpirun -np %d '%(int(self.__ranks[0].ranks(con))))
            #out.write(' '.join([r.bash_context(con)
            #                    for r in self.__ranks[0].args]))
            #out.write('\n')
        else:
            raise Exception("FIXME: Sam has not implemented MPMD yet.")
        out.write('# End of embedded process execution.\n')
        ret=out.getvalue()
        out.close()
        return ret

########################################################################

class EmbedBash(Scope):
    def __init__(self,defscopes):
        super(EmbedBash,self).__init__(defscopes)
        self.__template=None

    def validate_parameter(self,name):
        pass
        #if not re.match('(?s)^[a-zA-Z][a-zA-Z0-9_]*$',name):
            #raise ValueError('Invalid bash variable name $%s FIXME: use better exception here'%(name,))

    def bash_context(self,con):
        raise Exception("Cannot express a bash script in a bash string.")

    def _apply_rescope(self,scopemap=None,prepend=None):
        super(EmbedBash,self)._apply_rescope(scopemap,prepend)
        self.__template._apply_rescope(scopemap,prepend)
        self.__template.defscopes=[self]+self.defscopes

    def rescope(self,scopemap=None,prepend=None):
        s=super(EmbedBash,self).rescope(scopemap,prepend)
        s.__template=self.__template.rescope(scopemap,prepend)
        s.__template.defscopes=[self]+self.defscopes
        return s

    def __str__(self):
        return "bash script \"%s\" %s"%(
            elipses(repr(self.gettemplate())),
            super(EmbedBash,self).__str__())

    def __repr__(self):
        return "bash script \"%s\" %s"%(
            elipses(repr(self.gettemplate())),
            super(EmbedBash,self).__str__())

    def apply_parameters(self,scope,con):
        s=super(EmbedBash,self).apply_parameters(scope,con)
        s.__template=self.__template.rescope({self:s, scope:s})
        s.__template.defscopes=[self]+self.defscopes
        return s

    def is_valid_rvalue(self,con):
        return self.__template is not None

    def string_context(self,con): 
        return '%d'%(self.numeric_context(con),)

    def settemplate(self,template):
        assert(isinstance(template,String))
        self.__template=template
        self.__template.defscopes=[self]+self.defscopes

    def gettemplate(self):
        return self.__template

    def numeric_context(self,con):
        return self.run(con)

    def bash_context(self,con):
        template=self.gettemplate()
        template.defscopes=[self]+self.defscopes # workaround for unknown bug
        expanded=template.string_context(con)
        #template=template.string_context(con)
        #expanded=self.expand_string(template,con)

        stream=StringIO.StringIO()
        env=dict()
        unset_me=list()
        for k,v in self.iterlocal():
            if '%' in k or '.' in k:
                stream.write('# $%s: skip; invalid shell variable name\n'
                             %(k,))
            elif k[0:2] == '__':
                stream.write("# $%s: skip; name begins with __\n"%(k,))
            elif v is null_value:
                stream.write('# $%s: skip; has no value\n'%(k,))
            elif not v.is_scalar:
                stream.write('# $%s: skip; value is not scalar\n'%(k,))
            else:
                unset_me.append(k)
                stream.write('%s=%s\n'%(k,v.bash_context(con)))
        stream.write("# Embedded bash script:\n")
        stream.write(expanded)
        stream.write('\n# End of embedded bash script.\n')
        for k in unset_me:
            stream.write('unset %s\n'%(k,))
        if con.verbose:
            stream.write('set -xe\n\n')
        else:
            stream.write('set -e\n\n')
        script=stream.getvalue()
        stream.close()
        return script

    def run(self,con):
        script=self.bash_context(con)
        # yell('%-7s %-7s %s\n'%("RUN","BASH",script))
        cmd=produtil.run.exe("bash")
        if con.verbose:
            cmd=cmd<<'set -xue\n'+script
        else:
            cmd=cmd<<'set -ue\n'+script
        env=dict(self.iterlocal())
        if env: cmd.env(**env)
        return produtil.run.run(cmd)
            
    def logical_context(self,con): 
        return bool(self.numeric_context(con)==0)

    def new_empty(self):
        s=EmbedBash(self.defscopes)
        s.__template=self.__template.rescope({self:s})
        return s

########################################################################

class Task(Scope):
    def __init__(self,defscopes,name,runvar='run'):
        super(Task,self).__init__(defscopes)
        self.__deps=list()
        self.__name=str(name)
        self.runvar=runvar
    @property
    def name(self):
        return self.__name
    def bash_context(self,con):
        assert(self.haslocal(self.runvar))
        return self.getlocal(self.runvar).bash_context(con)

    def iterdeps(self):
        for dep in self.__deps:
            yield dep

    def add_dependency(self,dep):
        self.__deps.append(dep)

    def is_valid_rvalue(self,con):
        return self.haslocal(self.runvar) and \
            self.getlocal(self.runvar) is not null_value

    def string_context(self,con): 
        return self.getlocal(self.runvar).string_context(con)

    def numeric_context(self,scopes,con):
        return self.getlocal(self.runvar).numeric_context(con)

    def run(self,con):
        return self.getlocal(self.runvar).run(con)
            
    def logical_context(self,con): 
        return self.getlocal(self.runvar).numeric_context(con)

    def new_empty(self):
        return Task(self.defscopes,self.__name,self.runvar)

########################################################################

class Build(Task):
    def __init__(self,defscopes,name):
        super(Build,self).__init__(defscopes,name,'build')
    def new_empty(self):
        return Build(self.defscopes,self.name)

########################################################################

class Platform(Task):
    def __init__(self,defscopes,name):
        super(Platform,self).__init__(defscopes,name,'detect')
    def new_empty(self):
        return Platform(self.defscopes,self.name)

########################################################################

class Test(Scope):
    def __init__(self,defscopes,name,mode):
        super(Test,self).__init__(defscopes)
        assert(mode in [ BASELINE, EXECUTION ])
        self.mode=mode
        self.__name=str(name)
        self.__deps=list()
        self.is_test=True
    @property
    def name(self):
        return self.__name

    def bash_context(self,con):
        if self.mode==BASELINE:
            steps=['prep','input','execute','make_baseline']
        else:
            steps=['prep','input','execute','verify']

        name=self.resolve("TEST_NAME").bash_context(con)
        try:
            descr=self.resolve("TEST_DESCR").bash_context(con)
        except KeyError as ke:
            descr='no description'
        report=self.resolve("COM").bash_context(con)
        report=os.path.join(report,'report.txt')

        out=StringIO.StringIO()
        out.write("report_start %s Test %s starting at $( date ) '('%s')'\n"%(
                report,name,descr))
        for step in steps:
            try:
                stepobj=self.getlocal(step)
            except KeyError as ke:
                if step in [ 'make_baseline', 'verify' ]:
                    stepobj=self.getlocal('output')
                else:
                    raise
            out.write(stepobj.bash_context(con))
            out.write('\n')
        out.write("report_finish\n")
        ret=out.getvalue()
        out.close()
        return ret

    def iterdeps(self):
        for dep in self.__deps:
            yield dep

    def add_dependency(self,dep):
        self.__deps.append(dep)

    def is_valid_rvalue(self,con):
        if self.mode==BASELINE:
            steps=['prep','input','execute','make_baseline']
        else:
            steps=['prep','input','execute','verify']

        for step in steps:
            if self.haslocal(step):
                if self.getlocal(step) is not null_value:
                    continue
            elif step in ['make_baseline','verify'] and \
                    self.haslocal('output'):
                if self.getlocal('output') is not null_value:
                    continue
            raise KeyError(step)
        return True

    def string_context(self,con): 
        raise TypeError('A Test cannot be evaluated in a string context.')

    def numeric_context(self,con):
        raise TypeError('A Test cannot be evaluated in a numeric context.')

    def run(self,con):
        raise TypeError('A Test cannot be run directly.')
            
    def logical_context(self,con): 
        raise TypeError('A Test cannot be evaluated in a logical context.')

    def new_empty(self):
        return Test(self.defscopes,self.name,self.mode)

########################################################################

class AutoDetectPlatform(object):
    def __init__(self):
        super(AutoDetectPlatform,self).__init__()
        self.__platforms=list()
    def add(self,platform):
        self.__platforms.append(platform)
    def detect(self,con):
        matches=list()
        names=list()
        for platform in self.__platforms:
            detecter=platform.resolve('detect')
            name=platform.resolve('PLATFORM_NAME')
            name=name.string_context(con)
            con.info('%s: detection...'%(name,))
            detection=detecter.logical_context(con)
            if detection:
                con.info('%s: PLATFORM DETECTED'%(name,))
                matches.append(platform)
                names.append(name)
            else:
                con.info('%s: not detected'%(name,))
        con.info('List of platforms detected: '+
                 ' '.join([ repr(s) for s in names ]))
        return matches

########################################################################

class Numeric(BaseObject):
    def __init__(self,value):
        super(Numeric,self).__init__([])
        self.is_scalar=True
        self.__value=value
    def string_context(self,con):
        return '%g'%self.__value
    def bash_context(self,con):
        return '"%g"'%self.__value
    def numeric_context(self,con):
        return self.__value
    def logical_context(self,con):
        return self.numeric_context(con)!=0
    def __str__(self):
        return str(self.__value)
    def __repr__(self):
        return repr(self.__value)
    def rescope(self,scopemap=None,prepend=None): 
        return Numeric(self.__value)
    
########################################################################

class String(BaseObject):
    def __init__(self,defscopes,value,should_expand):
        super(String,self).__init__(defscopes)
        self.__value=str(value)
        self.is_scalar=True
        self.should_expand=bool(should_expand)
    def rescope(self,scopemap=None,prepend=None):
        if prepend is None: prepend=[]
        if scopemap is None: scopemap={}
        assert(isinstance(scopemap,dict))
        assert(isinstance(prepend,list))
        defscopes=list(prepend)
        for s in self.defscopes:
            assert(isinstance(s,Scope))
            if s in scopemap:
                defscopes.append(scopemap[s])
            else:
                defscopes.append(s)
                # prepend+[ scopemap[s] if s in scopemap else s
                #         for s in self.defscopes ],
        return String(defscopes,
                      self.__value,self.should_expand)
    def string_context(self,con):
        if self.should_expand:
            return self.defscopes[0].expand_string(
                self.__value,con,self.defscopes[1:])
        else:
            return self.__value
    def bash_context(self,con):
        return bashify_string(self.string_context(con))
    def logical_context(self,con):
        s=self.string_context(con)
        s=s[-30:].lower()
        if s in [ '.true.', 'true', 'yes', 't', 'y' ]: return True
        if s in [ '.false.', 'false', 'no', 'f', 'n' ]: return False
        try:
            i=float(s)
        except ValueError as ve:
            pass
        raise ValueError('Cannot parse %s as a logical value.'%(s,))
    def numeric_context(self,con):
        s=self.string_context(con)
        return float(s)
    def __str__(self): return self.__value
    def __repr__(self): return 'String(%s)'%(repr(self.__value),)

########################################################################

class Environ(Scope):
    def __init__(self):
        super(Environ,self).__init__([])
        self.can_be_used=False
    def new_empty(self): return Environ()
    def bash_context(self,con):
        raise Exception('Cannot express the environment in a bash context.')
    def string_context(self,con):
        raise Exception('Cannot evaluate the environment in a string context.')
    def no_nulls(self): return True
    def _set_parameters(self,update):
        raise Exception('Cannot set parameters in the environment.')
    def numeric_context(self,con):
        raise Exception('Cannot evaluate the environment in a numeric context.')
    def logical_context(self,con):
        raise Exception('Cannot evaluate the environment in a logical context.')
    def as_parameters(self,con):
        raise Exception('Cannot turn the environment into a parameter list.')
    def rescope(self,scopemap=None,prepend=None):
        return self
    def has_parameters(self):
        return False
    def use_from(self,used_scope,only_scalars=False):
        raise Exception('Cannot use other scopes within the environment.')
    def apply_parameters(self,scope,con):
        raise Exception('Cannot call the environment.')
    def __str__(self):
        return 'Environ()'
    def __repr__(self):
        return 'Environ()'
    def subscope(self,key):
        raise TypeError('The environment has no subscopes.')
    def getlocal(self,key):
        return String([self],os.environ[key],False)
    def setlocal(self,key,value):
        raise Exception('Refusing to modify the environment.')
    def haslocal(self,key):
        return key in os.environ
    def iterlocal(self):
        for k,v in os.environ.iteritems():
            yield k,String([self],v,False)
    def resolve(self,key,scopes=None):
        if '.' in key or '%' in key:
            raise ValueError('Invalid environment variable \"%s\".'%(key,))
        return os.environ[key]
    def force_define(self,key,value):
        raise Exception('Refusing to modify the environment.')
    def check_define(self,key,value):
        raise Exception('Refusing to modify the environment.')
