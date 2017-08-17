import sys, re, StringIO, collections, os, datetime, logging, math
import produtil.run, produtil.log, produtil.setup

# This module really does use everything public from utilities,
# parsetree and tokenize, hence the "import *"
from produtil.testing.utilities import *
from produtil.testing.parsetree import *
from produtil.testing.tokenize import *

from produtil.testing.setarith import arithparse

__all__=[ 'Parser' ]

########################################################################

class RunConPair(object):
    def __init__(self,runnable,context):
        self.runnable=runnable
        self.context=context
    @property
    def as_tuple(self):
        """!A tuple self.runnable,self.context"""
        return self.runnable,self.context
    def __iter__(self):
        yield self.runnable
        yield self.context
    def __hash__(self):
        return hash(self.runnable)
    def __eq__(self,other):
        if not isinstance(other,RunConPair):
            return self.runnable==other
        return self.runnable==other.runnable
    def __cmp__(self,other):
        if not isinstance(other,RunConPair):
            return cmp(self.runnable,other)
        return cmp(self.runnable,other.runnable)

class Parser(object):
    def __init__(self,run_mode=None,logger=None,verbose=True):
        super(Parser,self).__init__()
        if logger is None: logger=module_logger
        if run_mode is None: run_mode=EXECUTION
        if run_mode is not EXECUTION and run_mode is not BASELINE:
            raise ValueError(
                'The Parser.__init__ run_mode argument must be the '
                'special module constants EXECUTION or BASELINE.')
        self.__runsets=collections.defaultdict(ListableSet)
        self.__runobjs=dict()
        self.__run_mode=run_mode
        self.__logger=logger
        self.__verbose=bool(verbose)
    @property
    def run_mode(self):
        return self.__run_mode
    @property
    def verbose(self):
        return self.__verbose
    @property
    def logger(self):
        return self.__logger
    @property
    def allset(self):
        """!Returns the special "**all**" runset, which contains all
        runnables that had an explicit run statement."""
        return self.__runsets['**all**']
    def setarith(self,expr=None):
        """!Executes the specified expression via the setarith module,
        returning the resulting runset as a ListableSet.
        Automatically resolves dependencies."""
        if not expr:
            runme=ListableSet(self.allset)
        else:
            runme=arithparse(expr,self.__runsets,self.__runobjs)

        result=ListableSet()
        processed=set()
        for runcon in runme:
            self._resolve_deps_impl(result,processed,runcon)
        return result
    def iterrun(self,runset='**all**'):
        """!Iterates over all runnables in the runset.  

        @param runset The runset of interest.  If no runset is
        specified then the special "**all**" runset is used, which
        contains all runnables that had a "run" command, plus
        dependencies."""
        runset=runset or '**all**'
        for runcon in self.__runsets[runset]:
            yield runcon
    def itersets(self):
        """"!Iterates over all runset names in self, including the
        special "**all**" set, which contains all runnables that had a
        "run" command, plus dependencies."""
        for setname,runset in self.__runsets.iteritems():
            yield setname,runset
    def con(self,token=None,scopes=None):
        if scopes is None: scopes=[]
        if token is None:
            return produtil.testing.parsetree.fileless_context(
                scopes=scopes,verbose=self.verbose)
        return Context(scopes,token,self.__run_mode,self.__logger,
                       verbose=self.verbose)
    def add_run(self,runset,runme,con):
        if not isinstance(runme,Task) and not isinstance(runme,Test):
            raise TypeError('The runme argument to add_run must be a Test or Task (in produtil.testing.parsetree) or subclass.  You provided a %s %s'%(type(runme).__name__,elipses(repr(runme))))
        if runset is None:
            self.__runobjs[runme.name]=RunConPair(runme,con)
            return
        addme=RunConPair(runme,con)
        for xrunset in [ runset, '**all**' ]:
            self.__runsets[xrunset].add(addme)
    def _resolve_deps_impl(self,newset,processed,runcon):
        """!Adds a runnable object and its dependencies to a list of
        objects to run.

        Adds runcon to the processed set to mark that it is being
        processed, returning immediately if it was already there.
        Then recurses into each dependency of runcon, calling this
        function on each, and then adds runcon to newset.  The result
        is that the newset will contain all runnable objects in
        correct dependency order.
        
        @param newset The ListableSet to receive the runnable objects
        @param processed a set of all runnable objects already
        processed or being processed by recursive calls to this
        function. 
        @param runcon The runnable object to process, a RunConPair"""
        assert(isinstance(runcon,RunConPair))
        assert(isinstance(newset,ListableSet))
        if runcon in processed: return
        processed.add(runcon)
        for prereq in runcon.runnable.iterdeps():
            if prereq not in newset:
                self._resolve_deps_impl(
                    newset,processed,RunConPair(prereq,runcon.context))
        newset.add(runcon)
    def resolve_deps(self):
        """!Resolves dependencies in the runsets, replacing them with
        new runsets that contain all dependencies and have runnables
        listed in correct dependency order."""
        newrunsets=dict()
        for setname in self.__runsets.iterkeys():
            processed=set()
            newset=ListableSet()
            for runcon in self.__runsets[setname]:
                assert(isinstance(runcon,RunConPair))
                self._resolve_deps_impl(newset,processed,runcon)
            newrunsets[setname]=newset
        self.__runsets=newrunsets
    def parse(self,tokenizer,scope=None,unique_id=None,morevars=None):
        assert(scope is not None)
        if scope is None:
            scope=Scope()
        if unique_id is None:
            unique_id=os.getpid()
        if morevars is not None:
            for k,v in morevars.iteritems():
                scope.setlocal(str(k),String([scope],str(v),False))
        if not isinstance(unique_id,int):
            raise TypeError(
                'The unique_id argument to Parser.parse() must be an '
                'int, not a %s %s.'%( type(unique_id).__name__,
                                      elipses(repr(unique_id)) ))
        tokiter=peekable(tokenizer)
        scope.setlocal('ENV',Environ())
        assert(isinstance(unique_id,int))
        scope.setlocal('UNIQUE_ID',Numeric(unique_id))
        try:
            result=self.parse_subscope(
                tokiter,[scope],[end_of_text_type],
                self.parse_between_assignments,
                allow_overwrite=False,
                allow_resolve=True,
                allow_run=True,
                allow_null=False,
                allow_use=False,
                allow_load=True,
                scope_name='global scope')
            self.resolve_deps()
        except Exception as e: # FIXME: change exception type
            # try:
            #     peek=tokiter.peek()
            #     filename=peek.filename
            #     lineno=peek.lineno
            #     sys.stderr.write('%s:%d: uncaught exception: %s\n'%(
            #             filename,lineno,str(e)))
            # except StopIteration as se:
            #     sys.stderr.write('StopIteration while peeking: %s\n'%(
            #             str(se),))
            #     pass
            raise
    def parse_between_arguments(self,tokiter,ends=None):
        if ends is None: ends=[')']
        peek=tokiter.peek()
        # yell('%-7s peek type=%s value=%s\n'%(
        #         'BETWEEN',str(peek.token_type),elipses(str(
        #                 peek.token_value))))
        while True:
            if peek.token_type == ',':
                tokiter.next()
                # yell('%-7s consume ,\n'%('BETWEEN',))
                return True
            elif peek.token_type in ends:
                # yell('%-7s stop at %s\n'%('BETWEEN',peek.token_type))
                return False
            elif peek.token_type==end_of_line_type:
                # yell('%-7s saw \n'%('BETWEEN',))
                tokiter.next()
                peek=tokiter.peek()
            else:
                return False
    def skip_eoln(self,tokiter):
        peek=tokiter.peek()
        while peek.token_type==end_of_line_type:
            tokiter.next()
            peek=tokiter.peek()
    def parse_between_assignments(self,tokiter):
        peek=tokiter.peek()
        seen=False
        # yell('%-7s peek type=%s value=%s in parse_between_assignments\n'%(
        #         'BETWEEN',str(peek.token_type),
        #         elipses(str(peek.token_value))))
        while True:
            seen=True
            if peek.token_type in [ ',', ';' ]:
                tokiter.next()
                # yell('%-7s consume %s\n'%('BETWEEN',peek.token_type))
                return True
            elif peek.token_type in [ '}', end_of_text_type ]:
                # yell('%-7s stop at %s\n'%('BETWEEN',peek.token_type))
                return True
            elif peek.token_type==end_of_line_type:
                # yell('%-7s saw %s\n'%('BETWEEN',repr(peek.token_type)))
                seen=True
                tokiter.next()
                peek=tokiter.peek()
            else:
                break
        return seen
        self.error('between_assignments',peek)

    def parse_embed_script(self,tokiter,scopes,ends,parse_between=None):
        token=tokiter.next()
        if token.token_type != 'varname':
            self.error('embed',token)
        if token.token_value != 'bash':
            self.error('embed',token,'unknown language "%s"'%(
                    token.token_value,))
        nametoken=tokiter.next()
        if token.token_type != 'varname':
            self.error('embed script name',token)
        scope=EmbedBash(scopes)
        token=tokiter.next()

        while token.token_type==end_of_line_type: token=tokiter.next()
        if token.token_type=='(':
            self.parse_subscope(tokiter,[scope]+scopes,[')'],
                                self.parse_between_arguments,
                                allow_overwrite=False,
                                allow_resolve=False,
                                allow_null=True,
                                only_scalars=True,
                                scope_name='embed script parameters')
            scope=scope.as_parameters(self.con(token,scopes))
            token=tokiter.next()
        while token.token_type==end_of_line_type: token=tokiter.next()

        if token.token_type=='{':
            self.parse_subscope(tokiter,[scope]+scopes,['}'],
                                self.parse_between_assignments,
                                allow_overwrite=True,
                                allow_resolve=True,
                                allow_null=False,
                                allow_use=True,
                                only_scalars=True,
                                scope_name='embed script variables')
            token=tokiter.next()
        while token.token_type==end_of_line_type: token=tokiter.next()

        if token.token_type in [ 'qstring', 'dqstring', 'bracestring' ]:
            scope.settemplate(self.action_string([scope]+scopes,token))
        else:
            self.error('embed script contents',token)
        if parse_between: 
            parse_between(tokiter)
        return (nametoken.token_value,scope)

    def parse_set_list(self,tokiter,scopes):
        peek=tokiter.peek()
        while peek.token_type=='varname':
            yield peek.token_value
            lvaltoken=peek
            tokiter.next() # discard varname
            peek=tokiter.peek()
            while peek.token_type==end_of_line_type:
                tokiter.next()
                peek=tokiter.peek()
            if peek.token_type=='==':
                # this is an "if var==value" condition
                tokiter.next() # consume ==
                peek=tokiter.peek()
                if peek.token_type!='varname':
                    self.error('run setname @ ... var==',peek)
                rvaltoken=tokiter.next()
                lval=self.action_resolve(lvaltoken,scopes)
                rval=self.action_resolve(rvaltoken,scopes)
                if lval != rval:
                    yield False

                peek=tokiter.peek()
            if peek.token_type!=',':
                return # reached end of list.
            tokiter.next() # discard ","
            peek=tokiter.peek()
                
    def parse_deplist(self,tokiter,scopes,task,ends):
        allscopes=[task]+scopes
        peek=tokiter.peek()
        while not peek.token_type in ends:
            if peek.token_type=='varname':
                varname=peek
                tokiter.next()
                peek=tokiter.peek()
                if peek.token_type==end_of_line_type:
                    continue # ignore blank lines
                elif self.parse_between_arguments(tokiter,['{']):
                    # varname is followed by a comma
                    dep=self.action_resolve(varname,allscopes)
                    self.action_dependency(task,scopes,dep)
                    peek=tokiter.peek()
                    continue
                elif peek.token_type in ends:
                    dep=self.action_resolve(varname,allscopes)
                    self.action_dependency(task,scopes,dep)
                    return
                elif peek.token_type == '(':
                    # This is a function call.
                    tokiter.next()
                    subscope=Scope(scopes)
                    self.parse_subscope(tokiter,scopes,[')'],
                                        self.parse_between_arguments,
                                        allow_overwrite=False,
                                        allow_resolve=False,
                                        allow_null=True,
                                        scope_name='dependency argument list')
                    subscope=subscope.as_parameters(self.con(peek,scopes))
                    peek=tokiter.peek()
                    if self.parse_between_arguments(tokiter) \
                            or peek.token_type in ends:
                        dep=self.action_call(varname,peek,scopes,subscope)
                        self.action_dependency(task,scopes,dep)
                        peek=tokiter.peek()
                        continue
            self.error('dependency argument list',peek)
    def parse_op_list(self,tokiter,scopes,subscope):
        token=tokiter.next()
        strings=[ 'qstring', 'dqstring', 'bracestring' ]
        if token.token_type != '{':
            self.error('operation list',token)
        while True:
            # Get target of operation:
            token=tokiter.next()
            while token.token_type==end_of_line_type:
                token=tokiter.next()
 
            if token.token_type=='varname' and token.token_value=='use':
                peek=tokiter.peek()
                if peek.token_type!='varname':
                    self.error('operation list use statement',peek)
                tokiter.next()
                subscope.use_from(self.action_resolve(peek,scopes))
                self.parse_between_assignments(tokiter)
                peek=tokiter.peek()
                if peek.token_type=='}':
                    tokiter.next()
                    return subscope
                continue
            elif token.token_type=='}':
                return subscope
            elif token.token_type not in strings:
                self.error('operator target',token)
            tgt=self.action_string(scopes,token)

            # Get operator:
            token=tokiter.next()
            while token.token_type==end_of_line_type:
                token=tokiter.next()
            if token.token_type!='oper':
                self.error('oper',token)
            op=self.action_operator(scopes,token)

            # Get source of operation (input or baseline file)
            token=tokiter.next()
            while token.token_type==end_of_line_type:
                token=tokiter.next()
            if token.token_type not in strings:
                self.error('operator source (input or baseline)',token)
            src=self.action_string(scopes,token)

            # Add operator:
            # FIXME: CONTEXT
            subscope.add_binary_operator(
                tgt,op,src,self.con(token,scopes))

            self.parse_between_assignments(tokiter)
            peek=tokiter.peek()
            if peek.token_type=='}':
                tokiter.next()
                return subscope


    def parse_hash_define(self,tokiter,scopes,subscope,parse_between=None,
                          allow_deps=False,hash_type='hash'):
        token=tokiter.next()
        parameters=False

        # if token.token_type=='(':
        #     parameters=True
        #     self.parse_subscope(tokiter,[subscope]+scopes,[')'],
        #                         self.parse_between_arguments,
        #                         allow_overwrite=False,
        #                         allow_resolve=False,
        #                         allow_null=True,
        #                         scope_name=hash_type+' argument list')
        #     subscope=subscope.as_parameters()
        #     token=tokiter.next()

        if allow_deps and token.token_type==':':
            self.parse_deplist(
                tokiter,[subscope]+scopes,subscope,['{'])
            token=tokiter.next()

        if token.token_type=='{':
            tokiter.next()
            self.parse_subscope(tokiter,[subscope]+scopes,['}'],
                                self.parse_between_assignments,
                                allow_overwrite=True,
                                allow_resolve=True,
                                allow_null=False,
                                allow_use=True,
                                scope_name=hash_type+' definition')
        else:
            self.error(
                hash_type+' definition',token,
                'missing {...} block in '+hash_type+' definition')

        if parse_between:
            parse_between(tokiter)

        # yell('%-7s define Scope@%s with %sparameters\n'%(
        #         hash_type.upper(),str(id(subscope)),
        #         ' ' if parameters else 'no '))
        return subscope

    def parse_spawn_element(self,tokiter,scopes,spawn,ends):
        token=tokiter.next()
        args=list()
        opts=list()
        saw_vars=False
        strings=[ 'qstring', 'dqstring', 'bracestring' ]
        while token.token_type not in ends:
            if token.token_type in strings:
                if saw_vars:
                    self.error('spawn process',token,'var=value elements '
                               'must come after all arguments')
                args.append(token)
                self.parse_between_arguments(tokiter,ends)
                token=tokiter.next()
                continue
            elif token.token_type==end_of_line_type:
                token=tokiter.next()
                continue
            elif token.token_type=='varname':
                name=token.token_value
                peek=tokiter.peek()
                if peek.token_type != '=':
                    self.error('spawn process',token)
                tokiter.next()
            else:
                self.error('spawn process',token)
            # we're at the value in varname=value
            rvalue=self.parse_rvalue(
                tokiter,scopes,['}',','],
                only_scalars=True)
            self.parse_between_arguments(tokiter,ends)
            opts.append([name,rvalue])
            token=tokiter.next()
        scope=Scope(scopes)
        allscopes=[scope]+scopes
        for k,v in opts:
            scope.setlocal(k,v)
        if not args:
            self.error('spawn process',token,'no command nor arguments')
        for arg in args:
            assert(isinstance(arg,Token))
        argobjs=[ 
            self.action_string(allscopes,arg) for arg in args]
        return argobjs,scope

    def parse_spawn_block(self,tokiter,scopes,name,spawn,ends,parse_between):
        token=tokiter.next()
        while token.token_type==end_of_line_type: 
            token=tokiter.next()
        while token.token_type not in ends:
            if token.token_type!='{':
                self.error('spawned process',token)
            (args,opts)=self.parse_spawn_element(
                    tokiter,scopes,spawn,['}'])
            spawn.add_rank(args,opts)
            if parse_between: parse_between(tokiter)
            token=tokiter.next()

    def parse_spawn(self,tokiter,scopes,name,spawn):
        token=tokiter.next()
        if token.token_type!='{':
            self.error('spawn block',token)
        while self.parse_spawn_block(tokiter,scopes,name,spawn,['}'],
                                     self.parse_between_assignments):
            continue
        return spawn

    def parse_autodetect(self,tokiter,scopes,taskname,task):
        # Check for the (/ and skip it:
        token=tokiter.next()
        while token.token_type==end_of_line_type:
            token=tokiter.next()
        if token.token_type!='(/':
            self.error('autodetect platform list',token)

        while True:
            peek=tokiter.peek()
            while peek.token_type==end_of_line_type:
                tokiter.next()
                peek=tokiter.peek()
            if peek.token_type=='/)':
                tokiter.next()
                return
            rvalue=self.parse_rvalue(tokiter,scopes,['/)'],
                                     self.parse_between_arguments,False)
            task.add(rvalue)
            peek=tokiter.peek()
            self.parse_between_arguments(tokiter)
            if peek.token_type=='/)':
                tokiter.next()
                return

    def parse_load(self,tokiter,scope,seen_run):
        filetoken=tokiter.next()
        if filetoken.token_type!='qstring':
            self.error('load',token,"load statements can only include "
                       "'single-quote strings'")
        eoln=tokiter.peek()
        if eoln.token_type not in [ end_of_line_type, end_of_text_type ]:
            self.error('load',eoln,"a load statement must be followed "
                       "by an end of line or the end of the script.")
        newfile=filetoken.token_value
        if not os.path.isabs(newfile):
            newfile=os.path.join(os.path.dirname(filetoken.filename),newfile)
        tokenizer=tokiter.child
        with open(newfile,'rt') as fileobj:
            new_tokenizer=tokenizer.for_file(fileobj,newfile)
            new_tokiter=peekable(new_tokenizer)
            self.parse_subscope(
                new_tokiter,[scope],[end_of_text_type],
                self.parse_between_assignments,
                allow_overwrite=False,
                allow_resolve=True,
                allow_run=True,
                allow_null=False,
                allow_use=False,
                allow_load=True,
                scope_name='global scope',
                seen_run=seen_run)
        
    def parse_subscope(self,tokiter,scopes,ends,parse_between,
                       allow_overwrite=True,allow_resolve=True,
                       allow_run=False,allow_null=False,
                       allow_use=False,scope_name='subscope',
                       only_scalars=False,allow_load=False,
                       seen_run=False):
        go=True # set to False once an "ends" is seen
        seen_run=bool(seen_run) # Did we see an execution request yet?
        token=None
        strings=[ 'qstring', 'dqstring', 'bracestring' ]
        def define(con,key,val):
            if seen_run:
                self.error(
                    scope_name,token,
                    reason='Definitions must come before execution '
                    'requests.')
            if not val.is_valid_rvalue(con):  # FIXME: con
                self.error(scope_name,token,'not a valid rvalue: %s'%(
                        elipses(repr(val)),))
            # yell('%s:%s: define %s=%s\n'%(
            #         token.filename,str(token.lineno),
            #         str(key),repr(val)))
            if allow_overwrite:
                scopes[0].force_define(key,val)
            else:
                scopes[0].check_define(key,val)
        if allow_resolve:
            search_scopes=scopes
        else:
            search_scopes=scopes[1:]
        while go:
            token=tokiter.next()
            if token.token_type=='varname':
                peek=tokiter.peek()
                if peek.token_type=='=':
                    tokiter.next()
                    define(self.con(token,scopes),token.token_value,
                           self.parse_rvalue(tokiter,search_scopes,ends,
                                             parse_between,
                                             only_scalars=only_scalars))
                    parse_between(tokiter)
                    continue
                elif token.token_value=='load' and peek.token_type in strings:
                    if not allow_load:
                        self.error('subscope',token,'load statements are '
                                   'only allowed in the global scope.')
                    self.parse_load(tokiter,scopes[-1],seen_run)
                    if parse_between:  parse_between(tokiter)
                    continue
                elif token.token_value=='use' and peek.token_type=='varname' \
                        and allow_use:
                    tokiter.next() # consume the peeked value
                    self.action_use(scopes,peek,
                                    only_scalars=only_scalars)
                    if parse_between:  parse_between(tokiter)
                    continue
                elif allow_run and token.token_value=='run' \
                        and peek.token_type=='varname':
                    set_rvalue=self.parse_rvalue(tokiter,search_scopes,
                                                 ends,parse_between)
                    set_con=self.con(peek,scopes)
                    runsets=list()
                    peek=tokiter.peek()
                    keep=True
                    if peek.token_type=='@':
                        tokiter.next() # discard the @
                        for setname in self.parse_set_list(
                                tokiter,search_scopes):
                            if setname is False:
                                keep=False
                            elif setname not in runsets:
                                assert(isinstance(setname,basestring))
                                runsets.append(setname)
                    if keep:
                        if not runsets:
                            runsets.append('**all**')
                        seen_run=True
                        self.action_run_by_name(set_rvalue,set_con)
                        for setname in runsets:
                            runobj=self.action_run_in_set(
                                setname,set_rvalue,set_con)
                    peek=tokiter.peek()
                    if parse_between: parse_between(tokiter)
                    continue
                elif not only_scalars and token.token_value=='spawn' \
                        and peek.token_type=='varname':
                    taskname=peek.token_value
                    tokiter.next()
                    task=self.parse_spawn(tokiter,scopes,peek.token_value,
                                          SpawnProcess(scopes))
                    define(self.con(peek,scopes),taskname,task)
                    if parse_between:  parse_between(tokiter)
                    del taskname,task
                    continue
                elif not only_scalars and token.token_value in [
                    'filters', 'criteria' ] and peek.token_type=='varname':
                    taskname=peek.token_value
                    tokiter.next()
                    if token.token_value=='filters':
                        task=Filters(scopes)
                    elif token.token_value=='criteria':
                        task=Criteria(scopes)
                    task=self.parse_op_list(tokiter,scopes,task)
                    define(self.con(peek,scopes),taskname,task)
                    if parse_between:  parse_between(tokiter)
                    del taskname,task
                    continue
                elif not only_scalars and token.token_value=='autodetect' \
                        and peek.token_type=='varname':
                    taskname=peek.token_value
                    taskcon=self.con(peek,scopes)
                    tokiter.next() # Skip name token
                    task=AutoDetectPlatform()
                    self.parse_autodetect(tokiter,scopes,taskname,task)
                    task=self.action_autodetect(self.con(peek,scopes),
                                                tokiter,scopes,taskname,task)
                    define(taskcon,taskname,task)
                    del taskname, task, taskcon
                    continue
                elif not only_scalars and token.token_value in [
                        'build', 'task', 'test', 'platform' ] \
                        and peek.token_type=='varname':
                    taskname=peek.token_value
                    # yell('%-7s %-7s %s\n'%(
                    #         'PARSE',token.token_value,taskname))
                    tokiter.next() # consume the task name
                    # yell('%-7s %-7s %s\n'%(
                    #         'INIT',token.token_value,taskname))
                    if token.token_value=='task':
                        raise AssertionError('Should never make a Task')
                        task=Task(scopes,taskname)
                    elif token.token_value=='test':
                        task=Test(scopes,taskname,self.__run_mode)
                        task._set_constants({'TEST_NAME':
                                      String(scopes,taskname,False)})
                    elif token.token_value=='build':
                        task=Build(scopes,taskname)
                        task._set_constants({'BUILD_NAME':
                                      String(scopes,taskname,False)})
                    elif token.token_value=='platform':
                        task=Platform(scopes,taskname)
                        task._set_constants({'PLATFORM_NAME':
                                      String(scopes,taskname,False)})
                    else:
                        raise AssertionError(
                            'Unrecognized subscope type "%s".'%(
                                token.token_value,))
                    task=self.parse_hash_define(
                            tokiter,scopes,task,parse_between,
                            allow_deps=token.token_value!='platform')
                    # yell('%-7s %-7s %s\n'%('DEFINE',token.token_value,
                    #                        taskname))
                    define(self.con(peek,scopes),taskname,task)
                    del taskname, task
                    if parse_between:  parse_between(tokiter)
                    continue
                elif not only_scalars \
                        and token.token_value=='hash' \
                        and peek.token_type=='varname':
                    hashname=peek.token_value
                    tokiter.next() # consume the hash name
                    define(self.con(peek,scopes),
                           hashname,self.parse_hash_define(
                            tokiter,scopes,Scope(scopes),parse_between))
                    del hashname
                    if parse_between:  parse_between(tokiter)
                    continue
                elif not only_scalars \
                        and token.token_value=='embed' \
                        and peek.token_type=='varname':
                    (varname,script)=self.parse_embed_script(
                        tokiter,scopes,parse_between)
                    define(self.con(peek,scopes),varname,script)
                    if parse_between:  parse_between(tokiter)
                    continue
                elif allow_null and (
                    peek.token_value in ends or
                    parse_between and parse_between(tokiter)):
                    define(self.con(peek,scopes),
                           token.token_value,null_value)
                    if parse_between:  parse_between(tokiter)
                    continue
            elif token.token_type in ends:
                return scopes[0]
            elif token.token_type==end_of_line_type:
                continue # ignore blank lines.
            self.error(scope_name,token)
    def parse_rvalue(self,tokiter,scopes,ends,parse_between=None,
                     only_scalars=False):
        token=tokiter.next()
        if token.token_type in [ 'qstring', 'dqstring', 'bracestring' ]:
            ret=self.action_string(scopes,token)
            if parse_between: parse_between(tokiter)
            return ret
        elif token.token_type == 'number':
            ret=self.action_numeric(scopes,token)
            if parse_between: parse_between(tokiter)
            return ret
        elif not only_scalars and token.token_type in '{':
            subscope=Scope(scopes)
            ret=self.parse_subscope(
                tokiter,[subscope]+scopes,['}'],
                self.parse_between_assignments,
                allow_overwrite=True,
                allow_resolve=True,
                allow_run=False,
                allow_null=False,
                allow_use=True,
                scope_name="hash")
            if parse_between: parse_between(tokiter)
            return ret
        elif not only_scalars and token.token_type=='varname':
            varname=token.token_value
            peek=tokiter.peek()
            if peek.token_type=='(':
                # We are at the ( in varname(arguments...
                tokiter.next() # consume (
                subscope=Scope(scopes)
                scopesplus=[subscope]+scopes
                self.parse_subscope(tokiter,scopesplus,[')'],
                                    self.parse_between_arguments,
                                    allow_overwrite=False,
                                    allow_resolve=False,
                                    allow_null=True,
                                    scope_name='argument list')
                peek=tokiter.peek()
                if peek.token_type in ends or \
                        parse_between and parse_between(tokiter):
                    # This is a function call varname(arg,arg,...)
                    return self.action_call(varname,peek,scopes,subscope)
            elif peek.token_type in ends :
                return self.action_resolve(token,scopes)
            elif parse_between and parse_between(tokiter):
                return self.action_resolve(token,scopes)
        self.error('rvalue',token)
    def action_autodetect(self,con,tokiter,scopes,taskname,task):
        matches=task.detect(con)
        if len(matches)==0:
            raise Exception(
                'You are using an unknown platform.  Fixme: need better '
                'exception here.')
        elif len(matches)>1:
            raise Exception(
                'This machine can submit to multiple platforms: '+(
                    ' '.join([
                            s.resolve('PLATFORM_NAME') \
                                .string_context(fileless_context(
                                                    verbose=self.verbose)) \
                            for s in matches
                ])))
        return matches[0]
    def action_dependency(self,task,scopes,dep):
        task.add_dependency(dep)
    def action_call(self,varname,token,scopes,parameters):
        # yell('%-7s %s in parameter scope %s\n'%(
        #         'CALL',repr(varname),repr(parameters)))
        callme=scopes[0].resolve(varname)
        # yell('CALL APPLY PARAMETERS\n')
        return callme.apply_parameters(parameters,self.con(token,scopes))
    def action_use(self,scopes,key_token,only_scalars=False):
        assert(isinstance(key_token,Token))
        assert(isinstance(scopes,list))
        assert(len(scopes)>=2)
        assert(isinstance(scopes[0],Scope))
        key=key_token.token_value
        got=scopes[1].resolve(key)
        found_non_scalars=scopes[0].use_from(got,only_scalars)
        if only_scalars and found_non_scalars:
            self.error('use',key_token,'found non-scalars when '
                       'using %s'%(key,))
        # for k,v in got.iterlocal():
        #     if only_scalars and not isinstance(v,String):
        #         self.error('use',key_token,'found non-scalars when '
        #                    'using %s'%(key,))
        #     scopes[0].setlocal(k,v)
    def action_operator(self,scopes,token):
        assert(isinstance(token,Token))
        if token.token_value=='.copy.':
            return Copy(scopes)
        elif token.token_value=='.copydir.':
            return CopyDir(scopes)
        elif token.token_value=='.bitcmp.':
            return BitCmp(scopes)
        elif token.token_value=='.md5cmp.':
            return Md5Cmp(scopes)
        elif token.token_value=='.link.':
            return Link(scopes)
        elif token.token_value=='.atparse.':
            return AtParse(scopes)
        else:
            self.error('operator name',token,'unknown operator '+
                       token.token_value)
    def action_numeric(self,scopes,token):
        value=float(token.token_value)
        return Numeric(value)
    def action_string(self,scopes,token):
        assert(isinstance(token,Token))
        if token.token_type=='qstring':
            s=String(scopes,token.token_value,False)
        elif token.token_type=='dqstring':
            s=String(scopes,dqstring2bracestring(token.token_value),True)
        elif token.token_type=='bracestring':
            s=String(scopes,token.token_value,True)
        else:
            raise ValueError('Invalid token for a string: '+repr(token))
        # yell('%-7s %s = %s\n'%('STRING',repr(token.token_value),repr(s)))
        return s
    def action_resolve(self,varname_token,scopes):
        varname=varname_token.token_value
        for scope in scopes:
            try:
                return scope.resolve(varname)
            except KeyError as ke:
                pass
        raise KeyError(varname)
    def action_null_param(self,varname,scope):
        if '%' in varname:
            raise ValueError('%s: cannot have "%" in a parameter name.'%(
                    varname,))
        scope.check_define(varname,null_value)
    def action_assign_var(self,toscope,tovar,fromvar,fromscopes,
                          allow_overwrite):
        if fromscopes:
            value=fromscopes[0].resolve(fromvar)
        else: # Global scope assignment
            value=toscope.resolve(fromvar)
        self.action_assign(toscope,tovar,value,allow_overwrite)
    def action_assign(self,scope,varname,value,allow_overwrite):
        assert(isinstance(scope,Scope))
        assert(isinstance(varname,basestring))
        assert(isinstance(value,BaseObject))
        if '%' in varname:
            raise ValueError('Cannot assign to %s; subscope definitions must be of syntax "var1 = { var2= { ...."'%(
                    varname,))
        # yell('%-7s %s = %s IN %s\n'%(
        #     'ASSIGN', varname, repr(value),repr(scope) ))
        if allow_overwrite:
            scope.force_define(varname,value)
        else:
            scope.check_define(varname,value)
    def action_run_in_set(self,setname,obj,con):
        # yell('%-7s %s\n'%(
        #     'RUN', repr(obj)))
        return self.add_run(setname,obj,con)
    def action_run_by_name(self,obj,con):
        return self.add_run(None,obj,con)
    def error(self,mode,token,reason=None):
        if token is None:
            raise Exception('Unexpected end of file.')
        elif reason:
            raise Exception('%s:%s: %s (%s token with value %s)'%(
                    token.filename,token.lineno,str(reason),
                    repr(token.token_type),
                    repr(elipses(str(token.token_value)))))
        else:
            raise Exception(
                '%s:%s: unexpected %s in %s (token value %s)'%(
                    token.filename, token.lineno, repr(token.token_type),
                    str(mode), repr(elipses(str(token.token_value)))))
