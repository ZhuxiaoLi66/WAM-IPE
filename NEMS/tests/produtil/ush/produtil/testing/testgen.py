"""!Utility package for making programs that are wrapped around the
produtil.testing package."""

import os, logging

import produtil.fileop

from produtil.testing.utilities import BASELINE, EXECUTION
from produtil.testing.tokenize import Tokenizer, TokenizeFile
from produtil.testing.parse import Parser
from produtil.testing.rocoto import RocotoRunner
from produtil.testing.script import BashRunner
from produtil.testing.parsetree import fileless_context
from produtil.testing.setarith import arithparse

__all__=[ 'TestGen' ]

class TestGen(object):
    def __init__(self, run_mode, OutputType, outloc, inloc, dry_run, 
                 unique_id, logger=None, verbose=True, PWD=None,
                 setarith=None):
        if PWD is None:
            # Default PWD is produtil:
            PWD=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        if logger is None:
            logger=logging.getLogger('testgen')
        self.setarith=setarith
        self.logger=logger
        self.run_mode=run_mode
        self.OutputType=OutputType
        self.outloc=outloc
        self.inloc=inloc
        self.dry_run=dry_run
        self.unique_id=unique_id
        self.verbose=bool(verbose)
        self.PWD=PWD
        self.scope=None
        self.parser=None
        self.parse_result=None
        assert(isinstance(self.unique_id,int))

    def get_string(self,varname):
        var=self.scope.resolve(varname)
        con=fileless_context(verbose=self.verbose)
        return var.string_context(con)

    def make_vars(self):
        here=self.PWD
        PWD_UP1=os.path.dirname(here)
        PWD_UP2=os.path.dirname(PWD_UP1)
        PWD_UP3=os.path.dirname(PWD_UP2)
        PWD_UP4=os.path.dirname(PWD_UP3)
        PWD_UP5=os.path.dirname(PWD_UP4)
        return { 'PWD_UP1':PWD_UP1, 'PWD_UP2':PWD_UP2,
                 'PWD_UP3':PWD_UP3, 'PWD_UP4':PWD_UP4,
                 'PWD_UP5':PWD_UP5, 'OUTPUT_PATH':self.outloc,
                 'PWD': here }
    def make_more(self,result,con):
        """!This routine is used by subclasses to make any additional
        files from within testgen().  The default implementation does
        nothing.

        @param result the produtil.testing.parsetree.Scope for the
        global scope of the parsed files

        @param con the produtil.testing.parsetree.Context to use when
        expanding strings or other objects from result."""
    def override(self,scope):
        """!Overrides variables in the scope before parsing."""
    def parse(self):
        logger=self.logger
        tokenizer=Tokenizer()
        self.scope=produtil.testing.parsetree.Scope()
        self.override(self.scope)
        self.parser=Parser(self.run_mode,logger,self.verbose)
        morevars=self.make_vars()
        with open(self.inloc,'rt') as fileobj:
            self.parse_result=self.parser.parse(
                TokenizeFile(tokenizer,fileobj,self.inloc,1),self.scope,
                unique_id=self.unique_id,morevars=morevars)
    def generate(self):
        logger=self.logger
        outputter=self.OutputType()
        outputter.make_runner(parser=self.parser,dry_run=self.dry_run,
                              setarith=self.setarith)
        con=fileless_context(
            scopes=[self.parse_result],verbose=self.verbose,logger=logger,
            run_mode=self.run_mode)
        self.make_more(self.parse_result,con)
    def testgen(self):
        self.parse()
        self.generate()
            
